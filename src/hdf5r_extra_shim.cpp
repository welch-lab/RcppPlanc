#include <functional>
#include <highfive/highfive.hpp>
#include <Rcpp.h>

namespace Rcpp {
  class SpMat {
  public:
    IntegerVector i, p, Dim;
    NumericVector x;
    SpMat(S4 rec) {
      i = rec.slot("i");
      p = rec.slot("p");
      Dim = rec.slot("Dim");
      x = rec.slot("x");
    }
  };
  template <> SpMat as(SEXP rec) {return SpMat(rec);}
}

// [[Rcpp::export(.rcpp_mat_to_h5mat)]]
Rcpp::CharacterVector RcppToH5Mat(Rcpp::NumericMatrix x, Rcpp::String filename, Rcpp::String datapath, bool overwrite, Rcpp::Nullable<Rcpp::IntegerVector> chunk_size = R_NilValue) {
  HighFive::File::AccessMode props;
  if (overwrite) {
    props = HighFive::File::Overwrite;
  }
  else {
    props = HighFive::File::ReadWrite|HighFive::File::Excl;
  }
  try {
    HighFive::File file(filename, props);
    HighFive::DataSetCreateProps cparms;
    if (!chunk_size.isNull()) {
      Rcpp::IntegerVector chunk_size_vec = chunk_size.get();
      if (chunk_size_vec.size() != 2) {
        Rcpp::stop("chunk_size must be a vector of length 2.");
      }
      std::vector<hsize_t> chunk_dims = {static_cast<hsize_t>(chunk_size_vec[1]), static_cast<hsize_t>(chunk_size_vec[0])};
      cparms.add(HighFive::Chunking(chunk_dims));
    } else {
      std::vector<hsize_t> chunk_dims = {static_cast<hsize_t>(x.ncol()), static_cast<hsize_t>(x.nrow())};
      cparms.add(HighFive::Chunking(chunk_dims));
    }
    HighFive::DataSet dataset =  file.createDataSet(datapath, HighFive::DataSpace({x.ncol(), x.nrow()}, {HighFive::DataSpace::UNLIMITED, HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<double>() ,cparms);
    dataset.write_raw(x.begin());
  }
  catch (const HighFive::FileException& e) {
    Rcpp::stop(e.what());
  }
  return Rcpp::CharacterVector::create(filename, datapath);
  }

// [[Rcpp::export(.rcpp_spmat_to_h5spmat)]]
Rcpp::CharacterVector RcppToH5Spmat(const SEXP& x, Rcpp::String filename, Rcpp::String datapath, bool overwrite) {
  Rcpp::SpMat xp = Rcpp::as<Rcpp::SpMat>(x);
  HighFive::File::AccessMode props;
  if (overwrite) {
    props = HighFive::File::Overwrite;
  }
  else {
    props = HighFive::File::ReadWrite|HighFive::File::Excl;
  }
  HighFive::DataSetCreateProps datacparms;
  std::vector<hsize_t> chunk_dimdata = {static_cast<hsize_t>(xp.x.size())};
  datacparms.add(HighFive::Chunking(chunk_dimdata));
  HighFive::DataSetCreateProps indcparms;
  std::vector<hsize_t> chunk_dimind = {static_cast<hsize_t>(xp.i.size())};
  indcparms.add(HighFive::Chunking(chunk_dimind));
  HighFive::DataSetCreateProps ptrcparms;
  std::vector<hsize_t> chunk_dimptr = {static_cast<hsize_t>(xp.p.size())};
  int dimArray[2] = {xp.Dim[1], xp.Dim[0]};
  ptrcparms.add(HighFive::Chunking(chunk_dimptr));
  try {
    HighFive::File file(filename, props);
    if (datapath == "/") {
      HighFive::DataSet data = file.createDataSet("data", HighFive::DataSpace({xp.x.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<double>(), datacparms);
      HighFive::DataSet ind = file.createDataSet("indices", HighFive::DataSpace({xp.i.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<int>(), indcparms);
      HighFive::DataSet ptr = file.createDataSet("indptr", HighFive::DataSpace({xp.p.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<int>(), ptrcparms);
      HighFive::Attribute shape = file.createAttribute("shape", HighFive::DataSpace({2}), HighFive::create_datatype<int>());
      data.write_raw(xp.x.begin());
      ind.write_raw(xp.i.begin());
      ptr.write_raw(xp.p.begin());
      shape.write(dimArray);
    }
    else {
      HighFive::Group group = file.createGroup(datapath);
      HighFive::DataSet data = group.createDataSet("data", HighFive::DataSpace({xp.x.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<double>(), datacparms);
      HighFive::DataSet ind = group.createDataSet("indices", HighFive::DataSpace({xp.i.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<int>(), indcparms);
      HighFive::DataSet ptr = group.createDataSet("indptr", HighFive::DataSpace({xp.p.size()}, {HighFive::DataSpace::UNLIMITED}), HighFive::create_datatype<int>(), ptrcparms);
      HighFive::Attribute shape = group.createAttribute("shape", HighFive::DataSpace({2}), HighFive::create_datatype<int>());
      data.write_raw(xp.x.begin());
      ind.write_raw(xp.i.begin());
      ptr.write_raw(xp.p.begin());
      shape.write(dimArray);
    }
  }
  catch (const HighFive::FileException& e) {
    Rcpp::stop(e.what());
  }
  std::string cdatapath(datapath);
  return Rcpp::CharacterVector::create(filename, cdatapath += "/data", cdatapath + "/indices", cdatapath + "/indptr");
}
