#ifndef ARMA_DONT_PRINT_FAST_MATH_WARNING
#define ARMA_DONT_PRINT_FAST_MATH_WARNING
#endif

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export(.rcpp_mat_to_h5mat)]]
CharacterVector rcpp_mat_to_h5mat(const NumericMatrix& x, std::string filename,
    std::string dataPath) {
    CharacterVector res(2);
    res[0] = filename;
    res[1] = dataPath;
    // Existance of `filename` is already done in R
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    // Create a dataset of native double with the same shape as the matrix
    // // Set the dimensions of the dataset
    std::array<size_t, 2> tDims;
    tDims[0] = x.ncol();
    tDims[1] = x.nrow();
    HighFive::DataSpace fspace(tDims);
    // // Set chunking dimension for the new dataset
    std::vector<hsize_t> chunk_dims;
    unsigned int colChunkSize = 1000, rowChunkSize = 1000;
    if (colChunkSize > x.ncol()) {
        // Mainly happening in small unit test case, but worth checking
        chunk_dims.push_back(x.ncol());
    } else {
        chunk_dims.push_back(colChunkSize);
    }
    if (rowChunkSize > x.nrow()) {
        // Mainly happening in small unit test case, but worth checking
        chunk_dims.push_back(x.nrow());
    } else {
        chunk_dims.push_back(rowChunkSize);
    }
    HighFive::Chunking Chunks(chunk_dims);
    HighFive::DataSetCreateProps cparms;
    cparms.add(Chunks);
    HighFive::DataSet dataset = file.createDataSet<double>(dataPath, fspace, cparms);
    // Write the data to HDF5
    std::vector<size_t> offset;
    offset.push_back(0);
    offset.push_back(0);
    std::vector<size_t> count;
    count.push_back(x.ncol());
    count.push_back(x.nrow());
    HighFive::Selection selected = dataset.select(offset, count);
    selected.write_raw<double>(x.begin());
    return res;
}

// [[Rcpp::export(.rcpp_spmat_to_h5mat)]]
CharacterVector rcpp_spmat_to_h5mat(const NumericVector& data,
                                    const IntegerVector& rowind,
                                    const IntegerVector& colptr,
                                    unsigned int nrow, unsigned int ncol,
                                    std::string filename, std::string dataPath)
{
    CharacterVector res(2);
    res[0] = filename;
    res[1] = dataPath;

    // Existance of `filename` is already done in R
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    // Create a dataset of native double with the same shape as the matrix
    // // Set the dimensions of the dataset
    std::array<size_t, 2> tDims;
    tDims[0] = ncol;
    tDims[1] = nrow;
    HighFive::DataSpace fspace(tDims);
    // // Set chunking dimension for the new dataset
    std::vector<hsize_t> chunk_dims;
    int colChunkSize = 1000, rowChunkSize = 1000;
    if (colChunkSize > ncol) {
        // Mainly happening in small unit test case, but worth checking
        chunk_dims.push_back(ncol);
    } else {
        chunk_dims.push_back(colChunkSize);
    }
    if (rowChunkSize > nrow) {
        // Mainly happening in small unit test case, but worth checking
        chunk_dims.push_back(nrow);
    } else {
        chunk_dims.push_back(rowChunkSize);
    }
    HighFive::Chunking Chunks(chunk_dims);
    HighFive::DataSetCreateProps cparms;
    cparms.add(Chunks);
    HighFive::DataSet dataset = file.createDataSet<double>(dataPath, fspace, cparms);

    // Write the data to HDF5
    unsigned int numChunks = ncol / colChunkSize;
    if (numChunks * colChunkSize < ncol) numChunks++;
    for (unsigned int i = 0; i < numChunks; ++i) {
        unsigned int startCol = i * colChunkSize;
        unsigned int endCol = (i + 1) * colChunkSize - 1;
        if (endCol > ncol - 1) endCol = ncol - 1;
        unsigned int numCols = endCol - startCol + 1;
        // Create dense matrix from the given sparse information
        arma::mat x = arma::zeros(nrow, numCols);
        arma::uvec colptr_chunk = as<arma::uvec>(colptr[Rcpp::Range(startCol, endCol + 1)]);
        arma::uvec rowind_chunk = as<arma::uvec>(rowind[Rcpp::Range(colptr_chunk[0], colptr_chunk[numCols] - 1)]);
        arma::vec data_chunk = as<arma::vec>(data[Rcpp::Range(colptr_chunk[0], colptr_chunk[numCols] - 1)]);
        colptr_chunk -= colptr_chunk[0];
        arma::sp_mat chunk = arma::sp_mat(rowind_chunk, colptr_chunk, data_chunk, nrow, ncol, true);
        arma::mat chunk_dense = arma::mat(chunk);
        chunk.clear();
        // Write the data to HDF5
        std::vector<size_t> offset;
        offset.push_back(startCol);
        offset.push_back(0);
        std::vector<size_t> count;
        count.push_back(numCols);
        count.push_back(nrow);
        HighFive::Selection selected = dataset.select(offset, count);
        selected.write_raw<double>(chunk_dense.memptr());
    }
    return res;
}

// [[Rcpp::export(.rcpp_spmat_to_h5spmat)]]
CharacterVector rcpp_spmat_to_h5spmat(
    const NumericVector& value, const IntegerVector& rowind,
    const IntegerVector& colptr, unsigned int nrow, unsigned int ncol,
    std::string filename, std::string valuePath, std::string rowindPath,
    std::string colptrPath
) {
    CharacterVector res(4);
    res[0] = filename;
    res[1] = valuePath;
    res[2] = rowindPath;
    res[3] = colptrPath;

    // Create new HDF5 file
    HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
    // Processing colptr
    // // Create dataset for colptr
    std::vector<size_t> p_size;
    p_size.push_back(colptr.size());
    HighFive::DataSpace p_dataspace(p_size);
    std::vector<hsize_t> p_chunksize;
    if (colptr.size() < 2048) {
        p_chunksize.push_back(colptr.size());
    } else {
        p_chunksize.push_back(2048);
    }
    HighFive::Chunking p_Chunks(p_chunksize);
    HighFive::DataSetCreateProps p_cparms;
    p_cparms.add(p_Chunks);
    HighFive::DataSet p_dataset = file.createDataSet<int>(colptrPath, p_dataspace, p_cparms);
    // // Write colptr to HDF5
    std::vector<size_t> p_offset;
    p_offset.push_back(0);
    HighFive::Selection p_selected = p_dataset.select(p_offset, p_size);
    p_selected.write_raw<int>(colptr.begin());

    // Processing rowind
    // // Create dataset for rowind
    std::vector<size_t> i_size;
    i_size.push_back(rowind.size());
    HighFive::DataSpace i_dataspace(i_size);
    std::vector<hsize_t> i_chunksize;
    if (rowind.size() < 2048) {
        i_chunksize.push_back(rowind.size());
    } else {
        i_chunksize.push_back(2048);
    }
    HighFive::Chunking i_Chunks(i_chunksize);
    HighFive::DataSetCreateProps i_cparms;
    i_cparms.add(i_Chunks);
    HighFive::DataSet i_dataset = file.createDataSet<int>(rowindPath, i_dataspace, i_cparms);
    // // Write rowind to HDF5
    std::vector<size_t> i_offset;
    i_offset.push_back(0);
    HighFive::Selection i_selected = i_dataset.select(i_offset, i_size);
    i_selected.write_raw<int>(rowind.begin());

    // Processing data
    // // Create dataset for data
    std::vector<size_t> x_size;
    x_size.push_back(value.size());
    HighFive::DataSpace x_dataspace(x_size);
    std::vector<hsize_t> x_chunksize;
    if (value.size() < 2048) {
        x_chunksize.push_back(value.size());
    } else {
        x_chunksize.push_back(2048);
    }
    HighFive::Chunking x_Chunks(x_chunksize);
    HighFive::DataSetCreateProps x_cparms;
    x_cparms.add(x_Chunks);
    HighFive::DataSet x_dataset = file.createDataSet<double>(valuePath, x_dataspace, x_cparms);
    // // Write data to HDF5
    std::vector<size_t> x_offset;
    x_offset.push_back(0);
    HighFive::Selection x_selected = x_dataset.select(x_offset, x_size);
    x_selected.write_raw<double>(value.begin());

    return res;
}
