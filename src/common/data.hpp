#pragma once

#include "utils.hpp"
#include <filesystem>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

namespace planc {

    class H5Mat : public HighFive::File {
        // A contatiner only for a 2D dense matrix stored in an HDF5 file
        // with accessor function to columns of the matrix that reads and
        // returns a specified chunk of the matrix into memory
        protected:
        std::string filename, datapath;
        std::vector<hsize_t> chunk_dims;

        private:
        static std::string increUniqName(const std::string& base) {
            int suffix = 0;
            std::string tempPath = base + std::to_string(suffix) + ".h5";
            while (std::filesystem::exists(tempPath)) {
                suffix++;
                tempPath = base + std::to_string(suffix) + ".h5";
            }
            return tempPath;
        };

        public:
        //not thread safe
        H5Mat(const std::string& filename, const std::string& datapath) : HighFive::File(filename, HighFive::File::ReadOnly)
        {
            this->filename = filename;
            this->datapath = datapath;
            HighFive::DataSet H5D = this->getDataSet(datapath);
            HighFive::DataSpace dataspace = H5D.getSpace();
            // Get the rank (number of dimensions) of the H5D
            size_t i = dataspace.getNumberDimensions();
            unsigned int rank = i;

            // Check if the rank is 2
            if (rank != 2u) {
                Rcpp::Rcout << "The H5D does not have a rank of 2." << std::endl;
            }
            std::vector<size_t> dims = dataspace.getDimensions();
            // dataspace.close();
            this->n_cols = dims[0];
            this->n_rows = dims[1];

            HighFive::DataSetCreateProps cparms = H5D.getCreatePropertyList();
            this->chunk_dims = HighFive::Chunking(cparms).getDimensions();

            this->colChunkSize = chunk_dims[0];
            this->rowChunkSize = chunk_dims[1];
            // cparms.close();

//#ifdef _VERBOSE
            Rcpp::Rcout << "==H5Mat constructed==" << std::endl
                << "H5File:    " << this->filename << std::endl
                << "Mat path:  " << this->datapath << std::endl
                << "Dimension: " << this->n_rows << " x " << this->n_cols << std::endl;
//#endif
        }

        ~H5Mat() {
            this->flush();
            // this->H5D.close();
            // this->H5F.unlink(this->tempPath);
            // TODO: Have to find a way to unlink the temp transposed matrix
            // this->H5F.close();
        }

        arma::uword n_cols, n_rows, colChunkSize, rowChunkSize;

        arma::mat cols(arma::uword start, arma::uword end) {
            try {
                if (start < 0) {
                    throw std::invalid_argument(
                            "`start` must be an unsigned int, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (start > end) {
                    throw std::invalid_argument(
                            "`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (end >= this->n_cols) {
                    throw std::invalid_argument(
                            "`end` must be less than the number of columns, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
            }  catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);
#else
                throw ex;
#endif
            }
            arma::mat chunk(this->n_rows, end - start + 1);
            std::vector<size_t> offset;
            offset.push_back(start);
            offset.push_back(0);
            std::vector<size_t> count;
            count.push_back(end - start + 1);
            count.push_back(this->n_rows);
            #pragma omp critical
            {
            HighFive::DataSet H5D = this->getDataSet(datapath);
            HighFive::Selection selected = H5D.select(offset, count);
            selected.read<double>(chunk.memptr());
            }
            // dataspace.close();
            // memspace.close();
            return chunk;
        }

        arma::mat cols(arma::uvec index) {
            arma::mat out(this->n_rows, index.size());
            // Identify contiguous ranges from `index` and use .cols(start, end) for each range
            arma::uword start = index[0], end = index[0], curr = index[0];
            arma::uword outStart = 0, outEnd = 0;
            for (arma::uword i = 1; i < index.size(); ++i) {
                curr = index[i];
                try {
                    if (curr > this->n_cols - 1) {
                        throw std::invalid_argument("Index " + std::to_string(curr) + " is out of range.");
                    }
                }  catch(std::exception &ex) {
#ifdef USING_R
                    std::string ex_str = ex.what();
                    Rcpp::stop(ex_str);
#else
                    throw ex;
#endif
                }
                if (curr == end + 1) {
                    // Still contiguous
                    end = curr;
                } else {
                    out.cols(outStart, outEnd) = this->cols(start, end);
                    outStart = outEnd + 1;
                    start = curr;
                    end = curr;
                }
                outEnd++;
            }
            out.cols(outStart, outEnd) = this->cols(start, end);
            return out;
        }

        arma::mat rows(arma::uword start, arma::uword end) {
            try {
                if (start < 0) {
                    throw std::invalid_argument(
                            "`start` must be an unsigned int, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (start > end) {
                    throw std::invalid_argument(
                            "`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (end >= this->n_rows) {
                    throw std::invalid_argument(
                            "`end` must be less than the number of rows, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
            }  catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);
#else
                throw ex;
#endif
            }
            arma::mat chunk(end - start + 1, this->n_cols);
            std::vector<size_t> offset;
            offset.push_back(0);
            offset.push_back(start);
            std::vector<size_t> count;
            count.push_back(this->n_cols);
            count.push_back(end - start + 1);
            #pragma omp critical
            {
            HighFive::DataSet H5D = this->getDataSet(datapath);
            HighFive::Selection selected = H5D.select(offset, count);
            selected.read<double>(chunk.memptr());
            }
            // dataspace.close();
            // memspace.close();
            return chunk;
        }
// not thread-safe
        H5Mat t() {
            // Create new H5 FILE with only the transposed dataset
            std::string tmpfilename = planc::H5Mat::increUniqName(this->filename + ".dense_transposed.");
            HighFive::File tmpfile(tmpfilename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);
            // Set specified chunk dimension for the new dataset
            std::vector<hsize_t> chunk_dims_new;
            if (this->colChunkSize > this->n_rows) {
                // Mainly happening in small unit test case, but worth checking
                chunk_dims_new.push_back(this->n_rows);
            } else {
                chunk_dims_new.push_back(this->colChunkSize);
            }
            if (this->rowChunkSize > this->n_cols) {
                // Mainly happening in small unit test case, but worth checking
                chunk_dims_new.push_back(this->n_cols);
            } else {
                chunk_dims_new.push_back(this->rowChunkSize);
            }
            HighFive::Chunking newChunks(chunk_dims_new);
            HighFive::DataSetCreateProps cparms_new;
            cparms_new.add(newChunks);
            std::array<size_t, 2> tDims{};
            tDims[0] = this->n_rows;
            tDims[1] = this->n_cols;
            HighFive::DataSpace fspace(tDims);
            Rcpp::Rcout << "Creating transposed data at "
                << tmpfilename << ":data" << std::endl;
            HighFive::DataSet H5DT = tmpfile.createDataSet<double>("data", fspace, cparms_new);
            // HighFive::DataSet H5DT = this->createDataSet<double>(tempPath, fspace, cparms_new);
            //cparms_new.close();
            // H5::DataSpace dataspace = H5DT.getSpace();
            unsigned int nChunks = this->n_rows / this->colChunkSize;
            if (nChunks * this->colChunkSize < this->n_rows) nChunks++;
            for (unsigned int i=0; i < nChunks; ++i) {
                arma::uword start = i * this->colChunkSize;
                arma::uword end = (i + 1) * this->colChunkSize - 1;
                if (end > this->n_rows - 1) end = this->n_rows - 1;
                // Read row slices of the current H5D dataset and write to column chunks of the new dataset
                arma::mat chunk = this->rows(start, end).t();
                std::vector<size_t> offset;
                offset.push_back(start);
                offset.push_back(0);
                std::vector<size_t> count;
                count.push_back(end - start + 1);
                count.push_back(this->n_cols);
                // H5DT.select(offset, count).write<double>(*chunk.memptr());
                HighFive::Selection selected = H5DT.select(offset, count);
                selected.write_raw<double>(chunk.memptr());
                // memspace.close();
            }
            // fspace.close();
            // H5DT.close();
            tmpfile.flush();
            H5Mat transposedMat(tmpfilename, "data");
            return transposedMat;
        } // End of H5Mat.t()

        double normF() {
            arma::uword nChunks = this->n_cols / this->colChunkSize;
            if (nChunks * this->colChunkSize < this->n_cols) nChunks++;
            double norm = 0;
            for (arma::uword i = 0; i < nChunks; ++i) {
                arma::uword start = i * this->colChunkSize;
                arma::uword end = (i + 1) * this->colChunkSize - 1;
                if (end > this->n_cols - 1) end = this->n_cols - 1;
                arma::mat chunk = this->cols(start, end);
                norm += arma::accu(chunk % chunk);
            }
            return std::sqrt(norm);
        }
    }; // End of class H5Mat

    class H5SpMat : public HighFive::File  {
        protected:
        std::string filename, xPath, iPath, pPath;
        arma::uword x_chunksize, i_chunksize, p_chunksize;

        private:
        // The `start` and `end` refer to the start and end of the corresponding
        // arrays, but not the indices of the sparse matrix
        arma::uvec getPByRange(arma::uword start, arma::uword end) {
            arma::uvec p(end - start + 1);
            std::vector<size_t> p_start;
            p_start.push_back(start);
            std::vector<size_t> p_count;
            p_count.push_back(end - start + 1);
            #pragma omp critical
            {
            HighFive::DataSet H5D_P = this->getDataSet(pPath);
            HighFive::Selection selected_p = H5D_P.select(p_start, p_count);
            selected_p.read<arma::uword>(p.memptr());
            }
            // pDataspace.close();
            // pMemspace.close();
            return p;
        }

        arma::uvec getIByRange(arma::uword start, arma::uword end) {
            arma::uvec i(end - start + 1);
            std::vector<size_t> i_start;
            i_start.push_back(start);
            std::vector<size_t> i_count;
            i_count.push_back(end - start + 1);
            #pragma omp critical
            {
            HighFive::DataSet H5D_I = this->getDataSet(iPath);
            HighFive::Selection selected_i = H5D_I.select(i_start, i_count);
            selected_i.read<arma::uword>(i.memptr());
            }
            // iDataspace.close();
            // iMemspace.close();
            return i;
        }

        arma::vec getXByRange(arma::uword start, arma::uword end) {
            arma::vec x(end - start + 1);
            std::vector<size_t> x_start;
            x_start.push_back(start);
            std::vector<size_t> x_count;
            x_count.push_back(end - start + 1);
            #pragma omp critical
            {
            HighFive::DataSet H5D_X = this->getDataSet(xPath);
            HighFive::Selection selected_x = H5D_X.select(x_start, x_count);
            selected_x.read<double>(x.memptr());
            }
            // xDataspace.close();
            // xMemspace.close();
            return x;
        }

        static std::string increUniqName(const std::string& base) {
            int suffix = 0;
            std::string tempPath = base + std::to_string(suffix) + ".h5";
            while (std::filesystem::exists(tempPath)) {
                suffix++;
                tempPath = base + std::to_string(suffix) + ".h5";
            }
            return tempPath;
        }

        public:
// not thread safe
        H5SpMat(const std::string& filename, const std::string& iPath, const std::string& pPath,
                const std::string& xPath, arma::uword n_rows, arma::uword n_cols) :
                HighFive::File(filename, HighFive::File::ReadWrite) {
            this->filename = filename;
            this->iPath = iPath;
            this->pPath = pPath;
            this->xPath = xPath;
            this->n_rows = n_rows;
            this->n_cols = n_cols;
            HighFive::DataSet H5D_X = this->getDataSet(xPath);
            HighFive::DataSetCreateProps x_cparms = H5D_X.getCreatePropertyList();
            std::vector<hsize_t> x_chunkdim = HighFive::Chunking(x_cparms).getDimensions();
            this->x_chunksize = x_chunkdim[0];
            // x_cparms.close();

            HighFive::DataSpace xDataspace = H5D_X.getSpace();
            std::vector<size_t> xDims;
            xDims = xDataspace.getDimensions();
            // xDataspace.close();
            this->nnz = xDims[0];

            HighFive::DataSet H5D_I = this->getDataSet(iPath);
            HighFive::DataSetCreateProps i_cparms = H5D_I.getCreatePropertyList();
            std::vector<hsize_t> i_chunkdim = HighFive::Chunking(i_cparms).getDimensions();
            this->i_chunksize = i_chunkdim[0];
            // i_cparms.close();
            HighFive::DataSet H5D_P = this->getDataSet(pPath);
            HighFive::DataSetCreateProps p_cparms = H5D_P.getCreatePropertyList();
            std::vector<hsize_t> p_chunkdim = HighFive::Chunking(p_cparms).getDimensions();
            this->p_chunksize = p_chunkdim[0];
            // p_cparms.close();
// #ifdef _VERBOSE
            Rcpp::Rcout << "==H5SpMat constructed==" << std::endl
                << "H5File:    " << filename << std::endl
                << "colptr path:  " << pPath << std::endl
                << "rowind path:  " << iPath << std::endl
                << "value path:   " << xPath << std::endl
                << "Dimension: " << n_rows << " x " << n_cols << std::endl;
// #endif
        }
        ~H5SpMat()
        {
            this->flush();
            // this->H5D.close();
            // this->H5F.unlink(this->tempPath);
            // TODO: Have to find a way to unlink the temp transposed matrix
            // this->H5F.close();
        }
        arma::uword n_rows, n_cols, nnz;

        arma::sp_mat cols(arma::uword start, arma::uword end) {
            try {
                if (start < 0) {
                    throw std::invalid_argument(
                            "`start` must be an unsigned int, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (start > end) {
                    throw std::invalid_argument(
                            "`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
                if (end >= this->n_cols) {
                    throw std::invalid_argument(
                            "`end` must be less than the number of columns, got (" + std::to_string(start) + ", " +
                            std::to_string(end) + ").");
                }
            }  catch(std::exception &ex) {
#ifdef USING_R
                std::string ex_str = ex.what();
                Rcpp::stop(ex_str);
#else
                throw ex;
#endif
            }

            // Construct subsetted colptr
            arma::uvec colptr = this->getPByRange(start, end + 1);

            // Subset rowind and value according to colptr
            arma::uvec rowind = this->getIByRange(colptr[0], colptr[end - start + 1] - 1);
            arma::vec value = this->getXByRange(colptr[0], colptr[end - start + 1] - 1);

            colptr -= colptr[0];
            // Construct the subset sparse matrix
            arma::sp_mat chunk(rowind, colptr, value, this->n_rows, end - start + 1);

            return chunk;
        }

        arma::sp_mat cols(arma::uvec index) {
            arma::sp_mat out(this->n_rows, index.size());
            // Identify contiguous ranges from `index` and use .cols(start, end) for each range
            arma::uword start = index[0], end = index[0], curr = index[0];
            arma::uword outStart = 0, outEnd = 0;
            for (arma::uword i = 1; i < index.size(); ++i) {
                curr = index[i];
                try {
                    if (curr > this->n_cols - 1) {
                        throw std::invalid_argument("Index " + std::to_string(curr) + " is out of range.");
                    }
                }  catch(std::exception &ex) {
#ifdef USING_R
                    std::string ex_str = ex.what();
                    Rcpp::stop(ex_str);
#else
                    throw ex;
#endif
                }
                if (curr == end + 1) {
                    // Still contiguous
                    end = curr;
                } else {
                    out.cols(outStart, outEnd) = this->cols(start, end);
                    outStart = outEnd + 1;
                    start = curr;
                    end = curr;
                }
                outEnd++;
            }
            out.cols(outStart, outEnd) = this->cols(start, end);
            return out;
        }

        // arma::sp_mat cols(arma::uvec index) {
        //     arma::sp_mat out(this->n_rows, index.size());
        //     arma::uvec colptr = this->getPByRange(0, this->n_cols);
        //     arma::uword idx;
        //     for (arma::uword i = 0; i < index.size(); ++i) {
        //         idx = index[i];
        //         if (idx >= this->n_cols) {
        //             throw std::invalid_argument("Index " + std::to_string(idx) + " is out of range.");
        //         }
        //         arma::uvec rowind = this->getIByRange(colptr[idx], colptr[idx + 1] - 1);
        //         arma::vec value = this->getXByRange(colptr[idx], colptr[idx + 1] - 1);
        //         arma::uvec colptr_i = arma::zeros<arma::uvec>(2);
        //         colptr_i[1] = rowind.size();
        //         out.col(i) = arma::sp_mat(rowind, colptr_i, value, this->n_rows, 1);
        //     }
        //     return out;
        // }
// not thread safe
        H5SpMat t() {
            Rcpp::Rcout << "Creating on-disk transposition of the sparse matrix, which is currently poorly supported and slow" << std::endl;
            // Create new H5 FILE with only the transposed dataset
            std::string tmpfilename = planc::H5SpMat::increUniqName(this->filename + ".sparse_transposed.");
            HighFive::File tmpfile(tmpfilename, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate);

            // ================= Create the colptr of the transposed matrix =================
            // Read all colptr into memory, it's not big :)
            arma::uvec colptr = this->getPByRange(0, this->n_cols);
            // Go through chunks of rowind and initialize the transposed colptr
            arma::uvec colptrT(this->n_rows + 1);
            arma::uword nChunk = this->nnz / this->i_chunksize;
            if (nChunk * this->i_chunksize < this->nnz) nChunk++;
            for (arma::uword i = 0; i < nChunk; i++) {
                arma::uword start = i * this->i_chunksize;
                arma::uword end = (i + 1) * this->i_chunksize - 1;
                if (end > this->nnz - 1) end = this->nnz - 1;
                arma::uvec rowind_chunk = this->getIByRange(start, end);
                for (unsigned int j : rowind_chunk) {
                    // rowind_chunk[j] gives the number referring to which row of the original
                    // or which column of the transposed matrix have a non-zero value
                    colptrT[j + 1]++;
                }
            }
            // cumsum to get the final colptrT
            colptrT = arma::cumsum(colptrT);

            // Write colptrT to HDF5 file
            // std::string ptempPath = this->increUniqName(this->pPath + "_transposed_");
            std::vector<hsize_t> inputp_chunksize;
            if (this->p_chunksize > this->n_rows+1) {
                // Mainly happening in small unit test case, but worth checking
                inputp_chunksize.push_back(this->n_rows + 1);
            } else {
                inputp_chunksize.push_back(this->p_chunksize);
            }
            HighFive::Chunking p_newChunks(inputp_chunksize);
            HighFive::DataSetCreateProps p_cparms_new;
            p_cparms_new.add(p_newChunks);

            std::vector<size_t> p_new_size;
            p_new_size.push_back(this->n_rows + 1);
            HighFive::DataSpace p_new_dataspace(p_new_size);
            HighFive::DataSet H5D_PT = tmpfile.createDataSet<arma::uword>("colptr", p_new_dataspace, p_cparms_new);

            std::vector<size_t> offset;
            offset.push_back(0);
            HighFive::Selection p_selected = H5D_PT.select(offset, p_new_size);
            p_selected.write_raw<arma::uword>(colptrT.memptr());
            // p_new_dataspace.close();
            // p_new_memspace.close();
            // H5D_PT.close();
            // ================= self.t.colptr done ! =================

            // ================= Create the rowind and value of the transposed matrix =================
            // Pre-create the rowind and value of the transposed matrix
            // Create rowind.T
            // std::string itempPath = this->increUniqName(this->iPath + "_transposed_");
            std::vector<hsize_t> inputi_chunksize;
            inputi_chunksize.push_back(this->i_chunksize);
            HighFive::Chunking i_newChunks(inputi_chunksize);
            HighFive::DataSetCreateProps i_cparms_new;
            i_cparms_new.add(i_newChunks);

            std::vector<size_t> i_new_size;
            i_new_size.push_back(this->nnz);
            HighFive::DataSpace i_new_dataspace(i_new_size);
            // HighFive::DataSet H5D_IT = this->createDataSet<arma::uword>(itempPath,i_new_dataspace, i_cparms_new);
            HighFive::DataSet H5D_IT = tmpfile.createDataSet<arma::uword>("rowind", i_new_dataspace, i_cparms_new);
            // i_cparms_new.close();
            // Create value.T
            // std::string xtempPath = this->increUniqName(this->xPath + "_transposed_");
            std::vector<hsize_t> inputx_chunksize;
            inputx_chunksize.push_back(this->x_chunksize);
            HighFive::Chunking x_newChunks(inputx_chunksize);
            HighFive::DataSetCreateProps x_cparms_new;
            x_cparms_new.add(x_newChunks);

            std::vector<size_t> x_new_size;
            x_new_size.push_back(this->nnz);
            HighFive::DataSpace x_new_dataspace(x_new_size);
            // HighFive::DataSet H5D_XT = this->createDataSet<double>(xtempPath, x_new_dataspace, x_cparms_new);
            HighFive::DataSet H5D_XT = tmpfile.createDataSet<double>("value", x_new_dataspace, x_cparms_new);
            // x_cparms_new.close();

            arma::uvec colptrT_start = colptrT;
            Progress p(this->n_cols, true);
            for (arma::uword i = 0; i < this->n_cols; i++) {
                // Go through each column of the original matrix and fill in the
                // rowind and value of the transposed matrix
                // For each of the original column, get the corresponding fragment from original rowind and value
                arma::uvec rowind_ori_col = this->getIByRange(colptr[i], colptr[i + 1] - 1);
                arma::vec value_ori_col = this->getXByRange(colptr[i], colptr[i + 1] - 1);
                // With the frangment of rowind we get, we know which rows of the
                // original matrix or cols of the transposed matrix have non-zero values
                // Then it is mapped to the positions in colptrT. The values fetched here
                // then map to the positions in rowindT and valueT. The current col index, `i`,
                // will be the value to be filled in rowindT, and the fragment from `value` can
                // be simply copied to valueT at the same positions.
                arma::uvec nnz_idx = colptrT_start.elem(rowind_ori_col);
                arma::uword nnz_idx_size = nnz_idx.size();
                arma::uvec it_value(nnz_idx_size);
                it_value.fill(i);

                std::vector<size_t> coord;
                for (arma::uword j = 0; j < nnz_idx_size; j++) {
                    coord.push_back(nnz_idx[j]);
                }
                // size_t count = nnz_idx_size;
                // HighFive::DataSpace xt_writeDataSpace(1, count);
                // HighFive::DataSpace it_writeDataSpace(1, count);

                HighFive::Selection i_selected = H5D_IT.select(coord);
                HighFive::Selection x_selected = H5D_XT.select(coord);

                i_selected.write_raw<arma::uword>(it_value.memptr());
                x_selected.write_raw<double>(value_ori_col.memptr());

                // Increment it, so that next time when fetching the number `nnz_idx` of the same new-column/old-row,
                // it know the previous position has been filled
                colptrT_start.elem(rowind_ori_col) += 1;

                // delete[] coord;
                p.increment();
            }
            // i_new_dataspace.close();
            // x_new_dataspace.close();
            // H5D_IT.close();
            // H5D_XT.close();
            // ================= self.t.rowind and self.t.value done ! =================
            tmpfile.flush();
            H5SpMat transposedMat(tmpfilename, "rowind", "colptr", "value", this->n_cols, this->n_rows);
            // H5SpMat transposedMat(this->filename, itempPath, ptempPath, xtempPath, this->n_cols, this->n_rows);
            return transposedMat;
        } // End of H5SpMat.t()

        double normF() {
            arma::uword nChunks = this->nnz / this->x_chunksize;
            if (nChunks * this->x_chunksize < this->nnz) nChunks++;
            double norm = 0;
            for (arma::uword i = 0; i < nChunks; ++i) {
                arma::uword start = i * this->x_chunksize;
                arma::uword end = (i + 1) * this->x_chunksize - 1;
                if (end > this->nnz - 1) end = this->nnz - 1;
                arma::vec chunk = this->getXByRange(start, end);
                norm += arma::accu(chunk % chunk);
            }
            return std::sqrt(norm);
        }
    }; // End of class H5SpMat
} // End of namespace planc
