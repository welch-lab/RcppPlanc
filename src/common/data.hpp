#pragma once

#include "utils.hpp"

namespace planc {

    class H5Mat {
        // A contatiner only for a 2D dense matrix stored in an HDF5 file
        // with accessor function to columns of the matrix that reads and
        // returns a specified chunk of the matrix into memory
        protected:
        std::string filename, datapath;//, tempPath;
        H5::H5File H5F;
        H5::DataSet H5D;

        private:
        std::string increUniqName(std::string base) {
            int suffix = 0;
            std::string tempPath = base + std::to_string(suffix);
            while (this->H5F.exists(tempPath)) {
                suffix++;
                tempPath = base + std::to_string(suffix);
            }
            return tempPath;
        }

        public:
        H5Mat(std::string filename, std::string datapath)
            : filename(filename), datapath(datapath) {
            H5::H5File file(filename, H5F_ACC_RDWR);
            this->H5F = file;
            this->H5D = H5::DataSet(this->H5F.openDataSet(datapath));
            H5::DataSpace dataspace = this->H5D.getSpace();
            // Get the rank (number of dimensions) of the H5D
            int rank = dataspace.getSimpleExtentNdims();

            // Check if the rank is 2
            if (rank != 2) {
                std::cout << "The H5D does not have a rank of 2." << std::endl;
            }
            hsize_t dims[2];
            dataspace.getSimpleExtentDims(dims);
            dataspace.close();
            this->n_cols = dims[0];
            this->n_rows = dims[1];

            H5::DSetCreatPropList cparms = this->H5D.getCreatePlist();
            hsize_t chunk_dims[2];
            cparms.getChunk( 2, chunk_dims);
            this->colChunkSize = chunk_dims[0];
            this->rowChunkSize = chunk_dims[0];
            cparms.close();

//#ifdef _VERBOSE
            std::cout << "==H5Mat constructed==" << std::endl
                << "H5File:    " << filename << std::endl
                << "Mat path:  " << datapath << std::endl
                << "Dimension: " << n_rows << " x " << n_cols << std::endl;
//#endif
        }

        ~H5Mat() {
            this->H5D.close();
            // this->H5F.unlink(this->tempPath);
            // TODO: Have to find a way to unlink the temp transposed matrix
            this->H5F.close();
        }

        arma::uword n_cols, n_rows, colChunkSize, rowChunkSize;

        arma::mat cols(arma::uword start, arma::uword end) {
            if (start < 0) {
                throw std::invalid_argument("`start` must be an unsigned int, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (start > end) {
                throw std::invalid_argument("`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (end >= this->n_cols) {
                throw std::invalid_argument("`end` must be less than the number of columns, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            arma::mat chunk(this->n_rows, end - start + 1);
            H5::DataSpace dataspace = this->H5D.getSpace();
            hsize_t offset[2];
            offset[0] = start;
            offset[1] = 0;
            hsize_t count[2];
            count[0] = end - start + 1;
            count[1] = this->n_rows;
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
            H5::DataSpace memspace(2, count);
            this->H5D.read(chunk.memptr(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
            dataspace.close();
            memspace.close();
            return chunk;
        }

        arma::mat cols(arma::uvec index) {
            arma::mat out(this->n_rows, index.size());
            // Identify contiguous ranges from `index` and use .cols(start, end) for each range
            int start = index[0], end = index[0], curr = index[0];
            int outStart = 0, outEnd = 0;
            for (arma::uword i = 1; i < index.size(); ++i) {
                curr = index[i];
                if (curr > this->n_cols - 1) {
                    throw std::invalid_argument("Index " + std::to_string(curr) + " is out of range.");
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
            if (start < 0) {
                throw std::invalid_argument("`start` must be an unsigned int, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (start > end) {
                throw std::invalid_argument("`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (end >= this->n_rows) {
                throw std::invalid_argument("`end` must be less than the number of rows, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            arma::mat chunk(end - start + 1, this->n_cols);
            H5::DataSpace dataspace = this->H5D.getSpace();
            hsize_t offset[2];
            offset[0] = 0;
            offset[1] = start;
            hsize_t count[2];
            count[0] = this->n_cols;
            count[1] = end - start + 1;
            dataspace.selectHyperslab(H5S_SELECT_SET, count, offset);
            H5::DataSpace memspace(2, count);
            this->H5D.read(chunk.memptr(), H5::PredType::NATIVE_DOUBLE, memspace, dataspace);
            dataspace.close();
            memspace.close();
            return chunk;
        }

        H5Mat t() {
            // Create new H5 dataset at a unique temporary path, with incrementing suffix
            std::string tempPath = this->increUniqName(this->datapath + "_transposed_");
            // this->tempPath = tempPath;
            // Set specified chunk dimension for the new dataset
            hsize_t chunk_dims_new[2];
            chunk_dims_new[0] = this->colChunkSize;
            chunk_dims_new[1] = this->n_cols;
            H5::DSetCreatPropList cparms_new;
            cparms_new.setChunk(2, chunk_dims_new);

            hsize_t fdim[2];
            fdim[0] = this->n_rows;
            fdim[1] = this->n_cols;
            H5::DataSpace fspace(2, fdim);
            std::cout << "Creating transposed data at " << tempPath << std::endl;
            H5::DataSet H5DT = this->H5F.createDataSet(tempPath, H5::PredType::NATIVE_DOUBLE, fspace, cparms_new);
            cparms_new.close();

            // H5::DataSpace dataspace = H5DT.getSpace();
            for (int i=0; i < this->n_rows / this->colChunkSize; ++i) {
                arma::uword start = i * this->colChunkSize;
                arma::uword end = (i + 1) * this->colChunkSize - 1;
                if (end > this->n_rows - 1) end = this->n_rows - 1;
                // Read row slices of the current H5D dataset and write to column chunks of the new dataset
                arma::mat chunk = this->rows(start, end).t();
                hsize_t offset[2];
                offset[0] = start;
                offset[1] = 0;
                hsize_t count[2];
                count[0] = end - start + 1;
                count[1] = this->n_cols;
                fspace.selectHyperslab(H5S_SELECT_SET, count, offset);
                // Write the chunk to dataspace
                H5::DataSpace memspace(2, count);
                H5DT.write(chunk.memptr(), H5::PredType::NATIVE_DOUBLE, memspace, fspace);
                memspace.close();
            }
            fspace.close();
            H5DT.close();
            H5Mat transposedMat(this->filename, tempPath);
            return transposedMat;
        } // End of H5Mat.t()

    }; // End of class H5Mat

    class H5SpMat {
        protected:
        H5::H5File H5F;
        std::string filename, xPath, iPath, pPath;
        H5::DataSet H5D_X, H5D_I, H5D_P;
        arma::uword x_chunksize, i_chunksize, p_chunksize;

        private:
        // The `start` and `end` refer to the start and end of the corresponding
        // arrays, but not the indices of the sparse matrix
        arma::uvec getPByRange(arma::uword start, arma::uword end) {
            arma::uvec p(end - start + 1);
            hsize_t p_start[1] = {start};
            hsize_t p_count[1] = {end - start + 1};
            H5::DataSpace pDataspace = H5D_P.getSpace();
            pDataspace.selectHyperslab(H5S_SELECT_SET, p_count, p_start);
            H5::DataSpace pMemspace(1, p_count);
            H5D_P.read(p.memptr(), H5::PredType::NATIVE_UINT, pMemspace, pDataspace);
            pDataspace.close();
            pMemspace.close();
            return p;
        }

        arma::uvec getIByRange(arma::uword start, arma::uword end) {
            arma::uvec i(end - start + 1);
            hsize_t i_start[1] = {start};
            hsize_t i_count[1] = {end - start + 1};
            H5::DataSpace iDataspace = H5D_I.getSpace();
            iDataspace.selectHyperslab(H5S_SELECT_SET, i_count, i_start);
            H5::DataSpace iMemspace(1, i_count);
            H5D_I.read(i.memptr(), H5::PredType::NATIVE_UINT, iMemspace, iDataspace);
            iDataspace.close();
            iMemspace.close();
            return i;
        }

        arma::vec getXByRange(arma::uword start, arma::uword end) {
            arma::vec x(end - start + 1);
            hsize_t x_start[1] = {start};
            hsize_t x_count[1] = {end - start + 1};
            H5::DataSpace xDataspace = H5D_X.getSpace();
            xDataspace.selectHyperslab(H5S_SELECT_SET, x_count, x_start);
            H5::DataSpace xMemspace(1, x_count);
            H5D_X.read(x.memptr(), H5::PredType::NATIVE_DOUBLE, xMemspace, xDataspace);
            xDataspace.close();
            xMemspace.close();
            return x;
        }

        std::string increUniqName(std::string base) {
            int suffix = 0;
            std::string tempPath = base + std::to_string(suffix);
            while (this->H5F.exists(tempPath)) {
                suffix++;
                tempPath = base + std::to_string(suffix);
            }
            return tempPath;
        }

        public:
        H5SpMat(std::string filename, std::string iPath, std::string pPath,
                std::string xPath, arma::uword n_rows, arma::uword n_cols)
        : filename(filename), xPath(xPath), iPath(iPath), pPath(pPath),
            n_rows(n_rows), n_cols(n_cols) {
            H5::H5File file(filename, H5F_ACC_RDWR);
            this->H5F = file;
            H5D_X = H5::DataSet(file.openDataSet(xPath));
            H5::DSetCreatPropList x_cparms = this->H5D_X.getCreatePlist();
            hsize_t x_chunkdim[1];
            x_cparms.getChunk(1, x_chunkdim);
            this->x_chunksize = x_chunkdim[0];
            x_cparms.close();

            H5::DataSpace xDataspace = H5D_X.getSpace();
            hsize_t xDims[1];
            xDataspace.getSimpleExtentDims(xDims);
            xDataspace.close();
            this->nnz = xDims[0];

            H5D_I = H5::DataSet(file.openDataSet(iPath));
            H5::DSetCreatPropList i_cparms = this->H5D_I.getCreatePlist();
            hsize_t i_chunkdim[1];
            i_cparms.getChunk(1, i_chunkdim);
            this->i_chunksize = i_chunkdim[0];
            i_cparms.close();

            H5D_P = H5::DataSet(file.openDataSet(pPath));
            H5::DSetCreatPropList p_cparms = this->H5D_P.getCreatePlist();
            hsize_t p_chunkdim[1];
            p_cparms.getChunk(1, p_chunkdim);
            this->p_chunksize = p_chunkdim[0];
            p_cparms.close();
// #ifdef _VERBOSE
            std::cout << "==H5SpMat constructed==" << std::endl
                << "H5File:    " << filename << std::endl
                << "colptr path:  " << pPath << std::endl
                << "rowind path:  " << iPath << std::endl
                << "value path:   " << xPath << std::endl
                << "Dimension: " << n_rows << " x " << n_cols << std::endl;
// #endif
        }

        arma::uword n_rows, n_cols, nnz;

        arma::sp_mat cols(arma::uword start, arma::uword end) {
            if (start < 0) {
                throw std::invalid_argument("`start` must be an unsigned int, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (start > end) {
                throw std::invalid_argument("`start` must be less than or equal to `end`, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
            }
            if (end >= this->n_cols) {
                throw std::invalid_argument("`end` must be less than the number of columns, got (" + std::to_string(start) + ", " + std::to_string(end) + ").");
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
            int start = index[0], end = index[0], curr = index[0];
            int outStart = 0, outEnd = 0;
            for (arma::uword i = 1; i < index.size(); ++i) {
                curr = index[i];
                if (curr > this->n_cols - 1) {
                    throw std::invalid_argument("Index " + std::to_string(curr) + " is out of range.");
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

        H5SpMat t() {
            // ================= Create the colptr of the transposed matrix =================
            // Read all colptr into memory, it's not big :)
            std::cout << "Creating on-disk transposition of the sparse matrix, which is currently poorly supported and slow" << std::endl;
            arma::uvec colptr = this->getPByRange(0, this->n_cols);
            // Go through chunks of rowind and initialize the transposed colptr
            arma::uvec colptrT(this->n_rows + 1);
            int nChunk = this->nnz / this->i_chunksize;
            if (nChunk * this->i_chunksize < this->nnz) nChunk++;
            for (arma::uword i = 0; i < nChunk; i++) {
                arma::uword start = i * this->i_chunksize;
                arma::uword end = (i + 1) * this->i_chunksize - 1;
                if (end > this->nnz - 1) end = this->nnz - 1;
                arma::uvec rowind_chunk = this->getIByRange(start, end);
                for (arma::uword j = 0; j < rowind_chunk.size(); j++) {
                    // rowind_chunk[j] gives the number referring to which row of the original
                    // or which column of the transposed matrix have a non-zero value
                    colptrT[rowind_chunk[j] + 1]++;
                }
            }
            // cumsum to get the final colptrT
            colptrT = arma::cumsum(colptrT);

            // Write colptrT to HDF5 file
            std::string ptempPath = this->increUniqName(this->pPath + "_transposed_");
            hsize_t pChunk[1] = { this->p_chunksize };
            H5::DSetCreatPropList p_cparms_new;
            p_cparms_new.setChunk(1, pChunk);

            hsize_t p_new_memspacesize[1] = { this->n_rows + 1 };
            H5::DataSpace p_new_dataspace(1, p_new_memspacesize);
            H5::DataSet H5D_PT = this->H5F.createDataSet(ptempPath, H5::PredType::NATIVE_UINT, p_new_dataspace, p_cparms_new);
            p_cparms_new.close();
            hsize_t offset[1] = { 0 };
            p_new_dataspace.selectHyperslab(H5S_SELECT_SET, p_new_memspacesize, offset);
            // Write the chunk to dataspace
            H5::DataSpace p_new_memspace(1, p_new_memspacesize);
            H5D_PT.write(colptrT.memptr(), H5::PredType::NATIVE_UINT, p_new_memspace, p_new_dataspace);
            p_new_dataspace.close();
            p_new_memspace.close();
            H5D_PT.close();
            // ================= self.t.colptr done ! =================

            // ================= Create the rowind and value of the transposed matrix =================
            // Pre-create the rowind and value of the transposed matrix
            // Create rowind.T
            std::string itempPath = this->increUniqName(this->iPath + "_transposed_");
            hsize_t iChunk[1] = { this->i_chunksize };
            H5::DSetCreatPropList i_cparms_new;
            i_cparms_new.setChunk(1, iChunk);
            hsize_t i_new_memspacesize[1] = { this->nnz };
            H5::DataSpace i_new_dataspace(1, i_new_memspacesize); // Gonna be used later for writing
            H5::DataSet H5D_IT = this->H5F.createDataSet(itempPath, H5::PredType::NATIVE_UINT, i_new_dataspace, i_cparms_new);
            i_cparms_new.close();
            // Create value.T
            std::string xtempPath = this->increUniqName(this->xPath + "_transposed_");
            hsize_t xChunk[1] = { this->x_chunksize };
            H5::DSetCreatPropList x_cparms_new;
            x_cparms_new.setChunk(1, xChunk);
            hsize_t x_new_memspacesize[1] = { this->nnz };
            H5::DataSpace x_new_dataspace(1, x_new_memspacesize); // Gonna be used later for writing
            H5::DataSet H5D_XT = this->H5F.createDataSet(xtempPath, H5::PredType::NATIVE_DOUBLE, x_new_dataspace, x_cparms_new);
            x_cparms_new.close();

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

                hsize_t coord[nnz_idx_size];
                for (arma::uword j = 0; j < nnz_idx_size; j++) {
                    coord[j] = nnz_idx[j];
                }
                hsize_t count[1] = { nnz_idx_size };
                H5::DataSpace xt_writeDataSpace(1, count);
                H5::DataSpace it_writeDataSpace(1, count);
                i_new_dataspace.selectElements(H5S_SELECT_SET, nnz_idx_size, coord);
                x_new_dataspace.selectElements(H5S_SELECT_SET, nnz_idx_size, coord);

                H5D_IT.write(it_value.memptr(), H5::PredType::NATIVE_UINT, it_writeDataSpace, i_new_dataspace);
                H5D_XT.write(value_ori_col.memptr(), H5::PredType::NATIVE_DOUBLE, xt_writeDataSpace, x_new_dataspace);

                // Increment it, so that next time when fetching the number `nnz_idx` of the same new-column/old-row,
                // it know the previous position has been filled
                colptrT_start.elem(rowind_ori_col) += 1;

                // delete[] coord;
                p.increment();
            }
            i_new_dataspace.close();
            x_new_dataspace.close();
            H5D_IT.close();
            H5D_XT.close();
            // ================= self.t.rowind and self.t.value done ! =================

            H5SpMat transposedMat(this->filename, itempPath, ptempPath, xtempPath, this->n_cols, this->n_rows);
            return transposedMat;
        } // End of H5SpMat.t()
    }; // End of class H5SpMat
} // End of namespace planc
