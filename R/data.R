#' Example single-cell transcriptomic data in sparse form
#' 
#' @description
#' The two datasets, namingly \code{ctrl.sparse} and \code{stim.sparse}, are 
#' single-cell transcriptomic data preprocessed and subsampled from the study
#' of Hyun Min Kang and et al., Nat Biotech., 2018. The raw datasets were two
#' sparse matrices of integer values indicating the counts of genes (rows) per
#' cell (columns). We normalized each column of both matrices by its library 
#' size (sum), and selected common variable genes across the datasets. Finally,
#' we scaled the genes without centering them, in order to keep the 
#' non-negativity. The processed datasets were then randomly subsampled for a
#' minimal example. 
#' @source https://www.nature.com/articles/nbt.4042
#' @rdname data
#' @importFrom Matrix Matrix
"ctrl.sparse"

#' @rdname data
#' @importFrom Matrix Matrix
"stim.sparse"
