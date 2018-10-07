#' @title Solves the linear assignment problem using the Jonker-Vogenant algorithm
#'
#' @description Find a set of vertices pairs in the order of goodness of matching according to a
#' specified measure.
#'
#' @param cost A non-negative matrix-like object that can be coerced to a matrix
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the 
#' costs are maximized
#'
#' @return The assignment of rows to columns as an integer vector
#'
#' @export
lapjv <- function(cost, maximize = FALSE) {
    m <- max(cost)
    n <- nrow(cost)

    cost <- rbind(cbind(as.matrix(cost), m + m * runif(n)),
            m + m * runif(n +1))
    cost[n + 1, n + 1] <- m ^ 3

    ind <- cpp_lapjv(
        ,
        maximize)
    ind[1:n]
}


#' @title Solves the linear assignment problem using the LAPMOD algorithm
#'
#' @description Find a set of vertices pairs in the order of goodness of matching according to a
#' specified measure.
#'
#' @param n number of rows in the cost matrix
#' @param cc vector of all finite elements of the assignement cost matri
#' @param ii vector of indices of the zero indexed  row starts in cc. The following must hold
#'  ii[1] = 0 and ii[n+2] = length(cc).
#' @param kk 0-based column numbers for each finite cost in the matrix, 
#'  i.e., kk must be in 0:(nrow(.)-1).
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the 
#'  costs are maximized
#'    
#' @return The assignment of rows to columns as an integer vector
#'
#' @export
lapmod_index <- function(n, cc, ii, kk, maximize = FALSE) {
    # warning("Currently does not produce expected answers for all matrices.")
    cpp_lapmod(n, cc, ii, kk, maximize)
}

#' @title Solves the linear assignment problem using the LAPMOD algorithm
#'
#' @description Find a set of vertices pairs in the order of goodness of matching according to a
#' specified measure.
#'
#' @param cost A non-negative CsparseMatrix object from the Matrix package
#' @param maximize If FALSE (default) then costs are minimized and if TRUE the 
#' costs are maximized
#'
#' @return The assignment of rows to columns as an integer vector
#'
#' @export
lapmod <- function(spmat, maximize = FALSE){
    # warning("Currently does not produce expected answers for all matrices.")
    n <- nrow(spmat)
    m <- max(spmat@x)
    browser()
    spmat <- rbind2(cbind2(spmat, m + m * runif(n)),
        m + m * runif(n+1))
    ind <- cpp_lapmod(n, spmat@x,
        spmat@p, spmat@i, maximize)
    ind[1:n]
}

