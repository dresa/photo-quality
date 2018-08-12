# Image MDL: Minimum Description Length
#
# MDL is a model-selection method: simple model that still helps to
# compress data significantly. To use MDL, we need to estimate the
# number of bits it takes to encrypt the model and data.
#
# Esa Junttila, 2018-08-12


# RLE: Run-length encoding.
#####
# ?rle --> already implemented


# Quadtree
##########
# not implemented yet


# Improved 2D RLE:
#################
# Based on the following article:
#  "An Improved Two-Dimensional Run-Length Encoding Scheme and Its Application"
# --> not implemented yet


#' @rdname rle2D
#' @name rle2D
#' @aliases rle2d
#' @aliases invRle2d
#' @title 2D Run-Length Encoding (2D RLE).
#' @description
#' 2D RLE is a 2D matrix version of ordinary Run-Length Encoding (1D RLE),
#' which is a popular compression method for repetitive vectors.
#' 2D RLE uses the same idea to encode matrices, instead of vectors.
#' @details
#' To remind you, ordinary RLE compresses highly repetitive vectors
#' by encoding the values in the vector as a series of
#' pairs \code{<val,len>}, where each pair represents a batch of
#' \emph{len} consecutive occurrences of value \emph{val}.
#' 
#' In 2D RLE, we assume that a matrix frequently has identical values
#' that occur in square-like patterns. We aim at finding maximal squares,
#' each of which contains only one value. We encode the matrix
#' as a series of pairs \code{<val,width>}, where each pair represents
#' a maximal square that has width \emph{width} and contains only \emph{val}.
#' 
#' The encoding starts from top-left (1,1) corner of the matrix and
#' proceeds column-by-column. Each time it encounters a matrix item that has
#' not been encoded yet, it treats this item as top-left item of a
#' square, and finds a maximal square whose all values are identical.
#' All included items are encoded by a single \code{<val,width>} pair.
#' Both encoding and decoding algorithms skip
#' the matrix items that have already been encoded/decoded.
#' @param mat matrix to be encoded
#' @param nr number of rows in decoded matrix
#' @param nc number of columns in decoded matrix
#' @param values values of squares in 2D RLE encoding
#' @param widths widths of squares in 2D RLE encoding
#' @examples
#' # Setup test matrix
#' testmat <- matrix(c(
#'   1,1,1,2,2,1,3,3,
#'   2,2,1,1,1,1,1,1,
#'   2,2,2,1,1,1,1,1,
#'   1,2,2,2,1,1,1,3,
#'   2,2,2,2,2,3,3,3,
#'   3,2,2,2,2,3,3,3,
#'   3,3,2,2,1,3,3,3), nrow=7, byrow=TRUE)
#' nr <- nrow(testmat)
#' nc <- ncol(testmat)
#' 
#' @seealso \code{\link{rle}}
NULL


#' \code{rle2d} encodes a matrix with 2D Run-Length Encoding.
#' @return \code{rle2d} returns a 2D Run-Length Encoding of given matrix.
#' It is a list that contains vectors:
#' \itemize{
#'   \item \code{values}: square values in the 2D encoding
#'   \item \code{widths}: square widths in the 2D encoding
#'   \item \code{max.widths}: maximum square widths that were possible
#'     (helps to compress the widths efficiently).
#' }
#' @examples
#' # Encode matrix with 2D RLE:
#' enc <- rle2d(testmat)
#' all(enc$values == c(1,2,1,2,3,3,1,2,3,1,1,2,2,2,1,2,2,1,2,2,1,1,1,1,3,3,1,3,1,1,3))
#' all(enc$widths == c(1,2,1,1,1,1,1,3,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,3,1,1,1,1,1,1))
#' 
#' @rdname rle2D
#' @export
rle2d <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  NEXT_COL <- nr + 1
  reconstructed <- matrix(NA, nrow=nr, ncol=nc)
  values <- integer(nr*nc)
  widths <- integer(nr*nc)
  max.widths <- integer(nr*nc)
  square.idx <- 1
  for (col in 1:nc) {
    row <- 1
    while (row < NEXT_COL) {
      if (is.na(reconstructed[row,col])) {
        v <- mat[row,col]
        wd <- 1
        # FIXME, check only "new" elements?
        while (row+wd<=nr && col+wd<=nc && all(mat[row:(row+wd), col:(col+wd)] == v)) { wd <- wd + 1 }
        reconstructed[row:(row+wd-1), col:(col+wd-1)] <- v
        values[square.idx] <- v
        widths[square.idx] <- wd
        max.widths[square.idx] <- min(nr-row+1, nc-col+1)
        square.idx <- square.idx + 1
        row <- row + wd
      } else {
        if (row==nr) { row <- NEXT_COL }
        else row <- row + head(which(is.na(c(reconstructed[(row+1):nr, col]))), 1)
        if (!length(row)) row <- NEXT_COL
      }
    }
  }
  stopifnot(all(reconstructed == mat))
  s <- square.idx - 1
  return(list(values=values[1:s], widths=widths[1:s], max.widths=max.widths[1:s]))
}


#' \code{invRle2d} decodes a 2D Run-Length Encoding into a matrix.
#' @return \code{invRle2d} returns a decoded matrix.
#' @rdname rle2D
#' @examples
#' # Decode matrix from a 2D RLE encoding:
#' m <- invRle2d(nr, nc, enc$values, enc$widths)
#' all(m == testmat)
#' 
#' @export
invRle2d <- function(nr, nc, values, widths) {
  NEXT_COL <- nr + 1
  mat <- matrix(NA, nrow=nr, ncol=nc)
  square.idx <- 1
  row <- 1
  for (col in 1:nc) {
    row <- 1
    while (row < NEXT_COL) {
      if (is.na(mat[row,col])) {
        v <- values[square.idx]
        wd <- widths[square.idx]
        max.width <- min(nr-row+1, nc-col+1)
        mat[row:(row+wd-1), col:(col+wd-1)] <- v
        square.idx <- square.idx + 1
        row <- row + wd
      } else {
        if (row==nr) row <- NEXT_COL
        else row <- row + head(which(is.na(c(mat[(row+1):nr, col]))), 1)
        if (!length(row)) row <- NEXT_COL
      }
    }
  }
  return(mat)
}

# Note that the squares may overwrite values (the value does not change, though).
# Therefore the following may not always hold:
#   sum(enc.1$widths^2) == length(testmat.1)
# Example:
#  a b b b
#  b b b b
#  b b b b
# Gives (val,width) pairs: (a,1), (b,2), (b,3)
# They encode writing 14 entries (=1^2+2^2+3^2), while matrix has 12 (=3*4) entries.


#' @rdname mdl
#' @name mdl
#' @aliases eliasDelta
#' @aliases mdlRle2d
#' @title MDL: Model selection with Minimum Description Length.
#' @description Methods for encoding and compressing data and model
#' descriptions so that an appropriate model is chosen by the MDL principle.
#' @details Minimum Description Length (MDL) is a principle that can be used
#' in model selection. In general, we should choose a model
#' that helps us to compress the data, but we should also favor simple
#' models over complicated ones. As the principle goes,
#' we should choose the model for which the sum of encoded data plus
#' encoded model description is minimum. The encoded versions should be
#' measured in bits.
#' @param n Positive integer to be encoded in bits.
#' @param mat Matrix to be encoded, contains \code{k} unique values.
#' @param k Number of unique values in \code{mat}.
#' @seealso MDL: \url{https://en.wikipedia.org/wiki/Minimum_description_length}
#' 
NULL


#' \code{eliasDelta} gives the length of encoding a positive integer
#' with \emph{Elias Delta coding}.
#' @return \code{eliasDelta} returns the length of integer encoding in bits.
#' @examples
#' # Encodings as 1:"1", 2:"0100", 3:"0101", 4:"01100", 15:"00100111", and 16:"001010000".
#' all(sapply(c(1,2,3,4,15,16), eliasDelta) == c(1,4,4,5,8,9))
#' 
#' @rdname mdl
#' @seealso Elias Delta coding:
#'   \url{https://en.wikipedia.org/wiki/Elias_delta_coding}
#' @export
eliasDelta <- function(n) floor(log2(n)) + 2*floor(log2(floor(log2(n))+1)) + 1


#' \code{mdlRle2d} gives the length of encoding a matrix that has
#'   \code{k} unique values.
#' @details
#' In \code{mdlRle2d} encoding, we use the following pieces of information:
#' \enumerate{
#'   \item number of rows in the matrix (Elias Delta coding), integer model description
#'   \item number of columns in the matrix (Elias Delta coding), integer model description
#'   \item number of unique values \code{k} (uniform from 1 to matrix length),
#'     integer model description
#'   \item matrix values from \code{k} possibilities
#'     (as \emph{values} and \emph{widths}
#'     from 2D RLE, pairs of positive integers for encoding data
#' }
#' @return \code{mdlRle2d} returns the encoding lengths of different components,
#' organized as a list with the following fields:
#' \itemize{
#'   \item \code{dims}: Bits for matrix dimensions.
#'   \item \code{range}: Bits for the range of possible (\code{k}) values.
#'   \item \code{values}: Bits for values of 2D RLE sub-squares within matrix.
#'   \item \code{widths}: Bits for widths of 2D RLE sub-squares within matrix.
#' }
#' The total encoding length of given \code{mat} and  \code{k} is
#' the sum of these lengths.
#' @examples
#' testmat <- matrix(c(
#'   1,1,1,2,2,1,3,3,
#'   2,2,1,1,1,1,1,1,
#'   2,2,2,1,1,1,1,1,
#'   1,2,2,2,1,1,1,3,
#'   2,2,2,2,2,3,3,3,
#'   3,2,2,2,2,3,3,3,
#'   3,3,2,2,1,3,3,3), nrow=7, byrow=TRUE)
#' k <- 3  # max(testmat)
#' mdl <- mdlRle2d(testmat, k)
#' abs(sum(unlist(mdl)) - 109.7863491) < 1e-6  # Total encoding length in bits.
#' 
#' # Reference encoding lengths computed as:
#' #   dims: eliasDelta(7) + eliasDelta(8) = 5 + 8
#' #   range: log2(7*8)
#' #   values: 31*log2(3)
#' #   widths: 9*0+4*1+6*log2(3)+4*2+3*log2(5)+3*log2(6)+2*log2(7), sum(log(max.widths)
#' # abs(unlist(mdl) - c(13, 5.80735492, 49.1338375, 41.84515664)) < 1e-6
#' 
#' @rdname mdl
#' @export
mdlRle2d <- function(mat, k) {
  mdl.dims <- eliasDelta(nrow(mat)) + eliasDelta(ncol(mat))
  mdl.range <- log2(length(mat))  # encoding the range from 1 to k
  encoding <- rle2d(mat)
  s <- length(encoding$values)
  # No need to encode the number of squares: just decode until matrix is filled.
  mdl.square.values <- s * log2(k)
  mdl.square.widths <- sum(log2(encoding$max.widths))
  return(list(dims=mdl.dims, range=mdl.range, values=mdl.square.values, widths=mdl.square.widths))
}


# Testing
#########

test <- function() {
  testmat.1 <- matrix(c(
    1,1,1,1,6,3,3,3,4,4,3,3,4,4,
    1,1,1,1,1,6,3,5,4,4,4,4,4,4,
    1,1,1,1,1,5,6,3,5,4,4,4,4,3,
    1,1,1,1,1,5,5,6,3,5,4,4,3,3,
    2,2,2,2,6,5,7,5,6,3,5,3,3,3,
    3,2,2,2,2,7,3,3,1,6,3,5,5,5,
    3,3,2,2,4,4,3,3,8,1,6,5,5,5,
    3,3,3,2,4,1,1,1,1,1,1,1,1,1
  ), nrow=8, byrow=TRUE)
  enc.1 <- rle2d(testmat.1)
  stopifnot(enc.1$values == c(1,2,3,3,2,2,3,2,2,2,2,6,1,1,1,6,2,4,4,3,6,5,5,5,7,4,1,3,3,6,5,7,3,1,3,5,3,6,5,1,4,5,3,6,1,8,1,4,5,3,6,1,1,3,4,4,5,3,6,1,3,4,3,5,1,4,4,3,1,3,5,5,1))
  stopifnot(enc.1$widths == c(4,1,1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,2,1,2,1,2,1,1,1,1,1))
  stopifnot(head(enc.1$max.widths, 12) == c(8,4,3,2,4,2,1,4,3,2,1,8))
  stopifnot(all(testmat.1 == invRle2d(nrow(testmat.1), ncol(testmat.1), enc.1$values, enc.1$widths)))

  testmat.2 <- matrix(c(
    1,1,1,2,2,1,3,3,
    2,2,1,1,1,1,1,1,
    2,2,2,1,1,1,1,1,
    1,2,2,2,1,1,1,3,
    2,2,2,2,2,3,3,3,
    3,2,2,2,2,3,3,3,
    3,3,2,2,1,3,3,3
  ), nrow=7, byrow=TRUE)
  enc.2 <- rle2d(testmat.2)
  stopifnot(enc.2$values == c(1,2,1,2,3,3,1,2,3,1,1,2,2,2,1,2,2,1,2,2,1,1,1,1,3,3,1,3,1,1,3))
  stopifnot(enc.2$widths == c(1,2,1,1,1,1,1,3,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,3,1,1,1,1,1,1))
  stopifnot(head(enc.2$max.widths, 14) == c(7,6,4,3,2,1,7,4,1,6,6,5,1,5))
  stopifnot(all(testmat.2 == invRle2d(nrow(testmat.2), ncol(testmat.2), enc.2$values, enc.2$widths)))
  stopifnot(abs(sum(unlist(mdlRle2d(testmat.2, length(unique(c(testmat.2))))) - c(5 + 8, log2(7*8), 31*log2(3), sum(log2(enc.2$max.widths))))) < 1e-9)
  
  testmat.rnd <- matrix(sample(2, 600*601, replace=TRUE), nrow=600)
  enc.rnd <- rle2d(testmat.rnd)
  stopifnot(all(testmat.rnd == invRle2d(nrow(testmat.rnd), ncol(testmat.rnd), enc.rnd$values, enc.rnd$widths)))

  testmat <- testmat.rnd
  Rprof(filename="Rprof.out", append=FALSE, interval=0.01, line.profiling=TRUE)
  enc <- rle2d(testmat)
  reconstructed <- invRle2d(nrow(testmat), ncol(testmat), enc$values, enc$widths)
  Rprof(NULL)
  print(summaryRprof('Rprof.out', lines='show'))
}
#test()

