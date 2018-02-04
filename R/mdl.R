# Image MDL: Minimum Description Length
# WORK IN PROGRESS


# RLE
#####
# ?rle

# Quadtree
##########

# 2D RLE: 2D Run-Length Encoding with maximal squares, column-by-column
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


# MDL: Minimum Description Length
#####

eliasDelta <- function(n) floor(log2(n)) + 2*floor(log2(floor(log2(n))+1)) + 1  # Elias delta encoding, bits

mdlRle2d <- function(mat, k) {
  mdl.dims <- eliasDelta(nrow(mat)) + eliasDelta(ncol(mat))
  mdl.range <- log2(length(mat))  # encoding the range from 1 to k
  encoding <- rle2d(mat)
  s <- length(encoding$values)
  # There is no need to encode the number of squares: just decode until matrix is filled.
  mdl.square.values <- s * log2(k)
  mdl.square.widths <- sum(log2(encoding$max.widths))
  return(list(dims=mdl.dims, range=mdl.range, values=mdl.square.values, widths=mdl.square.widths))
}




# Improved 2D RLE: "An Improved Two-Dimensional Run-Length Encoding Scheme and Its Application"
#################


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

