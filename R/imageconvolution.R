##
## Image and photo 2D convolution
##
## Esa Junttila, 2016-03-23


#' Extend image dimensions by existing values.
#' 
#' Extend the row and column dimensions of an image array by given amount.
#' All the available methods use the existing values in the array to fill
#' in the values. All channels in the array are affected the same way.
#' The method for choosing the fill-values depends on the \emph{padding}
#' method (or value) chosen.
#' @details This function is needed to expand an array so that
#'   2D filtering (convolution) works on the edge areas of the image.
#' @param img An image construct with \code{m x n x 3} dimensions:
#'   rows, columns, and channels.
#' @param amounts Vector of length 2 that specifies both the number
#'   of padded rows (top and down) and the number of padded columns
#'   (left and right).
#' @param pad Name of a method for expanding the row and column dimensions of 
#' \code{img}, or a numerical value to fill. The selection of \code{pad}
#' method has an impact on results on edge areas. Available values are:
#' \itemize{
#'   \item if numerical value, then uses that value.
#'   \item \code{'replicate'}: use the edge value (same as
#'     'nearest neighbor' or 'reflect').
#'   \item \code{'circular'}: use the values from the other side of the
#'     matrix as if the matrix was circular. \strong{NOT IMPLEMENTED YET.}
#'   \item \code{'symmetric'}: use the values near the edge like
#'     a symmetric mirror. \strong{NOT IMPLEMENTED YET.}
#' }
#' @examples
#' img <- array(1:36, dim=c(3,4,3))
#' all(img[ , , 1] == 1:12)
#' padded <- padarray(img, c(1,2), 'replicate')
#' 
#' first.channel <- matrix(c(
#'   1,1, 1,4,7,10, 10,10,
#'   
#'   1,1, 1,4,7,10, 10,10,
#'   2,2, 2,5,8,11, 11,11,
#'   3,3, 3,6,9,12, 12,12,
#'   
#'   3,3, 3,6,9,12, 12,12), nrow=5, byrow=TRUE)
#' 
#' all(first.channel == padded[ , , 1])
#' 
#' pad.by.val <- padarray(img, c(1,2), -8)
#' ref.channel <- matrix(c(
#'   -8,-8, -8,-8,-8,-8, -8,-8,
#'   
#'   -8,-8, 13,16,19,22, -8,-8,
#'   -8,-8, 14,17,20,23, -8,-8,
#'   -8,-8, 15,18,21,24, -8,-8,
#'   
#'   -8,-8, -8,-8,-8,-8, -8,-8), nrow=5, byrow=TRUE)
#' all(ref.channel == pad.by.val[ , , 2])
#' dim(pad.by.val) == dim(padded) && all(pad.by.val[pad.by.val > 0] == 1:36)
#' 
#' @export
padarray <- function(img, amounts, pad) {
  # Original array variables:
  dims <- dim(img)
  h <- dims[1]
  w <- dims[2]
  nchannels <- dims[3]
  dh <- amounts[1]  # num padded rows (added both sides)
  dw <- amounts[2]  # num padded cols (added both sides)
  # Padded array variables:
  ph <- h + 2*dh  # padded array height
  pw <- w + 2*dw  # padded array width
  r.in <- (dh+1):(dh+h)   # original rows in padded array
  r.up <- 1:dh            # up-padded rows
  r.down <- (ph-dh+1):ph
  c.in <- (dw+1):(dw+w)   # original columns
  c.left <- 1:dw          # left-padded columns
  c.right <- (pw-dw+1):pw
  # Do the work:
  if (is.numeric(pad)) {
    padded <- array(rep(pad, ph*pw*nchannels), dim=c(ph,pw,nchannels))  # fill vals
    padded[r.in, c.in, ] <- img[1:h, 1:w, ]  # overwrite original values
  } else if (is.character(pad)) {
    padded <- array(dim=c(ph,pw,nchannels))  # fill with NA's
    padded[r.in, c.in, ] <- img[1:h, 1:w, ]  # overwrite original values
    pad.method <- tolower(pad)
    if (pad.method == 'replicate') {  # same as 'nearest neighbor' or 'reflect'
      for (ch in 1:nchannels) {
        if (dh) padded[r.up  , c.in   , ch] <- matrix(rep(img[1  , 1:w, ch], dh), nrow=dh, byrow=TRUE)   # up
        if (dh) padded[r.down, c.in   , ch] <- matrix(rep(img[h  , 1:w, ch], dh), nrow=dh, byrow=TRUE)   # down
        if (dw) padded[r.in  , c.left , ch] <- matrix(rep(img[1:h, 1  , ch], dw), nrow=dw, byrow=FALSE)  # left
        if (dw) padded[r.in  , c.right, ch] <- matrix(rep(img[1:h, w  , ch], dw), nrow=dw, byrow=FALSE)  # right
        if (dh && dw) {
          padded[r.up  , c.left , ch] <- img[1, 1, ch]  # up-left corner
          padded[r.up  , c.right, ch] <- img[1, w, ch]  # up-right corner
          padded[r.down, c.left , ch] <- img[h, 1, ch]  # down-left corner
          padded[r.down, c.right, ch] <- img[h, w, ch]  # down-right corner
        }
      }
    }
    else if (pad.method == 'symmetric') {
      stop('pad method not implemented yet')  # how to handle corners?
    }
    else if (pad.method == 'circular') {
      stop('pad method not implemented yet')  # how to handle corners?
    }
    else stop(paste('unsupported pad method:', pad.method))
  }
  else stop(paste('unsupported pad type:', typeof(pad)))
  return(padded)
}


#' @rdname convolution2D
#' @name convolution2D
#' @aliases conv2dFFT
#' @aliases conv2
#' @title Convolution 2D operator on matrices and images.
#' @description Perform 2D convolution on target matrix with
#'   given \emph{filter} matrix. The different implementations produce (almost)
#'   identical results, except for small rounding and replication errors.
#' @param m matrix to apply 2D convolution on.
#' @param img image on whose channels (matrices) 2D convolution is applied.
#'   It is an array with three dimensions as \emph{rows}, \emph{columns},
#'   and \emph{channels}.
#' @param f filter matrix used in the 2D convolution operations, the dimensions
#'   of which must be odd and smaller than those of matrix \code{m}.
#' @param dims required dimensions of the output.
#' @param shape name for the shape of required output dimensions:
#'   \itemize{
#'     \item \code{'same'} (default): use dimensions from \code{m},
#'       returning the "original" part of the padded (convoluted) matrix.
#'     \item \code{'full'}: return the full padded (convoluted) matrix.
#'     \item \code{'valid'}: return the inner part of the convoluted
#'       matrix; that is, the values not affected by padding the border.
#'   }
#' @param pad parameter for function \code{\link{padarray}}: what value
#'   or method is used for padding the matrix when enlargened during
#'   2D convolution? (default 0)
#' @param fft choice of implementation: [boolean] should we use Fast Fourier
#'   Transform (FFT) for 2D convolution or direct matrix operations?
#' @details Extension of ordinary (1D) convolution method on sequences,
#'   like \code{fft}. The 2D extension of FFT returns the real-number
#'   parts of the complex numbers.
#'   
#'   A collection of pre-formatted and commonly used filter matrices
#'   can be found in  \code{convolutionfilter}.
#' @seealso
#'   \url{http://www.songho.ca/dsp/convolution/convolution2d_example.html}
#'   
#'   \url{http://matlabtricks.com/post-26/tutorial-on-2d-convolution-of-images}
#'   
#'   \url{https://en.wikipedia.org/wiki/Multidimensional_discrete_convolution}
#'   
#'   \code{\link{fft}}
#'   
#'   \code{\link{convolutionfilter}}
#' @examples
#' img <- array(matrix(1:9, nrow=3, byrow=TRUE), dim=c(3,3,1))
#' m <- img[,,1, drop=TRUE]
#' f <- matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3, byrow=TRUE)
#' p <- padarray(img, c(1,1), 0)
#' pm <- p[,,1, drop=TRUE]
#' target.vals <- c(-13,-20,-17, -18,-24,-18, 13,20,17)
#' target <- array(matrix(target.vals, byrow=TRUE, nrow=3), dim=c(3,3,1))
#' tol <- 1e-13
NULL


#' \code{conv2dFFT} computes the 2D convolution, using Fast Fourier Transform,
#'   on matrix \code{m} with filter matrix \code{f}.
#' @return \code{conv2dFFT} returns a convoluted matrix that has
#'   dimensions as specified in \code{shape} parameter.
#' @examples
#' all(abs(conv2dFFT(pm, f, dims=dim(m)) - target[,,1, drop=TRUE]) < tol)
#' @rdname convolution2D
#' @export
conv2dFFT <- function(m, f, dims) {
  nr <- dim(m)[1]
  nc <- dim(m)[2]
  fr <- dim(f)[1]
  fc <- dim(f)[2]
  if (nr<fr || nc<fc) stop("filter matrix 'f' cannot be larger than matrix 'm' (at the moment)")
  if (!(fr%%2) || !(fc%%2)) stop('cannot support even filter dimension lengths at the moment')
  ar <- floor(fr/2)
  ac <- floor(fc/2)
  mat <- array(0, dim=dim(m))  # make filter matrix 'f' as large as 'm'
  mat[(-ar:ar)%%nr+1, (-ac:ac)%%nc+1] <- f  # center point of filter matrix is at mat(1,1)
  y <- fft(fft(m) * fft(mat), inverse=TRUE) / (nr*nc)
  if (!all(dim(y) == dims)) {
    y <- array(y[(ar+1):(nr-ar), (ac+1):(nc-ac)], dim=dims)
  }
  return(Re(y))
}


#' \code{conv2} computes the 2D convolution, using a naive matrix
#'   implementation, on matrix \code{m} with filter matrix \code{f}.
#' @return \code{conv2} returns a convoluted matrix that has
#'   dimensions as specified in \code{dims} parameter.
#' @examples
#' all(abs(conv2(pm, f) - target[,,1, drop=TRUE]) < tol)
#' all(conv2(matrix(c(1,5,2,3,8,7,3,6,3,3,9,1),byrow=TRUE,nrow=3), matrix(c(1,2,3,0,0,0,6,5,4),nrow=3,byrow=TRUE)) == c(65,76))
#' @rdname convolution2D
#' @export
conv2 <- function(m, f, shape='same') {
  nr=nrow(m)  #number of rows
  nc=ncol(m)
  fr=nrow(f)  # number of filter rows
  fc=ncol(f)
  if (!(fr%%2) || !(fc%%2)) stop('No support for even dimensions at the moment.')
  ar <- floor(fr/2)  # added border rows (both sides)
  ac <- floor(fc/2)  # added border columns (both sides)
  start.r <- ar + 1  # first row index of the center cell
  start.c <- ac + 1  # first column index of the center cell
  end.r <- nr - ar   # last center
  end.c <- nc - ac   # last center
  # Using two for loops is of course slow beyond belief.
  filtr <- f[fr:1, fc:1, drop=FALSE]  # flipped and mirrored kernel
  res <- array(0, dim=dim(m))  # initialized with zeros
  res.rows <- (1+ar):(nr-ar)  # these rows are affected
  res.cols <- (1+ac):(nc-ac)  # these columns are affected
  for (i in 1:fr) {  # row indices from filter matrix (smaller matrix)
    for (j in 1:fc) {  # column indices from filter matrix
      src.rows <- i:(nr-fr+i)
      src.cols <- j:(nc-fc+j)
      res[res.rows, res.cols] <- res[res.rows, res.cols] + filtr[i,j] * m[src.rows, src.cols]
    }
  }
  # Shape output:
  if (shape == 'same') {
    rows <- start.r:end.r
    cols <- start.c:end.c
  }
  else if (shape == 'full')  {  # append missing border values
    left <- 1:ac           # leftmost columns indices
    right <- (nc-ac+1):nc  # rightmost columns
    up <- 1:ar             # topmost rows
    down <- (nr-ar+1):nr   # bottommost rows
    if (ar) res[up, res.cols] <- m[up, res.cols]        # up
    if (ar) res[down, res.cols] <- m[down, res.cols]    # down
    if (ac) res[res.rows, left] <- m[res.rows, left]    # left
    if (ac) res[res.rows, right] <- m[res.rows, right]  # right
    if (ar && ac) {
      res[up, left] <- m[up, left]        # up-left corner
      res[up, right] <- m[up, right]      # up-right corner
      res[down, left] <- m[down, left]    # down-left corner
      res[down, right] <- m[down, right]  # down-right corner
    }
    rows <- 1:nr
    cols <- 1:nc
  }
  else if (shape == 'valid') {
    rows <- (2*ar+1):(nr-2*ar)  # drop part of computed values
    cols <- (2*ac+1):(nc-2*ac)
  }
  else stop(paste('Unknown shape:', shape))
  return(res[rows, cols])
  #return(res[start.r:end.r, start.c:end.c])
}

# Slower alternative:
# FIXME: computations of 2D convolution should be fast instead of super slow
#flipped <- c(f[fr:1, fc:1])  # flipped kernel
#res <- array(dim=dim(m))
#for (i in start.r:end.r) {
#  for (j in start.c:end.c) {  # Using two for loops is of course slow beyond belief.
#    res[i,j] <- sum(flipped * c(m[(i-ar):(i+ar), (j-ac):(j+ac)]))
#  }
#}


#' \code{imfilter} computes the image filtering (2D convolution),
#'   on each channel of an image \code{img} with filter matrix \code{f}.
#' @return \code{imfilter} returns a filtered image that has
#'   2D convoluted matrices (channels).
#' @details \code{imfilter} is a simplified
#'   version of \url{http://www.makalab.org/ojb/imfilter.m}.
#' @examples
#' all(abs(imfilter(img, f, pad=0) - target) < tol)
#' # Example from (http://www.songho.ca/dsp/convolution/convolution2d_example.html):
#' #
#' #         | 1 2 3 |        | -1 -2 -1 |             | -13 -20 -17 |
#' #  img := | 4 5 6 |,  f := |  0  0  0 |,  result := | -18 -24 -18 |
#' #         | 7 8 9 |        |  1  2  1 |             |  13  20  17 |
#' #
#' @rdname convolution2D
#' @export
imfilter <- function(img, f, pad=0, fft=FALSE) {
  dims <- dim(img)
  h <- dims[1]
  w <- dims[2]
  arr <- if (is.matrix(img)) array(img, dim=c(h,w,1)) else img
  nchannels <- dim(arr)[3]
  m <- nrow(f)
  n <- ncol(f)
  if (m%%2==0 || n%%2==0) stop('filter dimensions mod(2) not supported at the moment in imfilter')
  src = padarray(arr, floor(c(m/2, n/2)), pad);
  res <- array(dim=c(h,w,nchannels))
  for (ch in 1:nchannels) {
    if (fft) {  # Fast Fourier Transform (2D convolution)
      mat <- conv2dFFT(src[ , , ch, drop=TRUE], f, dim(res[,,ch]))
    } else {  # direct convolution
      mat <- conv2(src[ , , ch, drop=TRUE], f)
      tol <- 1e-12
      mat[abs(mat) < tol] <- 0      # suppress floating-point inaccuracies
      mat[abs(mat - 1) < tol] <- 1  # same here
    }
    res[ , , ch] <- mat
  }
  return(res)
}

