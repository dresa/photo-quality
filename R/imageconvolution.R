##
## Image and photo 2D convolution
##
## Esa Junttila, 2016-03-23


# Extend arrays by given amounts, with values generating according to 'pad'
# pad: numeric value, 'circular', 'replicate', or 'symmetric'
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


# 2D convolution with Fast Fourier Transform.
# note
conv2dFFT <- function(x, m, dims) {
  nr <- dim(x)[1]
  nc <- dim(x)[2]
  mr <- dim(m)[1]
  mc <- dim(m)[2]
  if (nr<mr || nc<mc) stop("filter matrix 'm' cannot be larger than matrix 'x' (at the moment)")
  if (!(mr%%2) || !(mc%%2)) stop('cannot support even filter dimension lengths at the moment')
  ar <- floor(mr/2)
  ac <- floor(mc/2)
  mat <- array(0, dim=dim(x))  # make filter matrix 'm' as large as 'x'
  mat[(-ar:ar)%%nr+1, (-ac:ac)%%nc+1] <- m  # center point of filter matrix is at mat(1,1)
  y <- fft(fft(x) * fft(mat), inverse=TRUE) / (nr*nc)
  if (!all(dim(y) == dims)) {
    y <- array(y[(ar+1):(nr-ar), (ac+1):(nc-ac)], dim=dims)
  }
  return(Re(y))
}


# A naive and slow 2D convolution
# Matrix 'mat' and filter 'f' are both matrices, but f is assumed to be a smaller one.
conv2 <- function(mat, f, shape='same') {
  nr=nrow(mat)  #number of rows
  nc=ncol(mat)
  dr=nrow(f)  # number of filter rows
  dc=ncol(f)
  if (!(dr%%2) || !(dc%%2)) stop('No support for even dimensions at the moment.')
  ar <- floor(dr/2)  # added border rows (both sides)
  ac <- floor(dc/2)  # added border columns (both sides)
  start.r <- ar + 1  # first row index of the center cell
  start.c <- ac + 1  # first column index of the center cell
  end.r <- nr - ar   # last center
  end.c <- nc - ac   # last center
  # Using two for loops is of course slow beyond belief.
  filtr <- f[dr:1, dc:1, drop=FALSE]  # flipped and mirrored kernel
  res <- array(0, dim=dim(mat))  # initialized with zeros
  res.rows <- (1+ar):(nr-ar)  # these rows are affected
  res.cols <- (1+ac):(nc-ac)  # these columns are affected
  for (i in 1:dr) {  # row indices from filter matrix (smaller matrix)
    for (j in 1:dc) {  # column indices from filter matrix
      src.rows <- i:(nr-dr+i)
      src.cols <- j:(nc-dc+j)
      res[res.rows, res.cols] <- res[res.rows, res.cols] + filtr[i,j] * mat[src.rows, src.cols]
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
    if (ar) res[up, res.cols] <- mat[up, res.cols]        # up
    if (ar) res[down, res.cols] <- mat[down, res.cols]    # down
    if (ac) res[res.rows, left] <- mat[res.rows, left]    # left
    if (ac) res[res.rows, right] <- mat[res.rows, right]  # right
    if (ar && ac) {
      res[up, left] <- mat[up, left]        # up-left corner
      res[up, right] <- mat[up, right]      # up-right corner
      res[down, left] <- mat[down, left]    # down-left corner
      res[down, right] <- mat[down, right]  # down-right corner
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
#flipped <- c(f[dr:1, dc:1])  # flipped kernel
#res <- array(dim=dim(mat))
#for (i in start.r:end.r) {
#  for (j in start.c:end.c) {  # Using two for loops is of course slow beyond belief.
#    res[i,j] <- sum(flipped * c(mat[(i-ar):(i+ar), (j-ac):(j+ac)]))
#  }
#}


# Image filter
# Simplified from: http://www.makalab.org/ojb/imfilter.m
# Example (http://www.songho.ca/dsp/convolution/convolution2d_example.html):
#   imfilter(array(matrix(1:9, nrow=3, byrow=TRUE), dim=c(3,3,1)), matrix(c(-1,-2,-1,0,0,0,1,2,1), nrow=3, byrow=TRUE), pad=0)
#     , , 1
#          [,1] [,2] [,3]
#     [1,]  -13  -20  -17
#     [2,]  -18  -24  -18
#     [3,]   13   20   17
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
    if (fft) {
      mat <- conv2dFFT(src[ , , ch, drop=TRUE], f, dim(res[,,ch]))  # Fast Fourier Transform (2D convolution)
    } else {
      mat <- conv2(src[ , , ch, drop=TRUE], f)  # direct convolution
      tol <- 1e-12
      mat[abs(mat) < tol] <- 0      # suppress small inaccuracies in floating point calculations
      mat[abs(mat - 1) < tol] <- 1  # same here
    }
    res[ , , ch] <- mat
  }
  return(res)
}

