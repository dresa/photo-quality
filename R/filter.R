#
# Quality measures for photographs: no-reference analysis of aesthetic features
# Author: Esa Junttila
# 

library(jpeg)  # function 'readJPG' to read JPEG files
library(png)   # function 'readPNG' to read PNG files
library(grid)  # function 'grid.raster' to view images
library(circular)
library(ggplot2)  # hue-usage in polar coordinates


###########
## IMAGE ##
###########

# Color channel identifiers:
RED <- list(idx=1, name='RED')
GREEN <- list(idx=2, name='GREEN')
BLUE <- list(idx=3, name='BLUE')
RGB <- list(RED, GREEN, BLUE)

HUE <- list(idx=1, name='HUE')
SATURATION <- list(idx=2, name='SATURATION')
VALUE <- list(idx=3, name='VALUE')
HSV <- list(HUE, SATURATION, VALUE)

# Read image file (either a valid PNG or JPEG extension required)
readImage <- function(filename) {
  refname <- tolower(trimws(filename))
  # Does text end with any of given suffixes? Sensitive to regex chars!
  endsWithAny <- function(txt, suffixes) {
    any(sapply(suffixes, function(s) grepl(paste(s,'$',sep=''), txt)))
  }
  if (endsWithAny(refname, c('\\.jpg', '\\.jpeg'))) {  # . is special regex char
    img <- readJPEG(filename)
  } else if (endsWithAny(refname, c('\\.png'))) {
    img <- readPNG(filename)
  } else stop('Image file must be a .PNG or .JPG file')
  dimnames(img) <- list(NULL, NULL, c('Red', 'Green', 'Blue'))
  return(img)
}

# Return a matrix that represents channel 'col.channel' of image 'img'.
extractRGBChannel <- function(img.rgb, col.channel) img.rgb[ , , col.channel$idx]

extractHSVChannel <- function(img.hsv, hsv.channel) img.hsv[ , , hsv.channel$idx]

# Convert RGB image to an HSV image. Assumption: RGB values are between 0 and 1.
# Hues are degrees between 0 and 360; Saturations and Values are between 0 and 1.
toHSV <- function(img.rgb, radians=FALSE) {
  # Extract R,G, and B color channels
  red <- extractRGBChannel(img.rgb, RED)
  green <- extractRGBChannel(img.rgb, GREEN)
  blue <- extractRGBChannel(img.rgb, BLUE)
  # Temporary variables
  dims <- dim(red)
  max.values <- pmax(red, green, blue)
  diffs <- max.values - pmin(red, green, blue)
  # Compute H(ue), S(saturation), and V(alue):
  # S(aturation):
  saturations <- diffs / max.values
  saturations[max.values == 0] <- 0  # replace Inf by zero (originates from division by zero)
  # H(ue):
  R<-0; G<-2; B<-4  # additive color-shift values in the hue formula
  max.layer <- array(R, dim=dims)  # which color layer has maximum value? R, G, or B?
  max.layer[green > red] <- G
  max.layer[blue > pmax(red, green)] <- B
  red.mask <- max.layer == R
  green.mask <- max.layer == G
  blue.mask <- max.layer == B
  hues <- array(dim=dims)
  angle <- if (radians) pi/3 else 60
  hues[red.mask]   <- angle * ((green[red.mask]  - blue[red.mask])   / diffs[red.mask]   + R)%%6
  hues[green.mask] <- angle * ((blue[green.mask] - red[green.mask])  / diffs[green.mask] + G)%%6
  hues[blue.mask]  <- angle * ((red[blue.mask]   - green[blue.mask]) / diffs[blue.mask]  + B)%%6
  hues[max.values==0] <- NA  # undefined hue is not needed; red (zero degrees) as default
  hues[diffs==0] <- NA  # undefined hue is not needed; red (zero degrees) as default
  # V(alue):
  values <- max.values
  # Return an HSV image as an (M x N x 3) array
  dim.names <- list(NULL, NULL, c('Hue', 'Saturation', 'Value'))
  return(array(c(hues, saturations, values), dim=c(dims[1], dims[2], 3), dimnames=dim.names))
}


############
## COMMON ##
############
limit <- function(x, low, high) pmax(low, pmin(high, x))

binarySearch <- function(func, target, low, high, tol=1e-12) {
  if ( high < low ) { return(NULL) }
  else {
    mid <- (low + high) / 2
    func.val <- func(mid)
    if (abs(func.val - target) < tol) { return(mid) }
    else if (func.val > target) { binarySearch(func, target, low, mid) }
    else if (func.val < target) { binarySearch(func, target, mid, high) }
  }
}

# Misunderstood an article. Obsolete.
# Bessel function of the first kind.
# Approximately correct. It's hard to deduce a good number of terms, now 41.
# Highly restricted search area. The function is monotonic at least in range (-1.84, 1.84).
#besselJ <- function(alpha, x) sum(sapply(0:40, function(k) (-x^2/4)^k/factorial(k)/gamma(alpha+k+1)*(x/2)^alpha))
#invertedBesselJ <- function(alpha, y, low=-1.84, high=1.84) binarySearch(function(x) besselJ(alpha, x), y, low, high, tol=1e-6)

#Evaluates the first and zeroth order Bessel functions of the first kind at a specified non-negative real
#number, and returns the ratio. [Related to circularity]
#See: https://cran.r-project.org/web/packages/circular/circular.pdf
#https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1.R?view=markup&root=circular
A1 <- function(kappa) {
    result <- besselI(kappa, nu=1, expon.scaled = TRUE)/besselI(kappa, nu=0, expon.scaled = TRUE)
    return(result)
}
#Inverse function of the ratio of the first and zeroth order Bessel functions of the first kind.  This
#function is used to compute the maximum likelihood estimate of the concentration parameter of a
#von Mises distribution.
# https://r-forge.r-project.org/scm/viewvc.php/pkg/R/A1inv.R?view=markup&root=circular&pathrev=6
A1inv <- function(x) {
   ifelse (0 <= x & x < 0.53, 2 * x + x^3 + (5 * x^5)/6,
           ifelse (x < 0.85, -0.4 + 1.39 * x + 0.43/(1 - x), 1/(x^3 - 4 * x^2 + 3 * x)))
}

radianToDegree <- function(rad) rad*(360/(2*pi))
degreeToRadian <- function(degree) degree*(2*pi/360)


#############
## VIEWING ##
#############

# Open a new window that shows given RGB image 'img'.
#   <img>: a single color channel (matrix with dim (h,w)) OR
#          an RGB image with red, green and blue channels (array with dim (h,w,3))
#  <size>: size of the minimum window axis (height or width), for example 4
# <title>: string in the window title
# Examples:
#   view(matrix(1/(1:24), nrow=6))
#   view(array(runif(1/(1:72)), dim=c(6,4,3)))
view <- function(img, size=4, title='Image') {
  dims <- dim(img)
  img.h <- dims[1]
  img.w <- dims[2]
  win.h <- if (img.h <= img.w) size else size * img.h/img.w
  win.w <- if (img.w <= img.h) size else size * img.w/img.h
  dev.new(width=win.w, height=win.h, title=title)
  grid.raster(img, interpolate=FALSE)
}


rescaleImage <- function(img) {
  # Rescale values: set 0 as minimum value, and 1 as maximum value. Preserve dimensions.
  a <- min(img)
  b <- max(img)
  return(array(pmax(0, pmin(1, (img-a)/(b-a))), dim=dim(img)))
}

# View (or save) a plot of color hues in a polar histogram (DDC ~ Dispersion Dominant Colors).
# This is very slow for an unknown reason. You need to wait half a minute for rendering.
viewDDC <- function(filename, img.hsv, mu, kappa) {
  hues <- degreeToRadian(extractHSVChannel(img.hsv, HUE))  # convert from degrees to radians
  num.bins <- 120
  x.max <- 2*pi  # 2*pi radians is 360 degrees
  h.raw <- c(hues)
  clean.hues <- h.raw[!is.na(h.raw)]
  bin.factor <- x.max / num.bins
  d <- (clean.hues %/% bin.factor)*bin.factor
  cols <- hsv(h=(d%%x.max)/x.max)
  toCirc <- function(x) circular(x, type='angles', units='radians', rotation='clock')
  size.factor <- max(table(findInterval(clean.hues, seq(0, 2*pi, length.out=num.bins+1))))
  tick.labels <- c('0', expression(pi * '/2'), expression(pi), expression('3' * pi * '/2'))
  p <- qplot(d, fill=cols, xlim=c(0, x.max), bins=num.bins, xlab='hue (radians)') +
    scale_fill_identity() +
    coord_polar(start=pi/2) +
    scale_x_continuous(breaks=c(0, pi/2, pi, 3*pi/2), labels=tick.labels) +
    stat_function(fun=function(x) size.factor*dvonmises(toCirc(x), mu=toCirc(mu), kappa=kappa), colour='black') +
    theme(axis.text.x=element_text(face="plain", color='black', size=16, angle=0),
          axis.text.y=element_text(face="plain", color='black', size=12, angle=0))
  print(p)
  if (!is.null(filename)) ggsave(filename, width=12, height=12, units='cm')
}


#############
## FILTERS ##
#############

meanFilter <- function(size) {
  if (!(size%%2)) stop('even filter size not possible at the moment')
  f <- matrix(1/(size*size), nrow=size, ncol=size)
  return(f)
}

gaussianFilter1D <- function(n, sigma, as.row=FALSE) {
  variance <- sigma^2
  # Due to normalization, we can omit constants from Gaussian:
  #   1/(sqrt(2*pi)*sigma)*exp(-x^2/(2*sigma^2))
  # This method simply uses the midpoints of the "pixels", which neglects non-linear shapes.
  # Points: -n, -(n-1), ..., 0, 1, 2, ..., n
  pos.indices <- 1:n  # symmetric (and g(0)=1), so only positive points need to be computed
  vals <- exp(-pos.indices^2/(2*variance))  # note that function value at 0 is 1
  m <- matrix(c(rev(vals), 1, vals) / (2*sum(vals)+1))  # reversed, zero-point, values: normalize such that sum=1.0
  return(if (as.row) t(m) else m)
}

# For testing purposes. A more efficient method of applying two 1D filters exists.
gaussianFilter2D <- function(n, sigma) {
  variance <- sigma^2
  # Due to normalization, we can omit constants from Gaussian:
  #   1/(2*pi*sigma^2)*exp(-(x^2+y^2)/(2*sigma^2))
  # This method simply uses the midpoints of the "pixels", which neglects non-linear shapes.
  indices <- -n:n
  idx.grid <- expand.grid(indices, indices)
  x <- idx.grid$Var1
  y <- idx.grid$Var2
  vals <- exp(-(x^2+y^2)/(2*variance))
  vals <- vals / sum(vals)  # normalize such that sum=1.0
  m <- matrix(vals, nrow=length(indices))
  return(m)
}

sobelVertical <- function() {
  return(matrix(c(
    c(-1, 0, +1),
    c(-2, 0, +2),
    c(-1, 0, +1)
  ), nrow=3, byrow=TRUE))
}

#################
## CONVOLUTION ##
#################

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


####################
## DERIVED IMAGES ##
####################

# Nobuyuki Otsu, ?A threshold selection method from gray-level histograms".
# IEEE Trans. Sys., Man., Cyber. 9 (1): 62?66. (1979). 
# Ported from http://clickdamage.com/sourcecode/code/otsuThreshold.m
otsuThreshold <- function(histogram) {
  n <- length(histogram)
  idxs <- 1:n
  cumus <- cumsum(histogram)
  weighted <- cumsum(histogram*idxs)
  u <- weighted / cumus
  tmp <- cumus[n] - cumus
  v <- (weighted[n] - weighted) / (tmp + (tmp==0))
  f <- 1.0 * cumus * (cumus[n] - cumus) * (u-v)^2
  i <- which.max(f)
  return(i)
}


# Simplified blur image (might be incorrectly understood)
meanBlurred <- function(img, size=3) {
  f <- meanFilter(size)
  blurred <- imfilter(img, f, 'replicate')
  return(blurred)
}

gaussianBlurred <- function(img, sigma) {
  # Gaussian blur is separable: combining vertical and horizontal 1D convolutions
  # gives the same outcome as using a single 2D convolution.
  n <- as.integer(ceiling(3*sigma))  # n pixels to left, up, right, and down
  #gf2D <- gaussianFilter2D(n, sigma)  # Slower direct 2D method for testing
  #blurred <- imfilter(img, gf2D, 'replicate')
  #return(blurred)
  gf1D <- gaussianFilter1D(n, sigma)                # (vertical) 1D Gaussian blur
  img.vert <- imfilter(img, gf1D, 'replicate')      # apply vertical blur
  img.vert.horiz <- imfilter(img.vert, t(gf1D), 'replicate')  # apply horizontal blur
  return(img.vert.horiz)
}

# CIE 1931 linear luminance Y
# ---------------------------
# RGB<->XYZ conversion in sRGB space using reference white D65:
#   |X|       |R|              |X|   |R|
#   |Y| = M * |G|        M(-1)*|Y| = |G|
#   |Z|       |B|              |Z|   |B|
# Full matrices:
#         |  0.4124564  0.3575761  0.1804375 |
# M =     |  0.2126729  0.7151522  0.0721750 |
#         |  0.0193339  0.1191920  0.9503041 |
#
#         |  3.2404542 -1.5371385 -0.4985314 |
# M(-1) = | -0.9692660  1.8760108  0.0415560 |
#         |  0.0556434 -0.2040259  1.0572252 |
#
# In particular Y coefficients for R, G, and B are 0.2126729, 0.7151522, and 0.0721750
# Many other choices: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
# I'm not sure whether gamma correction should be applied here.
luminance <- function(img) 0.212673*img[,,RED$idx] + 0.715152*img[,,GREEN$idx] + 0.072175*img[,,BLUE$idx]

# Adjusts the amount of color-channel-wise blur.
# The effects of blur.factor:
# Negatives add blur, positives reduce blur up to 1.0; larger than 1.0 overcompensate.
#   < -1: image extra blurred
#   = -1: the result is the blurred filter image
#   =  0: the result is identical to original image
#   = +1: blur removed (with a blur filter)
#   > +1: overcompensated blur removal
adjustBlur <- function(img, blur.factor, size=3) {
  blurred <- meanBlurred(img, size=3)
  return(pmin(pmax(img + (img-blurred)*blur.factor, 0), 1))
}

# My own sharpness map that is based on the basic blur filter
sharpAmount <- function(img, size=3) {
  blurred <- meanBlurred(img,size)
  return(abs(img - blurred))
}

# My own blur map that is based on the basic blur filter
blurAmount <- function(img, size=3) {
  blurred <- meanBlurred(img, size)
  return(1 - sharpAmount(img))
}


######################
## QUALITY MEASURES ##
######################

# The Blur Effect: Perception and Estimation with a New No-Reference Perceptual Blur Metric"
# Crete F., Dolmiere T., Ladret P., Nicolas M. - GRENOBLE - 2007
# In SPIE proceedings - SPIE Electronic Imaging Symposium Conf Human Vision and Electronic Imaging, tats-Unis d'Amrique (2007)
# https://hal.archives-ouvertes.fr/hal-00232709/document
# Error is possible (div by 0 in normalization)?
# Minimum value 0, best quality in terms of blur perception.
# Maximum value 1.0 is reported. Is it reachable? Checkerboard image: sum.gray=1.0*h*w, sum.cmp=h*w*0.5 yields 0.5 (approx)
# http://www.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric/content//blurMetric.m

blurAnnoyanceQuality <- function(img, f.len=9) {
  h <- dim(img)[1]
  w <- dim(img)[2]
  # Sharpness requires only one channel: luminance
  img.gray <- luminance(img)
  # Vertical and horizontal filters:
  f.vert <- matrix(1/f.len, nrow=f.len)
  f.horiz <- t(f.vert)
  blurred.vert <- imfilter(img.gray, f.vert, pad='replicate')[ , , 1]
  blurred.horiz <-imfilter(img.gray, f.horiz, pad='replicate')[ , , 1]
  # Differences in neighboring pixels (both within original and blurred):
  diff.gray.vert     <- abs(    img.gray[2:h, 1:w] -     img.gray[1:(h-1), 1:w])
  diff.blurred.vert  <- abs(blurred.vert[2:h, 1:w] - blurred.vert[1:(h-1), 1:w])
  diff.gray.horiz    <- abs(     img.gray[1:h, 2:w] -      img.gray[1:h, 1:(w-1)])
  diff.blurred.horiz <- abs(blurred.horiz[1:h, 2:w] - blurred.horiz[1:h, 1:(w-1)])
  # Difference between original and blurred pixels:
  cmp.vert <- pmax(0, diff.gray.vert - diff.blurred.vert)  # including all columns (not just 2:w)
  cmp.horiz <- pmax(0, diff.gray.horiz - diff.blurred.horiz)  # including all rows (not just 2:h)
  # Comparison between original (gray) and blurred: sum of differences
  sum.gray.vert <- sum(diff.gray.vert)
  sum.cmp.vert <- sum(cmp.vert)
  sum.gray.horiz <- sum(diff.gray.horiz)
  sum.cmp.horiz <- sum(cmp.horiz)
  # Normalization
  blur.vert <- (sum.gray.vert - sum.cmp.vert) / sum.gray.vert
  blur.horiz <- (sum.gray.horiz - sum.cmp.horiz) / sum.gray.horiz
  blur <- max(blur.vert, blur.horiz)  # from 0 to 1.0, smaller is better
  # Convert into a quality parameter (from 1 to 5, larger is better)
  # Constants come from the interpolation model parameters, as reported in the article.
  quality <- 3.79/(1 + exp(10.72*blur - 4.55)) + 1.13  # from 1 to 5, larger is better
  ## My own normalization
  normalized <- (quality-1)/4  # convert to 0--1 (larger is better)
  return(normalized)
}


# Helper function for the 'mdwe' method -- not a universal method.
# TODO: how to derive the magic constant? How many bins do we consider?
detectVerticalEdgeAreas <- function(img.gray) {
  f.vert <- sobelVertical()
  filtered <- imfilter(img.gray, f.vert, pad='replicate')[ , , 1]
  edge.strengths <- abs(filtered)
  bins <- 20  # granularity: just my own magic constant
  breaks <- seq(min(edge.strengths), max(edge.strengths), length=bins)
  h <- hist(c(edge.strengths), breaks=breaks, plot=FALSE)
  thres.idx <- otsuThreshold(h$counts)
  thres.val <- h$mids[thres.idx]
  edges <- edge.strengths >= thres.val
  return(edges)
}
# Returns a vector of midpoints (rounding up) for consecutive sub-sequences of TRUE.
# Input: logical vector
# For example: TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE  --> 1,5,9
# TODO: consider whether consecutive edges should allow one index without
#       edge indication, such as ...,0,0,0,1,1,1,1,0,1,1,1,0,0,...
sequenceMidpoints <- function(vec) {
  idxs <- which(vec)  # indices that are TRUE
  starts <- idxs[which(diff(c(-1, idxs)) != 1)]    # start indices of consecutive TRUEs, first TRUE included
  ends <- idxs[which(diff(c(idxs, length(vec)+2)) != 1)]  # end indices of consecutive TRUEs, last TRUE included
  mids <- as.integer(ceiling(starts + (ends - starts)/2))  # mid indices of consecutive TRUE sub-sequences
  return(mids)
}
localMinimaIndices <- function(vec) {
  is.minima <- function(xs, i) i == 1 || i == length(xs) || (xs[i-1] >= xs[i] && xs[i+1] >= xs[i])
  minima.idxs <- which(sapply(1:length(vec), function(i) is.minima(vec, i)))
  return(minima.idxs)
}
localMaximaIndices <- function(vec) {
  is.maxima <- function(xs, i) i == 1 || i == length(xs) || (xs[i-1] <= xs[i] && xs[i+1] <= xs[i])
  maxima.idxs <- which(sapply(1:length(vec), function(i) is.maxima(vec, i)))
  return(maxima.idxs)
}

# Pina Marziliano, Frederic Dufaux, Stefan Winkler and Touradj Ebrahimi:
# A No-Reference Perceptual Blur Metric, 2002
# TODO: derive a model to produce a quality index between 1 and 5, for example (based on data in the article)
# TODO: horizontal-edge version of the method
mdweVertical <- function(img) {
  img.gray <- luminance(img)  # use only the luminosity image
  edges.areas <- detectVerticalEdgeAreas(img.gray)  # several consecutive pixels may be marked as belonging to an edge
  total.width <- 0
  num.edges <- 0
  for (r in 1:nrow(img)) {  # for each row
    # Find edge indices:
    edge.mids <- sequenceMidpoints(edges.areas[r, ])
    # Find all local maxima and minima for the luminosity vector r:
    horiz.vec <- img.gray[r, ]  # luminosities for row pixels
    minima.idxs <- localMinimaIndices(horiz.vec)
    maxima.idxs <- localMaximaIndices(horiz.vec)
    # For each edge mid-point, find a local minimum and maximum so that mid-point is enclosed:
    min.idxs <- findInterval(edge.mids, minima.idxs, all.inside=TRUE)
    max.idxs <- findInterval(edge.mids, maxima.idxs, all.inside=TRUE)
    prev.mins <- minima.idxs[min.idxs]  # preceding local minimum
    prev.maxs <- maxima.idxs[max.idxs]
    # From previous minimum and maximum, deduce the next minimum or maximum.
    # It is either the next index or the same index (when prev idx is last or edge-mid is at a local optimum).
    not.same.mins <- (edge.mids != prev.mins)*(length(minima.idxs)+1)  # impossibly large when edge.mid != prev.min
    not.same.maxs <- (edge.mids != prev.maxs)*(length(maxima.idxs)+1)  # impossibly large when edge.mid != prev.max
    next.mins <- minima.idxs[pmin(min.idxs + 1, length(minima.idxs), prev.mins + not.same.mins)]  # next, last, or same
    next.maxs <- maxima.idxs[pmin(max.idxs + 1, length(maxima.idxs), prev.maxs + not.same.maxs)]  # next, last, or same
    # Compute widths for all edge intervals on this row:
    widths <- pmin(next.mins - prev.maxs,  next.maxs - prev.mins)  # choose the minimum width: maxmin or minmax
    total.width <- total.width + sum(widths)
    num.edges <- num.edges + length(edge.mids)
  }
  score.raw <- total.width/num.edges  # Inf when num.edges==0
  ## My own normalization
  score.gaussian <- limit(score.raw^1.04 - 4.5, 1, 10)  # From no blur (1) to distractive blur (10)
  score.jpeg2k <- limit(17.5*max(0,score.raw-3)^0.3 - 19.5, 1, 10)
  norm <- function(x) -((score.gaussian - 1) / (10-1) - 1)  # 1--10 (smaller is better) --> 0--1 (larger is better)
  return(list(score=score.raw, gaussian=norm(score.gaussian), jpeg2k=norm(score.jpeg2k)))
}


# A New No-reference Method for Color Image Quality Assessment
# Sonia Ouni, Ezzeddine Zagrouba, Majed Chambah, 2012
# International Journal of Computer Applications (0975 ? 8887)
# Volume 40? No.17, February 2012

# Note that there are errors in the RGB->HSV conversion formula. See Wikipedia instead.
# The article contains numerous errors, the examples don't match with definitions,
# and the exact methods are sometimes vague. I've done my best to develop a method
# that is both practical and close to the spirit of the authors' ideas.

# COMMENTS:

## Kappa is the same as the following (measure of "compactness"):
#C <- sum(cos(hues[mask]))
#S <- sum(sin(hues[mask]))
#R <- sqrt(C^2 + S^2)
#n <- length(hues[mask])
#kappa <- A1inv(R/n)

## The exact computation of DS leads to a huge computational task.
#k <- length(row.idx)
#ds <- 0
#for (i in 1:k) {
#  local.dist <- sum(sqrt((row.idx - row.idx[i])^2 + (col.idx - col.idx[i])^2)) / k
#  ds <- ds + local.dist
#}
#spatial.dispersion <- ds / k
## OR
#spatial.dispersion <- sum(sapply(1:length(row.idx), function(i) sqrt((row.idx - row.idx[i])^2 + (col.idx - col.idx[i])^2))) / npdc^2

## This is an average distance between pixels that have dominant color.
## I still think the value should be normalized somehow, for example dividing
## by an average distance between two pixels in the whole image, not just
## in dominant color. Otherwise image size has an effect on the absolute
## magnitude of the number. In the original article the values
## like 0.0020936 don't quite match with the definitions.

dispersionDominantColor <- function(img.hsv) {
  hues <- degreeToRadian(extractHSVChannel(img.hsv, HUE))  # convert from degrees to radians
  nr <- nrow(hues)
  nc <- ncol(hues)
  mask <- !is.na(hues)  # existing hues, excluding white and black
  mu <- atan2(sum(sin(hues[mask])), sum(cos(hues[mask])))%%(2*pi)  # angles between hues
  # Estimate kappa parameter on a von Mises distribution in a circular domain.
  kappa <- A1inv(mean(cos(hues[mask] - mu)))

  # NPDC and pi measures
  ######################
  # distances between hues (in radians)
  dist.raw <- abs(hues[mask] - mu)
  dist <- pmin(dist.raw, 2*pi - dist.raw)  # distance between hues, mod 2*pi
  npdc.mask <- dist <= 1/kappa  # article has a different interpretation of 1/kappa ??
  npdc <- sum(npdc.mask)  # pixels in a dominant color (within range kappa)
  pi.measure <- npdc / length(hues[mask])
  px.row <- matrix(1:length(hues) %% nr, nrow=nr)
  px.col <- matrix(1:length(hues) %/% nr, nrow=nr)
  
  # Spatial dispersion, by approximation
  ######################################
  row.idx <- c(px.row[mask][npdc.mask])
  col.idx <- c(px.col[mask][npdc.mask])
  npdc.nr <- length(row.idx)
  npdc.nc <- length(col.idx)

  # Approximation for all Euclidean distances, for custom normalization
  num.samples <- 10000  # just enough to get a hint about the avg distance
  smpl <- function(values) sample(values, num.samples, replace=TRUE)
  approx.dist <- mean(sqrt((smpl(1:nr) - smpl(1:nr))^2 + (smpl(1:nc) - smpl(1:nc))^2))

  # NPDC Euclidean distance approximation by random sample
  px.1 <- smpl(1:npdc.nr)
  px.2 <- smpl(1:npdc.nr)
  npdc.approx.dist <- mean(sqrt((row.idx[px.1] - row.idx[px.2])^2 + (col.idx[px.1] - col.idx[px.2])^2))

  spatial.dispersion <- npdc.approx.dist
  custom.ds <- min(1, npdc.approx.dist / approx.dist)  # cap to 1 (exceeds by random variation)

  return(list(mu=mu, kappa=kappa, pi=pi.measure, ds=spatial.dispersion, custom.ds=custom.ds))
}



##########
## MAIN ##
##########

main <- function() {
  do.view <- FALSE
  print('=== START ===')
  filename <- '../examples/small_grid.png'
  #filename <- '../examples/blue_shift.png'
  #filename <- '../examples/no_shift.png'
  #filename <- '../examples/niemi.png'
  #filename <- '../examples/sharp_or_blur.png'  # Blur annoyance quality (1--5): 1.17416513963911"
  #filename <- '../examples/K5_10994.JPG'
  #filename <- '../examples/green_grass_blue_sky.png'
  #filename <- '../examples/dark_city.png'
  #filename <- '../examples/violetred.png'
  #filename <- '../examples/bluehue.png'
  #filename <- '../examples/temple_set/temple-a-original.png'
  #filename <- '../examples/temple_set/temple-b-blue.png'
  #filename <- '../examples/temple_set/temple-c-cyan.png'
  #filename <- '../examples/temple_set/temple-d-yellow.png'
  #filename <- '../examples/temple_set/temple-e-magenta.png'
  #filename <- '../examples/temple_set/temple-f-red.png'
  #filename <- '../examples/temple_set/temple-g-green.png'
  #filename <- '../examples/temple_set/temple-h-noise.png'
  #filename <- '../examples/temple_set/temple-i-contrast.png'
  #filename <- '../examples/temple_set/temple-j-colornoise.png'
  #filename <- '../examples/temple_set/temple-k-gaussianblur.png'
  
  img <- readImage(filename)
  print(paste('Processing image:', filename))
  for (channel in RGB) {
    if (do.view) view(extractRGBChannel(img,channel), title=paste(channel$name, 'color channel'))
  }
  if (do.view) view(luminance(img), title='Luminance')
  blurred <- meanBlurred(img)
  for (channel in RGB) {
    if (do.view) view(extractRGBChannel(blurred,channel), title=paste('Blurred in', channel$name, 'color channel'))
  }
  if (do.view) view(rescaleImage(blurred), title='Blurred')
  blurMap <- blurAmount(img, 5)
  if (do.view) view(blurMap, title='Blurriness')
  sharpMap <- sharpAmount(img, 5)
  if (do.view) view(sharpMap, title='Sharpness')
  if (do.view) view(adjustBlur(img, 1.0, 5), title='Blur plus 1.0')
  if (do.view) view(adjustBlur(img, -1.0, 5), title='Blur minus 1.0')
  gaussian.blurred <- gaussianBlurred(img, 1.0)
  if (do.view) view(gaussian.blurred, title='Gaussian blurred')
  if (do.view) view(img, title='Full RGB image')
  blur <- blurAnnoyanceQuality(img, f.len=9)
  print(paste('Blur annoyance quality (0--1):', blur))  # more is better
  mdwe.score <- mdweVertical(img)
  print(paste('MDWE horizontal blur width:', mdwe.score[['score']]))  # smaller is better
  print(paste('MDWE Gaussian quality (0--1):', mdwe.score[['gaussian']]))  # greater is better
  print(paste('MDWE JPEG2000 quality (0--1):', mdwe.score[['jpeg2k']]))  # greater is better
  ddc <- dispersionDominantColor(toHSV(img))

  # Dominant direction, spread, and portion of the dominant color.
  # Dominant hue as an angle
  print(paste('Color dispersion(mu):', ddc$mu, 'and in degrees', radianToDegree(ddc$mu)))
  # Compactness of hues, larger kappa means narrower hues
  print(paste('Color dispersion(kappa):', ddc$kappa))
  # Proportion of pixels whose hue is close to dominant
  print(paste('Color dispersion(pi):', ddc$pi))
  # Spatial dispersion: mean distance of pixels that are dominant
  print(paste('Color dispersion(ds):', ddc$ds, 'px'))
  # Custom normalized spatial dispersion: dominantdistances / alldistances
  print(paste('Color dispersion(custom.ds):', ddc$custom.ds))
  if (do.view) viewDDC(filename='polarcolor.png', toHSV(img), ddc$mu, ddc$kappa)
  print('===  END  ===')
}

main()



# Measuring 2D Gaussian blur:
#img <- readImage('../examples/small_grid.png')
#s <- Sys.time()
#replicate(1000, gaussianBlurred(img, 1.0))
#print(Sys.time() - s)

#print('Adjusted blur:')
#print(blurAnnoyanceQuality(readImage('../examples/small_grid.png'), f.len=3))  # conv2 & conv2dFFT: 2.002251 quality

#s <- Sys.time()
#replicate(50, adjustBlur(readImage('../examples/sharp_or_blur.png'), 1.0, 5))
#print(Sys.time() - s)
# With a very small filter matrix:
#   conv2:     1.245071 secs
#   conv2dFFT: 1.492085 secs
# I suppose it takes a larger filter matrix to pay off using FFT.

#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(c(1/3), nrow=3, ncol=1)))
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]    7    8    9   10   11   12
#[2,]   13   14   15   16   17   18
#[3,]   19   20   21   22   23   24
#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(c(1/3), nrow=1, ncol=3)))
#     [,1] [,2] [,3] [,4]
#[1,]    2    3    4    5
#[2,]    8    9   10   11
#[3,]   14   15   16   17
#[4,]   20   21   22   23
#[5,]   26   27   28   29

#n<-9; vals <- (1:n)/(n*(n+1)/2)
#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(vals, nrow=3, ncol=3)))
#     [,1] [,2] [,3] [,4]
#[1,]  6.8  7.8  8.8  9.8
#[2,] 12.8 13.8 14.8 15.8
#[3,] 18.8 19.8 20.8 21.8

# Marziliano data (estimated from plots):
#   Gaussian comparison (x: objective, y: subjective):
#   gx <- c(4.5, 4.5, 4.6, 4.7, 4.7, 4.8, 4.9, 5.1, 5.3, 5.4, 5.9, 6.0, 6.1, 6.4, 7.2, 7.3, 7.6, 8.0, 8.4, 8.6, 9.2, 9.4, 9.6, 9.8, 10.9, 11.2, 11.3, 11.4, 13.0, 13.4)
#   gy <- c(0.1, 0.2, 0.6, 0.8, 0.6, 0.3, 0.4, 0.3, 0.7, 1.2, 2.9, 2.8, 2.3, 2.5, 3.8, 2.3, 5.1, 4.9, 4.7, 6.0, 7.3, 5.5, 6.5, 7.3, 8.7, 7.0, 8.1, 7.6, 8.9, 8.0)
#   JPEG comparison (x: objective, y: subjective):
#   jx <- c(4.5, 4.5, 4.6, 4.8, 5.1, 5.2, 5.3, 5.3, 5.7, 5.7, 5.8, 5.9, 5.9, 6.2, 6.2, 6.3, 6.6, 6.7, 6.8, 6.8, 7.0, 7.0, 7.1, 7.1, 7.2, 7.7, 7.7, 7.8, 8.4, 9.2)
#   jy <- c(0.3, 0.6, 0.1, 0.8, 2.5, 0.9, 4.5, 3.5, 1.6, 4.8, 3.6, 7.1, 5.8, 6.0, 3.3, 8.5, 5.1, 9.2, 7.1, 6.6, 7.2, 6.1, 7.4, 8.7, 8.8, 8.1, 8.1, 8.0, 8.0, 8.9)
#  My own functions: objective x -> subjective y (limited to 1--10 disturbance, higher is worse)
#  Gaussian: y = x^1.04 - 4.5
#  JPEG: y = 17.5*(x-3)^(0.3)-19.5

# Von Mises test data:
# http://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/MiscellaniousTheses/Bentley-2006.pdf
#   degrees <- c(0,0,0,15,45,68,100,110,113,135,135,140,140,155,165,165,169,180,180,180,180,180,180,180,189,206,209,210,214,215,225,226,230,235,245,250,255,255,260,260,260,260,270,270)
#   rads <- degrees / 360 * (2*pi)
#   mu <- atan2(sum(sin(rads)), sum(cos(rads))) %% (2*pi)
#   kappa <- est.kappa(rads)
# Compare:
#   hist(degrees)
#   hist(as.numeric(rvonmises(10000, mu, kappa)) / (2*pi) * 360)
# Expected: maximum likelihood parameter estimate:
#   mu; stderr(mu); kappa; stderr(kappa)
#   199.4 degrees; 12.2 degrees; 1.07; 0.26
#
# Another test data:
#   degrees <- c(1.9, 12.4, 28.1, 28.9, 41.5, 46.0, 55.5, 56.6, 72.6, 75.5,
#                86.1, 109.6, 111.0, 115.3, 123.3, 127.6, 139.6, 140.8, 142.0, 147.5,
#                147.7, 149.8, 150.3, 154.1, 160.0, 161.9, 162.1, 162.4, 162.7, 163.1,
#                163.7, 168.1, 170.2, 170.4, 171.9, 172.2, 172.5, 175.4, 175.6, 175.7,
#                176.5, 177.1, 177.7, 179.0, 179.4, 179.7, 180.6, 180.7, 180.7, 181.1,
#                181.7, 182.0, 183.8, 184.1, 185.3, 188.5, 188.8, 189.0, 189.8, 192.6,
#                193.9, 194.9, 195.5, 195.7, 195.9, 196.0, 196.2, 196.4, 196.6, 198.0,
#                199.2, 199.8, 202.4, 202.9, 204.8, 206.7, 207.6, 210.5, 210.9, 212.4,
#                212.5, 213.1, 214.8, 218.0, 219.6, 220.6, 220.7, 224.5, 228.2, 232.4,
#                254.5, 255.0, 266.0, 277.4, 282.9, 289.4, 295.4, 301.1, 326.4, 354.9)
#   rads <- degrees / 360 * (2*pi)
#   mu <- atan2(sum(sin(rads)), sum(cos(rads))) %% (2*pi)
#   kappa <- est.kappa(rads)
# Expected: maximum likelihood parameter estimate:
#   mu; stderr(mu); kappa; stderr(kappa)
#   183.3 degrees; 5.9 degrees; 1.55; 0.21
