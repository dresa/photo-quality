##
## Derived images, processed from an original image.
##
## Esa Junttila, 2016-03-23


source('R/convolutionfilter.R')  # filters (kernels) for 2D convolution
source('R/imageconvolution.R')  # 2D convolution functions
source('R/colorspace.R')  # convert RGB to HSV

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
  return(1 - sharpAmount(img))
}


# Transform RGB colors so that image contains pure bright hues.
toBrightRGB <- function(img.rgb, na.color='gray') {
  na.color.hsv <- rgb2hsv(col2rgb(na.color))  # Pixels with NA hues are filled with this color
  img.hsv <- toHSV(img.rgb)
  hues <- extractHSVChannel(img.hsv, HUE)
  saturations <- matrix(1, nrow=nrow(hues), ncol=ncol(hues))
  saturations[is.na(hues)] <- na.color.hsv[SATURATION$idx]
  values <- matrix(1, nrow=nrow(hues), ncol=ncol(hues))
  values[is.na(hues)] <- na.color.hsv[VALUE$idx]
  img.rgb.bright <- toRGB(createImageHSV(hues, saturations, values))
  return(img.rgb.bright)
}

