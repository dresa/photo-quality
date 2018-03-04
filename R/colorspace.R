##
## Color space transformation between color spaces, such as
## RGB, HSV, XYZ, xyY, and LUV.
##
## Esa Junttila, 2016-03-23

source('R/photo.R')

#' Convert RGB image to an HSV image.
#'
#' Convert a 2D image from RGB (Red, Green, Blue) format into HSV
#' (Hue, Saturation, Value) format. In both formats, we assume that
#' the images are represented as three-dimensional arrays: height, width,
#' and three channels.
#'
#' @param img.rgb 2D image in RGB format: three-dimensional array with
#'   height, width, and channel. RGB values are within \code{[0;1]},
#'   unless maximum value is specified in \code{max.value}.
#' @param radians if TRUE, returns Hue as radians \code{[0;2*pi]};
#'   otherwise as degrees \code{[0;360]} (default)
#' @param max.value maximum valid RGB values, such as 1 or 255 (default 1)
#' @return Image in HSV formatted three-dimensional array, with height,
#'   width, and HSV channels. Saturation and Value are within \code{[0;1]}
#'   and Hues either in \code{[0;360]} degrees or \code{[0;2*pi]} radians.
#'   Whenever RGB values indicate black or white, HSV Hue will be \code{NA}.
#' @examples
#'   set.seed(1)
#'   d <- c(4,12,3)
#'   rgb.img <- array(sample(0:255, prod(d), replace=TRUE), dim=d) / 255
#'   hsv.img <- toHSV(rgb.img)
#'
#'   rgb.px <- array(c(200, 100, 150), dim=c(1,1,3)) / 255
#'   err.hsv <- max(abs(as.vector(toHSV(rgb.px)) - as.vector(rgb2hsv(200, 100, 150))*c(360,1,1)))
#'   err.rgb <- max(abs(as.vector(toRGB(toHSV(rgb.px))) - as.vector(rgb.px)))
#'   stopifnot(err.rgb < 1e-14)
#'
#'   stopifnot(all(abs(toHSV(rgb.px) - toHSV(rgb.px*255, max.value=255) < 1e-14)))
#'   stopifnot(abs(toHSV(rgb.px) - toHSV(rgb.px, radians=TRUE)*c(360/(2*pi), 1, 1)) < 1e-14)
#' @export
toHSV <- function(img.rgb, radians=FALSE, max.value=1) {
  stopifnot(is.logical(radians))
  stopifnot(max.value > 0)
  # Extract R,G, and B color channels
  red <- extractRGBChannel(img.rgb, RED)
  green <- extractRGBChannel(img.rgb, GREEN)
  blue <- extractRGBChannel(img.rgb, BLUE)
  # Temporary variables
  dims <- if (is.matrix(red)) dim(red) else c(length(red), 1)
  max.channels <- pmax(red, green, blue)
  diffs <- max.channels - pmin(red, green, blue)
  # Compute H(ue), S(saturation), and V(alue):
  # S(aturation):
  saturations <- diffs / max.channels
  saturations[max.channels == 0] <- 0  # replace Inf by zero (originates from division by zero)
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
  hues[max.channels==0] <- NA  # undefined hue is not needed; NA as default
  hues[diffs==0] <- NA  # undefined hue is not needed; NA as default
  # V(alue):
  values <- max.channels / max.value
  # Return an HSV image as an (M x N x 3) array
  dim.names <- list(NULL, NULL, c('Hue', 'Saturation', 'Value'))
  return(array(c(hues, saturations, values), dim=c(dims[1], dims[2], 3), dimnames=dim.names))
}


# Convert HSV image to a RGB image.
# Assumptions:
#   RGB values are between 0 and 1.
#   Hues are degrees between 0 and 360 (or marked as radians);
#   Saturations and Values are between 0 and 1.
toRGB <- function(img.hsv, radians=FALSE, na.hue=0) {
  # Extract H, S, and V color channels
  hue <- extractHSVChannel(img.hsv, HUE)
  saturation <- extractHSVChannel(img.hsv, SATURATION)
  value <- extractHSVChannel(img.hsv, VALUE)
  # Preprocess hues: replace undefined hues with default red (zero in radians and degrees).
  # Actual hue should not matter since value or saturation should be zero anyway.
  hue[is.na(hue)] <- na.hue
  # Convenience variables
  dims <- if (is.matrix(hue)) dim(hue) else c(length(hue), 1)
  chroma <- value * saturation
  angle <- if (radians) pi/3 else 60  # 60 degrees
  x <- chroma * (1 - abs(((hue / angle) %% 2)- 1))
  m <- value - chroma

  # Conversion pre-processing
  red.add <- array(dim=dims)
  green.add <- array(dim=dims)
  blue.add <- array(dim=dims)
  hue.group <- hue %/% angle
  ## proceed by groups
  for (group.id in 0:5) {
    group.px <- hue.group == group.id
    C <- chroma[group.px]
    X <- x[group.px]
    red.add[group.px] <-   switch(group.id + 1, C, X, 0, 0, X, C)
    green.add[group.px] <- switch(group.id + 1, X, C, C, X, 0, 0)
    blue.add[group.px] <-  switch(group.id + 1, 0, 0, X, C, C, X)
  }
  # Wrap it up
  return(createImageRGB(m + red.add, m + green.add, m + blue.add))
}

# Convert an RGB image into the "CIE XYZ" color space
# Reference: http://ninedegreesbelow.com/photography/xyz-rgb.html
# http://dougkerr.net/Pumpkin/articles/CIE_XYZ.pdf
# https://en.wikipedia.org/wiki/SRGB
# Checked with http://colormine.org/color-converter ?
# In ICC, the RGB color white (1,1,1) has the XYZ coordinates (0.9642, 1.0000, 0.8249).
# We are using the D65 XYZ coordinates, where RGB(1,1,1) translates to XYZ(0.9505 1.0000 1.0891).
# For comparison, values for D50:
#                    Red = RGB (1,0,0)      Green = RGB (0,1,0)      Blue = RGB (0,0,1)
#                    X       Y        Z      X        Y      Z       X        Y       Z
#   sRGB (D65)    0.4358  0.2224   0.0139  0.3853  0.7170  0.0971  0.1430  0.0606  0.7139
#   CIE-RGB (E)   0.4685  0.1699  -0.0007  0.3274  0.8242  0.0131  0.1683  0.0059  0.8125
toXYZ <- function(img.rgb, make.linear=TRUE) {
  adjustLinear <- function(x) {
    if (!make.linear) return(x)
    K0 <- 0.04045  # threshold; notation coming from Wikipedia definition
    phi <- 12.92
    a <- 0.055
    gamma <- 2.4
    lin <- ((x + a) / (1 + a))^gamma
    lin[x<=K0] <- x[x<=K0]/phi  # correction for small RGB values
    return(lin)
  }
  M <- matrix(c(0.4124564, 0.3575761, 0.1804375,
                0.2126729, 0.7151522, 0.0721750,
                0.0193339, 0.1191920, 0.9503041), 3, byrow=TRUE)
  dim.names <- list(NULL, NULL, c('X', 'Y', 'Z'))
  img.xyz <- array(matrix(adjustLinear(img.rgb), ncol=3) %*% t(M), dim=dim(img.rgb), dimnames=dim.names)
  return(img.xyz)
  # Alternative, almost as fast method (should add adjustLinear here, too):
  #  red <- extractRGBChannel(img.rgb, RED)
  #  green <- extractRGBChannel(img.rgb, GREEN)
  #  blue <- extractRGBChannel(img.rgb, BLUE)
  #  X <- 0.412453*red + 0.35758 *green + 0.180423*blue
  #  Y <- 0.212671*red + 0.71516 *green + 0.072169*blue
  #  Z <- 0.019334*red + 0.119193*green + 0.950227*blue
  #  img.xyz <- array(c(X, Y, Z), dim=dim(img.rgb), dimnames=dim.names)
  #
  # Inversion to be used in XYZ-->RGB conversions.
  # Minv <- matrix(c(3.2404542, -1.5371385, -0.4985314,
  #                  -0.9692660, 1.8760108, 0.0415560,
  #                  0.0556434, -0.2040259, 1.0572252), 3, byrow=TRUE)
}


# RGB --> "CIE xyY" color space conversion
# http://ninedegreesbelow.com/photography/xyz-rgb.html
# http://colormine.org/color-converter
toxyY <- function(img.rgb) {
  img.xyz <- toXYZ(img.rgb)
  X <- extractXYZChannel(img.xyz, CIEXYZ_X)
  Y <- extractXYZChannel(img.xyz, CIEXYZ_Y)
  Z <- extractXYZChannel(img.xyz, CIEXYZ_Z)
  s <- X+Y+Z
  x <- X / s
  y <- Y / s
  dim.names <- list(NULL, NULL, c('x', 'y', 'Y'))
  return(array(c(x, y, Y), dim=dim(img.rgb), dimnames=dim.names))
}

# Convert RGB image to an LUV image. RGB values are assumed to be between 0 and 1.
# http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Luv.html
# http://www.brucelindbloom.com/index.html?LContinuity.html
# http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html
# http://colormine.org/color-converter
# Note that for blacks the U and V components are undefined.
toLUV <- function(img.rgb) {
  # Define D65 white point: CIE chromaticity coordinates and CIE luminance
  xn <- 0.312713
  yn <- 0.329016
  Yn <- 1.0
  # Junction point correction constants:
  eps <- 216/24389
  kappa <- 24389/27
  # Other constants:
  un <- 4*xn / (-2*xn + 12*yn + 3)
  vn <- 9*yn / (-2*xn + 12*yn + 3)
  # Compute LUV
  img.xyz <- toXYZ(img.rgb, make.linear=TRUE)  # First convert RGB --> XYZ
  X <- extractXYZChannel(img.xyz, CIEXYZ_X)
  Y <- extractXYZChannel(img.xyz, CIEXYZ_Y)
  Z <- extractXYZChannel(img.xyz, CIEXYZ_Z)
  u <- 4*X / (X + 15*Y + 3*Z)
  v <- 9*Y / (X + 15*Y + 3*Z)
  u[is.nan(u)] <- NA  # black color has undefined u and v; use NA instead of NaN
  v[is.nan(v)] <- NA
  r <- Y/Yn
  L <- 116 * r^(1/3) - 16
  mask <- r <= eps
  L[mask] <- kappa*r[mask]  # correction
  U <- 13*L*(u-un)
  V <- 13*L*(v-vn)
  img.luv <- array(c(L,U,V), dim=dim(img.rgb), dimnames=list(NULL, NULL, c('L','U','V')))
  return(img.luv)
  # Computed L component values are in the range [0 to 100].
  # Computed U component values are in the range [-124 to 220], unless black.
  # Computed V component values are in the range [-140 to 116], unless black.
}

test <- function() {
  set.seed(1)
  d <- c(4,12,3)
  rgb.vec <- sample(0:255, prod(d), replace=TRUE)
  rgb <- array(rgb.vec, dim=d) / 255

  # reference data: use R's own RGB->HSV->RGB conversion
  hsv.ref <- rgb2hsv(255*matrix(c(rgb[,,1], rgb[,,2], rgb[,,3]), nrow=3, byrow=TRUE))
  rgb.ref <- col2rgb(hsv(hsv.ref[1, ], hsv.ref[2, ], hsv.ref[3, ]))
  roundtrip.ok <- all(rgb == array(c(rgb.ref[1, ], rgb.ref[2, ], rgb.ref[3, ]), dim=d) / 255)
  if (!roundtrip.ok) stop('internal test failed: RGB->HSV->RGB conversion')

  # custom implementations
  hsv.val <- toHSV(rgb)
  err.hsv <- max(abs(hsv.ref - matrix(c(hsv.val[,,1]/360, hsv.val[,,2], hsv.val[,,3]), nrow=3, byrow=TRUE)))
  print(paste('Matrix RGB->HSV error', err.hsv))
  rgb.roundtrip <- toRGB(hsv.val)
  err.rgb <- max(abs(rgb.roundtrip - rgb))
  print(paste('Matrix HSV->RGB error', err.rgb))

  # just one
  rgb <- array(c(200, 100, 150), dim=c(1,1,3)) / 255
  err.hsv <- max(abs(as.vector(toHSV(rgb)) - as.vector(rgb2hsv(200, 100, 150))*c(360,1,1)))
  print(paste('Single HSV error', err.hsv))
  err.rgb <- max(abs(as.vector(toRGB(toHSV(rgb))) - as.vector(rgb)))
  print(paste('Single RGB error', err.rgb))
}

speedTest <- function() {
  set.seed(1)
  d <- c(1000,1500,3)
  rgb.vec <- sample(0:255, prod(d), replace=TRUE)
  rgb <- array(rgb.vec, dim=d) / 255

  refFunc <- function(rgb) {
    hsv.ref <- rgb2hsv(matrix(c(rgb[,,1], rgb[,,2], rgb[,,3]), nrow=3, byrow=TRUE))
    rgb.ref <- col2rgb(hsv(hsv.ref[1, ], hsv.ref[2, ], hsv.ref[3, ]))
    return(array(c(rgb.ref[1, ], rgb.ref[2, ], rgb.ref[3, ]), dim=d))
  }

  st <- Sys.time(); cmp <- max(abs(toRGB(toHSV(rgb)) - rgb)); et <- Sys.time();
  print(paste('Speed test: maxerror', cmp, ', time', et - st))

  st.ref <- Sys.time(); cmp.ref <- max(abs(refFunc(255*rgb) - 255*rgb)); et.ref <- Sys.time();
  print(paste('Reference speed test: maxerror', cmp.ref, ', time', et.ref - st.ref))
}

# Test color space conversion on a single point
# Reference values are from http://colormine.org/color-converter and R library
example.pixel <- createImageRGB(matrix(165/255), matrix(107/255), matrix(11/255))
stopifnot(abs(255*example.pixel - c(165, 107, 11)) < 1e-5)
stopifnot(abs(toHSV(example.pixel) - c(37.4025965757828, 0.933333333333333, 0.647058823529412)) < 1e-5)
stopifnot(abs(toXYZ(example.pixel)*100 - c(20.8351499726374, 18.53888482291, 2.7968391383824)) < 1e-2)  # approximately
stopifnot(abs(toxyY(example.pixel) - c(0.494064932239259, 0.439613484225046, 18.53888482291/100)) < 1e-4)
stopifnot(abs(toLUV(example.pixel) - c(50.1432998776774, 47.8174152455571, 48.6306348377325)) < 1e-2)
stopifnot(abs(c(toLUV(example.pixel)) - convertColor(matrix(c(example.pixel), 1, byrow=TRUE), from="sRGB", to="Luv")) < 0.4) # poor accuracy for an unknown reason

