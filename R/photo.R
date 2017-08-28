##
## Photo and image representations
##
## Esa Junttila, 2016-03-23

# Color channel identifiers:
RED <- list(idx=1, name='RED')
GREEN <- list(idx=2, name='GREEN')
BLUE <- list(idx=3, name='BLUE')
RGB <- list(RED, GREEN, BLUE)

HUE <- list(idx=1, name='HUE')
SATURATION <- list(idx=2, name='SATURATION')
VALUE <- list(idx=3, name='VALUE')
HSV <- list(HUE, SATURATION, VALUE)

CIEXYZ_X <- list(idx=1, name='CIE_X')  # "red--green"
CIEXYZ_Y <- list(idx=2, name='CIE_Y')  # "luminosity"
CIEXYZ_Z <- list(idx=3, name='CIE_Z')  # "blue--yellow"
CIEXYZ_XYZ <- list(CIEXYZ_X, CIEXYZ_Y, CIEXYZ_Z)

CIExyY_x <- list(idx=1, name='CIE_x')  # "color coordinate"
CIExyY_y <- list(idx=2, name='CIE_y')  # "luminosity"
CIExyY_Z <- list(idx=3, name='CIE_Z')  # "color coordinate"
CIExyY_xyY <- list(CIExyY_x, CIExyY_y, CIExyY_Z)

CIELUV_L <- list(idx=1, name='L')
CIELUV_U <- list(idx=2, name='U')
CIELUV_V <- list(idx=3, name='V')
CIELUV <- list(CIELUV_L, CIELUV_U, CIELUV_V)

# Image construction from separate channels
createImageRGB <- function(reds, greens, blues) {
  if (!all(dim(reds) == dim(blues)) | !all(dim(blues) == dim(greens))) {
    stop('Mismatch in matrix sizes from R, G, and B channels.')
  }
  dim.names <- list(NULL, NULL, c('Red', 'Green', 'Blue'))
  img <- array(dim=c(nrow(reds), ncol(reds), length(RGB)), dimnames=dim.names)
  img[ , , RED$idx] <- reds
  img[ , , GREEN$idx] <- greens
  img[ , , BLUE$idx] <- blues
  return(img)
}

createImageHSV <- function(hues, saturations, values) {
  if (!all(dim(hues) == dim(saturations)) | !all(dim(saturations) == dim(values))) {
    stop('Mismatch in matrix sizes from H, S, and V channels.')
  }
  dim.names <- list(NULL, NULL, c('Hue', 'Saturation', 'Value'))
  img <- array(dim=c(nrow(hues), ncol(hues), length(HSV)), dimnames=dim.names)
  img[ , , HUE$idx] <- hues
  img[ , , SATURATION$idx] <- saturations
  img[ , , VALUE$idx] <- values
  return(img)
}


# Return a matrix that represents channel 'col.channel' of image 'img'.
extractRGBChannel <- function(img.rgb, col.channel) img.rgb[ , , col.channel$idx]

extractHSVChannel <- function(img.hsv, hsv.channel) img.hsv[ , , hsv.channel$idx]

extractXYZChannel <- function(img.xyz, xyz.channel) img.xyz[ , , xyz.channel$idx]

extractxyYChannel <- function(img.xyY, xyY.channel) img.xyY[ , , xyY.channel$idx]

extractLUVChannel <- function(img.luv, luv.channel) img.luv[ , , luv.channel$idx]
