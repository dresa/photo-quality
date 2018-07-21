##
## Photo and image representations in color spaces, such as RGB, HSV, ...
##
## Currently we use (m x n x 3) arrays with three "color" channels as matrices.
##
## Esa Junttila, 2016-03-23 (originally)


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
CIExyY_Y <- list(idx=3, name='CIE_Y')  # "color coordinate"
CIExyY_xyY <- list(CIExyY_x, CIExyY_y, CIExyY_Y)

CIELUV_L <- list(idx=1, name='L')
CIELUV_U <- list(idx=2, name='U')
CIELUV_V <- list(idx=3, name='V')
CIELUV <- list(CIELUV_L, CIELUV_U, CIELUV_V)


#' Image construction from separate RGB channels.
#' 
#' Create an internal representation of a color image,
#' given as three RGB (Red, Green, Blue) channel matrices.
#' All channel matrices must have equal dimensions.
#' Accepts both numeric \code{[0;1]} and integer \code{{0,1,...,255}} ranges.
#' 
#' @param reds matrix of the red channel values.
#' @param greens matrix of the green channel values.
#' @param blues matrix of the blue channel values.
#' @return RGB image construct; use \code{\link{extractRGBChannel}} to extract
#'         separate \code{RED}, \code{GREEN}, and \code{BLUE} channels.
#' @examples
#'   r <- matrix(c(0.1, 0.2, 0.3,  0.4, 0.5, 0.6), byrow=TRUE, nrow=2)
#'   g <- matrix(c(0.0, 0.2, 0.4,  0.6, 0.8, 1.0), byrow=TRUE, nrow=2)
#'   b <- matrix(c(0.40, 0.44, 0.48,  0.52, 0.56, 0.60), byrow=TRUE, nrow=2)
#'   img.rgb <- createImageRGB(r, g, b)
#' @seealso \code{\link{extractRGBChannel}}
#' @export
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


#' Image construction from separate HSV channels.
#' 
#' Create an internal representation of a color image,
#' given as three HSV (Hue, Saturation, Value) channel matrices.
#' All channel matrices must have equal dimensions.
#' 
#' @param hues matrix of hues as numeric \code{[0;360]}
#'             degrees or \code{[0;2*pi]} radians.
#' @param saturations matrix of saturations a within numeric \code{[0;1]}.
#' @param values matrix of HSV Values within numeric \code{[0;1]}.
#' @return HSV image construct; use \code{\link{extractHSVChannel}} to extract
#'         separate \code{HUE}, \code{SATURATION}, and \code{VALUE} channels.
#' @examples
#'   h <- matrix(c(30, 90, 150,  210, 270, 330), byrow=TRUE, nrow=2)
#'   s <- matrix(c(0.0, 0.2, 0.4,  0.6, 0.8, 1.0), byrow=TRUE, nrow=2)
#'   v <- matrix(c(0.40, 0.44, 0.48,  0.52, 0.56, 0.60), byrow=TRUE, nrow=2)
#'   img.hsv <- createImageHSV(h, s, v)
#' @seealso \code{\link{extractHSVChannel}}
#' @export
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


#' Return a color channel (matrix) from an RGB image.
#' 
#' @param img.rgb an RGB (Reg, Green, Blue) image construct
#' @param col.channel a channel to extract:
#'                    either \code{RED}, \code{GREEN}, or \code{BLUE} constant.
#' @return matrix that contains the required channel from given RGB image.
#' @examples
#'   r <- matrix(runif(6), nrow=2)
#'   g <- matrix(runif(6), nrow=2)
#'   b <- matrix(runif(6), nrow=2)
#'   img.rgb <- createImageRGB(r, g, b)
#'   ext.r <- extractRGBChannel(img.rgb, RED)
#'   ext.g <- extractRGBChannel(img.rgb, GREEN)
#'   ext.b <- extractRGBChannel(img.rgb, BLUE)
#'   stopifnot(all(ext.r == r), all(ext.g == g), all(ext.b == b))
#' @seealso \code{\link{createImageRGB}}
#' @export
extractRGBChannel <- function(img.rgb, col.channel) img.rgb[ , , col.channel$idx]


#' Return a channel (matrix) from an HSV image.
#' 
#' @param img.hsv an HSV (Hue, Saturation, Value) image construct
#' @param hsv.channel a channel to extract:
#'   either \code{HUE}, \code{SATURATION}, or \code{VALUE} constant.
#' @return matrix that contains the required channel from given HSV image.
#' @examples
#'   saturations <- extractHSVChannel(img.hsv, SATURATION)
#' @seealso \code{\link{createImageHSV}}
#' @export
extractHSVChannel <- function(img.hsv, hsv.channel) img.hsv[ , , hsv.channel$idx]


#' Return a channel (matrix) from an XYZ image.
#' 
#' @param img.xyz an XYZ image construct (CIE 1931 XYZ color space)
#' @param xyz.channel a channel to extract:
#'   either \code{CIEXYZ_X} ("red--green"), \code{CIEXYZ_Y} ("luminosity"), or
#'   \code{CIEXYZ_Z} ("blue--yellow") constant.
#' @return matrix that contains the required channel from given XYZ image.
#' @examples
#'   z.blueyellow <- extractXYZChannel(img.xyz, CIEXYZ_Z)
#' @seealso \code{\link{toXYZ}}
#' @export
extractXYZChannel <- function(img.xyz, xyz.channel) img.xyz[ , , xyz.channel$idx]


#' Return a channel (matrix) from an xyY image.
#' 
#' @param img.xyY an xyY image construct (CIE xyY color space)
#' @param xyY.channel a channel to extract:
#'   either \code{CIExyY_x} ("color coordinate"),
#'   \code{CIExyY_y} ("luminosity"), or
#'   \code{CIExyY_Y} ("color coordinate") constant.
#' @return matrix that contains the required channel from given xyY image.
#' @examples
#'   Y.coords <- extractxyYChannel(img.xyY, CIExyY_Y)
#' @seealso \code{\link{toxyY}}
#' @export
extractxyYChannel <- function(img.xyY, xyY.channel) img.xyY[ , , xyY.channel$idx]


#' Return a channel (matrix) from a CIE LUV image.
#' 
#' @param img.luv an LUV image construct (CIE L*u*v* color space)
#' @param xyY.channel a channel to extract:
#'   either \code{CIELUV_L}, \code{CIELUV_U}, \code{CIELUV_V} constant.
#' @return matrix that contains the required channel from given LUV image.
#' @examples
#'   L <- extractLUVChannel(img.luv, CIELUV_L)
#' @seealso \code{\link{toLUV}}
#' @export
extractLUVChannel <- function(img.luv, luv.channel) img.luv[ , , luv.channel$idx]
