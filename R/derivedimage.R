##
## Derived images, processed from an original image.
##
## Esa Junttila, 2016-03-23


source('R/convolutionfilter.R')  # filters (kernels) for 2D convolution
source('R/imageconvolution.R')  # 2D convolution functions
source('R/colorspace.R')  # convert RGB to HSV


#' Otsu threshold for image segmentation.
#' 
#' Compute the \emph{Otsu threshold}, which gives a non-parametric threshold
#' for segmenting the colors in a greyscale image into black and white.
#' 
#' @param freq frequencies of color intensity classes in a gray-scale image.
#' @return Returns a threshold index \code{k}: classify grayscale colors
#' from indices \code{[1,2,..,k]} into \emph{black} and the rest as \emph{white}.
#' @examples
#' x <- sample(1 - exp(-3*abs(seq(0,1,0.001))))
#' h <- hist(x, breaks=10)
#' all(h$count == c(36,39,44,52,61,74,96,135,231,233))
#' k <- otsuThreshold(h$count)
#' k == 6
#' threshold.value <- h$breaks[k+1]
#' classes <- x <= threshold.value
#' sum(classes) == sum(h$count[1:6])
#' 
#' @seealso 
#' Nobuyuki Otsu, "A threshold selection method from gray-level histograms".
#' \emph{IEEE Trans. Sys., Man., Cyber. 9 (1): 62--66. (1979).}
#' 
#' Link to the article:
#' \url{https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4310076}.
#' 
#' Ported from \url{http://clickdamage.com/sourcecode/code/otsuThreshold.m}.
#' @export
otsuThreshold <- function(freq) {
  n <- length(freq)
  idxs <- 1:n
  cumus <- cumsum(freq)
  weighted <- cumsum(freq*idxs)
  u <- weighted / cumus
  tmp <- cumus[n] - cumus
  v <- (weighted[n] - weighted) / (tmp + (tmp==0))
  f <- 1.0 * cumus * (cumus[n] - cumus) * (u-v)^2
  i <- which.max(f)
  return(i)
}



#' @rdname imageSharpBlur
#' @name imageSharpBlur
#' @aliases meanBlurred
#' @aliases gaussianBlurred
#' @aliases adjustBlur
#' @aliases sharpAmount
#' @aliases blurAmount
#' @title Blurring and sharpening methods.
#' @description Methods for measuring different types of blur and sharpening
#' and adding them to images.
#' @details On edge pixels, we use the "replication" method in 2D convolution.
#' @param img image construct of RGB type
#' @param size size of a filter window around pixel (odd integer, default 3)
#' @param sigma standard deviation for a Gaussian blur filter to weight
#' nearby pixels (units as pixel widths)
#' @param blur.factor Adjustment to the impact of default blur effect.
#' Negative values add blur, positive values reduce blur up to 1.0;
#' larger than 1.0 overcompensate blue removal.
#' \itemize{
#'   \item \code{< -1}: add extra blur to the image (on top of default blur)
#'   \item \code{= -1}: result is a blurred image through default blur filter
#'   \item \code{=  0}: result is identical to original image
#'   \item \code{= +1}: blur removed (with a blur filter)
#'   \item \code{> +1}: overcompensated blur removal
#' }
#' @examples
#' set.seed(1)
#' d <- c(4,12,3)
#' rgb.img <- array(sample(0:255, prod(d), replace=TRUE), dim=d) / 255
#' 
#' @seealso \code{\link{createImageRGB}}
NULL


#' \code{meanBlurred} adds blur to an image by means of 2D convolution:
#' assign to each pixel a color that is a linear average of the colors within
#' a \emph{size x size} window around the pixel.
#' @return \code{meanBlurred} returns an adjusted image, with blur added
#' via linear filter.
#' @examples
#' blurred <- meanBlurred(rgb.img)
#' shift <- blurred - rgb.img
#' all(abs(shift[ , 2, 1] - c(0.19084967, -0.43747277, -0.41176471, -0.08278867)) < 1e-6)
#' 
#' @rdname imageSharpBlur
#' @export
meanBlurred <- function(img, size=3) {
  f <- meanFilter(size)
  blurred <- imfilter(img, f, 'replicate')
  return(blurred)
}


#' \code{gaussianBlurred} adds blur to an image by means of weighted
#' 2D convolution: assign to each pixel a color that is a Gaussian-weighted
#' average of the colors of other pixels aroud it, weighted by Gaussian
#' distribution with \code{sigma} standard deviation.
#' @return \code{gaussianBlurred} returns an adjusted image, with blur added
#' via Guassian filter.
#' @examples
#' blurred <- gaussianBlurred(rgb.img, 1.0)
#' shift <- blurred - rgb.img
#' all(abs(shift[ , 2, 1] - c(0.20684652, -0.39785871, -0.36144287, -0.05093839)) < 1e-6)
#' 
#' @rdname imageSharpBlur
#' @export
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


#' \code{adjustBlur} adds a custom-amplified amount of blur -- or blur
#' removal -- to an image. The baseline is the \code{meanBlurred} blur
#' function, using a linear blur filter with given window \code{size}.
#' @return \code{adjustBlur} returns an adjusted image, that might have
#' blur added, blur removed, sharpened, or nothing at all.
#' @examples
#' half.blurred <- adjustBlur(rgb.img, -0.5)
#' b.shift <- half.blurred - rgb.img
#' all(abs(b.shift[ , 2, 1] - c(0.09542484, -0.21873638, -0.20588235, -0.04139434)) < 1e-6)
#' sharpened <- adjustBlur(rgb.img, +1.0)
#' s.shift <- sharpened - rgb.img
#' all(abs(s.shift[ , 2, 1] - c(-0.19084967, 0.10196078, 0.05490196, 0.08278867)) < 1e-6)
#' 
#' @rdname imageSharpBlur
#' @export
adjustBlur <- function(img, blur.factor, size=3) {
  blurred <- meanBlurred(img, size)
  return(pmin(pmax(img + (img-blurred)*blur.factor, 0), 1))
}


#' \code{sharpAmount} returns an "intensity image" that shows for
#' each pixel (per channel) the intensity of sharpening it requires, as
#' determined by baseline \code{meanBlurred} blur function. The blur function
#' uses a linear blur filter with given window \code{size}.
#' @return \code{sharpAmount} returns an "intensity image" of sharpening
#' intensity; that is, positive numbers within \code{[0;1]} about the intensity
#' of sharpening needed for each pixel. Note that sharpening intensity may not
#' be fully added to pixel values because they are bounded by minimum or
#' maximum values.
#' @examples
#' sharp.intensity <- sharpAmount(rgb.img)
#' all(abs(sharp.intensity[ , 2, 1] - c(0.19084967, 0.43747277, 0.41176471, 0.08278867)) < 1e-6)
#' 
#' @rdname imageSharpBlur
#' @export
sharpAmount <- function(img, size=3) {
  blurred <- meanBlurred(img, size)
  return(abs(img - blurred))
}

#' \code{blurAmount} returns an "intensity image" that shows for
#' each pixel (per channel) the intensity of blur that is present.
#' This is the same as \code{1 - sharpAmount}.
#' @return \code{blurAmount} returns an "intensity image" of blur effects
#' ; that is, positive numbers within \code{[0;1]} about the intensity
#' of blur that occurs for each pixel.
#' @examples
#' blur.intensity <- blurAmount(rgb.img)
#' all(abs(blur.intensity[ , 2, 1] - (1 - c(0.19084967, 0.43747277, 0.41176471, 0.08278867))) < 1e-6)
#' 
#' @rdname imageSharpBlur
#' @export
blurAmount <- function(img, size=3) {
  return(1 - sharpAmount(img))
}


#' Luminance of an image.
#' 
#' Estimate the perceptual luminance of each pixel in an RGB image.
#' 
#' @details
#' We use the \emph{Y} variable in CIE XYZ color space as a proxy
#' for human perception of luminance.
#' 
#' \preformatted{
#' CIE 1931 linear luminance Y
#' ---------------------------
#' RGB<->XYZ conversion in sRGB space using reference white D65:
#'   |X|       |R|              |X|   |R|
#'   |Y| = M * |G|        M(-1)*|Y| = |G|
#'   |Z|       |B|              |Z|   |B|
#'   
#' Full matrices:
#'         |  0.4124564  0.3575761  0.1804375 |
#' M =     |  0.2126729  0.7151522  0.0721750 |
#'         |  0.0193339  0.1191920  0.9503041 |
#'
#'         |  3.2404542 -1.5371385 -0.4985314 |
#' M(-1) = | -0.9692660  1.8760108  0.0415560 |
#'         |  0.0556434 -0.2040259  1.0572252 |
#' }
#' 
#' In particular \emph{Y} coefficients for \emph{R}, \emph{G},
#' and \emph{B} are 0.2126729, 0.7151522, and 0.0721750.
#' 
#' Many other choices can be found here:
#' \url{http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html}
#' 
#' I'm not sure whether gamma correction should be applied here.
#' 
#' @param img an RGB image construct
#' @return Returns a matrix (a "luminance channel") that contains luminance
#' values \code{[0;1]} for each pixel.
#' @examples
#' set.seed(1)
#' d <- c(4,12,3)
#' rgb.img <- array(sample(0:255, prod(d), replace=TRUE), dim=d) / 255
#' lum <- luminance(rgb.img)       # per-pixel luminance
#' (mean(lum) - 0.5102388) < 1e-6  # average luminance in image
#' 
#' @seealso \code{\link{createImageRGB}}
#' @export
luminance <- function(img) 0.212673*img[,,RED$idx] + 0.715152*img[,,GREEN$idx] + 0.072175*img[,,BLUE$idx]


#' Bright-colored image with pure hues.
#' 
#' Transform the RGB colors of an image so that resulting image contains
#' only pure bright hues with maximum saturation.
#' @details The pixels for which there is no \emph{hue} information,
#' such as black pixels, use the color specified in argument \code{na.color}
#' (default is "\emph{gray}").
#' 
#' Having maximum saturation means that for each pixel in the returned image
#' at least one of \emph{R}, \emph{G}, or \emph{B} is zero.
#' @param img.rgb an RGB image construct
#' @param na.color RGB color whose information is used on those
#' pixels that are missing hues, such as \emph{black}.
#' Example values: \code{"gray"} or \code{"#FFCCAA"}.
#' @return Returns an RGB image that shares the hues with original image,
#' but has maximum saturation in each pixel.
#' @examples
#' set.seed(1)
#' d <- c(4,12,3)
#' rgb.img <- array(sample(0:255, prod(d), replace=TRUE), dim=d) / 255
#' bright <- toBrightRGB(rgb.img)
#' all(abs(bright[ , 2, 1] - c(0, 1, 1, 0.628821)) < 1e-6)
#' 
#' @seealso \code{\link{createImageRGB}}
#' @export
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

