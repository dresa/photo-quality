##
## Image IO-related functions, such as reading a JPEG/PNG file.
##
## Esa Junttila, 2016-03-23 (originally)


library(jpeg)  # function 'readJPG' to read JPEG files
library(png)   # function 'readPNG' to read PNG files


#' Read image file from filesystem.
#' 
#' Read an image file from filesystem and construct an internal RGB
#' representation of the image.
#' Both \emph{PNG} and \emph{JPEG} encodings for image files are supported.
#' The filename must have an appropriate extension,
#' such as \code{.jpg}, \code{.png}, or \code{.jpeg}.
#' 
#' Grayscale images with only \emph{one} channel are converted
#' into \emph{three}-channel (gray) RGB images. Possible \emph{alpha}
#' channel is dropped if fully opaque, otherwise produces an error.
#' 
#' @param filename path to an image file in the filesystem,
#'   with proper filename extension.
#' @return An RGB image. The construct is an \code{m x n x 3} array,
#'   with image height \code{m}, image width \code{n}, and with
#'   three RGB channels. All values are numeric within \code{[0;1]}.
#' @examples
#'   img.rgb <- readImage('../examples/small_grid.png')
#' @seealso
#' \code{\link{extractRGBChannel}}, \code{\link{createImageRGB}}
#' @export
readImage <- function(filename) {
  ALPHA <- 4
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

  if (length(dim(img)) == 2) {  # input has only one value per pixel
    # Image is greyscale; repeat values on three channels to create RGB image.
    img <- array(img, dim=c(dim(img), 3))
  }
  if (dim(img)[3] == ALPHA) {  # extra "alpha" channel detected?
    if (all(img[ , , ALPHA] == 1.0)) {  # if fully opaque, then ignore alpha.
      img <- img[ , , -ALPHA]  # drop the ALPHA channel without effect
    } else {
      stop(paste('Error: Alpha channel in img', filename, 'has non-opaque values (<1.0).'))
      # If this happens, we need to adjust the image by mixing in
      # some background color. Choice of background color is open.
    }
  }
  # Give names to channel dimension (three RGB color channels),
  # but leave row and column dimensions nameless.
  dimnames(img) <- list(NULL, NULL, c('Red', 'Green', 'Blue'))
  return(img)
}
