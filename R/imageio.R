##
## Image IO-related functions, such as reading a JPEG/PNG file.
##
## Esa Junttila, 2016-03-23

library(jpeg)  # function 'readJPG' to read JPEG files
library(png)   # function 'readPNG' to read PNG files

# Read image file (either a valid PNG or JPEG extension required)
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

  if (length(dim(img)) == 2) {  # only one value per pixel
    # Image contains only greyscale values; repeat values to create RGB channels.
    img <- array(img, dim=c(dim(img), 3))
  }
  if (dim(img)[3] == ALPHA) {
    # alpha channel detected
    if (all(img[ , , ALPHA] == 1.0)) {  # are we fully opaque? If yes, then alpha can be ignored.
      img <- img[ , , -ALPHA]  # drop the ALPHA channel without effect
    } else {
      stop(paste('Error: Alpha channel in img', filename, 'has non-opaque values (<1.0).'))
      # If this happens, we need to adjust theimage  by mixing in some background color.
    }
  }
  dimnames(img) <- list(NULL, NULL, c('Red', 'Green', 'Blue'))
  return(img)
}
