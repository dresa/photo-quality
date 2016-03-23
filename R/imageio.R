##
## Image IO-related functions, such as reading a JPEG/PNG file.
##
## Esa Junttila, 2016-03-23

library(jpeg)  # function 'readJPG' to read JPEG files
library(png)   # function 'readPNG' to read PNG files

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
