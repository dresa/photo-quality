##
## Viewer for images and photos, on.screen
##
## Esa Junttila, 2016-03-23


library(grid)  # function 'grid.raster' to view images
library(ggplot2)  # hue-usage in polar coordinates
library(circular) # von Mises distribution densities

source('R/photo.R')  # use extractHSVChannel
source('R/common.R')  # use degreeToRadian


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
  # Replace with a window function call of your own OS:
  #dev.new(width=win.w, height=win.h, title=title)
  windows(width=win.w, height=win.h, title=title)
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
  mask <- !is.na(h.raw)
  if (sum(mask) <= 0.01*length(hues)) {
    print('Skipping DDC hue plot because the image is (almost) greyscale.')
    return()
  }
  clean.hues <- h.raw[mask]
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

