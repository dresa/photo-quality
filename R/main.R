##
## Main program: Quality measures for photographs: no-reference analysis of aesthetic features
##
## Esa Junttila, 2016-03-23


# This software has dependendies on the following libraries:
#   jpeg     -- reading jpeg image files
#   png      -- reading png image files
#   grid     -- viewing raster images on-screen
#   circular -- von Mises distribution densities
#   ggplot2  -- circular hue histogram

source('photo.R')         # photo data structures
source('imageio.R')       # read photo from a file
source('derivedimage.R')  # image conversions
source('viewer.R')        # viewing images on-screen
source('quality.R')       # computing quality measurements


##########
## MAIN ##
##########

main <- function() {
  do.view <- TRUE  #FALSE
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
  #filename <- '../examples/penguin.jpg'
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
  print(paste('Basic exposure:', basicExposureLevel(img)))
  print(paste('Basic RMS contrast:', basicRmsContrast(img)))
  print(paste('Basic interval contrast:', basicIntervalContrast(img)))
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
  if (do.view) view(toBrightRGB(img), title='Pure hues')
  print('===  END  ===')
}

main()

