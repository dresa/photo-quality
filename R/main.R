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
#   emdist
#   waveslim
#   NbClust
#   mmand

source('R/photo.R')         # photo data structures
source('R/imageio.R')       # read photo from a file
source('R/derivedimage.R')  # image conversions
source('R/viewer.R')        # viewing images on-screen
source('R/quality.R')       # computing quality measurements
source('R/qualitydatta.R')  # quality measurements by Datta et al


##########
## MAIN ##
##########

main <- function() {
  do.view <- FALSE #TRUE  #FALSE
  print('=== START ===')
  #filename <- 'examples/small_grid.png'
  #filename <- 'examples/blue_shift.png'
  #filename <- 'examples/no_shift.png'
  #filename <- 'examples/niemi.png'
  #filename <- 'examples/sharp_or_blur.png'  # Blur annoyance quality (1--5): 1.17416513963911"
  #filename <- 'examples/K5_10994.JPG'
  #filename <- 'examples/green_grass_blue_sky.png'
  #filename <- 'examples/dark_city.png'
  #filename <- 'examples/violetred.png'
  #filename <- 'examples/bluehue.png'
  filename <- 'examples/puffin.jpg'
  #filename <- 'examples/temple_set/temple-a-original.png'
  #filename <- 'examples/temple_set/temple-b-blue.png'
  #filename <- 'examples/temple_set/temple-c-cyan.png'
  #filename <- 'examples/temple_set/temple-d-yellow.png'
  #filename <- 'examples/temple_set/temple-e-magenta.png'
  #filename <- 'examples/temple_set/temple-f-red.png'
  #filename <- 'examples/temple_set/temple-g-green.png'
  #filename <- 'examples/temple_set/temple-h-noise.png'
  #filename <- 'examples/temple_set/temple-i-contrast.png'
  #filename <- 'examples/temple_set/temple-j-colornoise.png'
  #filename <- 'examples/temple_set/temple-k-gaussianblur.png'
  #filename <- 'examples/almost_black.png'
  #filename <- 'examples/grainy.jpg'
  #filename <- 'examples/uniform-buckets.png'
  #filename <- 'examples/many_colors.png'
  #filename <- 'examples/colorfulness-test.png'
  #filename <- 'examples/bluehue.png'
  #filename <- 'examples/pure-red.png'
  #filename <- 'examples/flower.jpg'

  # Example images from the article do not behave like the authors claim:
  #filename <- 'examples/datta-colorfulness-high-1.png'
  #filename <- 'examples/datta-colorfulness-high-2.png'
  #filename <- 'examples/datta-colorfulness-low-1.png'
  #filename <- 'examples/datta-colorfulness-low-2.png'

  #filename <- 'training-data/originals/33761.jpg'  # 10338, 11402, 10170, 10048

  img <- readImage(filename)
  img.hsv <- toHSV(img)
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
  print(paste('Average Intensity (Datta 1):', avgIntensity(img.hsv)))
  print(paste('Colorfulness (Datta 2):', colorfulness(img)))
  print(paste('Colorfulness-Grey (Datta 2 grey, n=6):', colorfulness(img, 'grey', n=6)))
  print(paste('Average saturation (Datta 3):', avgSaturation(img.hsv)))
  print(paste('Average hue (Datta 4):', avgHue(img.hsv)))
  print(paste('Average central hue (Datta 5):', avgCentralHue(img.hsv)))
  print(paste('Average central saturation (Datta 6):', avgCentralSaturation(img.hsv)))
  print(paste('Average central intensity (Datta 7):', avgCentralIntensity(img.hsv)))
  tx <- texture(img.hsv)
  print(paste('Texture, hue        (Datta 10,11,12,19):', paste(tx$hue.1, tx$hue.2, tx$hue.3, paste('sum', tx$hue.sum), sep=', ')))
  print(paste('Texture, saturation (Datta 13,14,15,20):', paste(tx$sat.1, tx$sat.2, tx$sat.3, paste('sum', tx$sat.sum), sep=', ')))
  print(paste('Texture, value      (Datta 16,17,18,21):', paste(tx$val.1, tx$val.2, tx$val.3, paste('sum', tx$val.sum), sep=', ')))
  print(paste('Size feature (Datta 22):', paste(sizeFeature(img), collapse=', ')))
  print(paste('Aspect ratio (Datta 23):', paste(aspectRatioFeature(img), collapse=', ')))
  datta.seg <- regionCompositionFeatures(img)
  print(paste('Number of large patches (Datta 24):', datta.seg$num.large.patches))
  print(paste('Number of clusters (Datta 25):', datta.seg$num.clusters))
  print(paste('Average patch HSV values (Datta 26--40):', paste(datta.seg$avg.patch.hsv, collapse=', ')))
  print(paste('Relative patch sizes (Datta 41--45):', paste(datta.seg$rel.patch.sizes, collapse=', ')))
  # Excluding Datta 46 & 47 measures, as they make no sense.
  print(paste('Segment position codes (Datta 48--52):', paste(datta.seg$segment.positions, collapse=', ')))
  print(paste('Segment distances from center (Esa proxy 48--52):', paste(datta.seg$segment.distances, collapse=', ')))
  print(paste('Depth-of-field for H, S, V (Datta 53--55):', paste(dof(img.hsv), collapse=', ')))
  print(paste('Shape convexity feature (Datta 56):', datta.seg$shape.convexity))

  ddc <- dispersionDominantColor(img.hsv)

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
  #if (do.view) viewDDC(filename='polarcolor.png', toHSV(img), ddc$mu, ddc$kappa)
  if (do.view) view(toBrightRGB(img), title='Pure hues')
  print('===  END  ===')
}

##
## LAUNCH MAIN PROGRAM
##

#set.seed(1)

#Rprof(filename="Rprof.out", append=FALSE, line.profiling=TRUE, interval=0.01)  # Start profiling

#main()

#Rprof(NULL)  # End profiling

# Use summaryRprof('Rprof.out', lines='show') to view the profiling results.
