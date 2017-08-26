# Run photo quality measures on photos in a training dataset.


# Key:
#   Photo -- Photo ID, like originals/12345.jpg --> 12345
#   Blur -- Blur annoyance quality (0--1)
#   MdweBlur -- MDWE horizontal blur width (smaller is better)
#   MdweGaussian -- MDWE Gaussian quality (0--1) (greater is better)
#   MdweJpeg2k -- MDWE JPEG2000 quality (0--1) (greater is better)
#   Exposure -- Basic exposure
#   RmsContrast -- Basic RMS contrast
#   IntervalContrast -- Basic interval contrast
#
# Dominant direction, spread, and portion of the dominant color.
#   DdcMean -- Dominant hue as an angle, Color dispersion(mu)
#   DdcCompactness -- Compactness of hues, larger kappa means narrower hues, Color dispersion(kappa)
#   DdcDominance -- Proportion of pixels whose hue is close to dominant, Color dispersion(pi)
#   DdcSpatial -- Spatial dispersion: mean distance of pixels that are dominant, Color dispersion(ds)
#   DdcNormSpatial -- Custom normalized spatial dispersion: dominantdistances / alldistances, Color dispersion(custom.ds)



source('photo.R')         # photo data structures
source('imageio.R')       # read photo from a file
source('derivedimage.R')  # image conversions
source('quality.R')       # computing quality measurements

PHOTO_DIR <- '../training-data/originals/'
TARGET_PATH <- '../training-data/measures.csv'

measure <- function(filename) {
  print(paste('Measuring:', filename))
  img.id <- as.numeric(unlist(strsplit(basename(filename), '\\.'))[1])
  img <- readImage(filename)
  blur <- blurAnnoyanceQuality(img, f.len=9)
  mdwe.score <- mdweVertical(img)
  exposure <- basicExposureLevel(img)
  rms.contrast <- basicRmsContrast(img)
  interval.contrast <- basicIntervalContrast(img)
  ddc <- dispersionDominantColor(toHSV(img))
  return(c(img.id, blur,
           mdwe.score[['score']], mdwe.score[['gaussian']], mdwe.score[['jpeg2k']],
           exposure, rms.contrast, interval.contrast,
           ddc$mu, ddc$kappa, ddc$pi, ddc$ds, ddc$custom.ds
  ))
}

training <- function(photo.dir, target.path) {
  photos <- list.files(photo.dir, full.names=TRUE) #[1:100]
  headers <- c('Photo', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k',
               'Exposure', 'RmsContrast', 'IntervalContrast', 'DdcMean',
               'DdcCompactness', 'DdcDominance', 'DdcSpatial', 'DdcNormSpatial')
  m <- matrix(sapply(photos, measure), ncol=length(headers), byrow=TRUE)
  df <- as.data.frame(m)
  colnames(df) <- headers
  write.table(df, file=target.path, row.names=FALSE, quote=FALSE, sep=';')
}

training(PHOTO_DIR, TARGET_PATH)

