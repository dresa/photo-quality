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
source('qualitydatta.R')  # computing "Datta" collection of quality measurements

PHOTO_DIR <- '../training-data/originals/'
TARGET_PATH <- '../training-data/measures.csv'

measure <- function(filename) {
  print(paste('Measuring:', filename))
  img.id <- as.numeric(unlist(strsplit(basename(filename), '\\.'))[1])
  img <- readImage(filename)
  img.hsv <- toHSV(img)
  blur <- blurAnnoyanceQuality(img, f.len=9)
  mdwe.score <- mdweVertical(img)
  exposure <- basicExposureLevel(img)
  rms.contrast <- basicRmsContrast(img)
  interval.contrast <- basicIntervalContrast(img)
  avg.intensity <- avgIntensity(img.hsv)
  colorfulness.score <- colorfulness(img)
  colorfulness.grey.score <- colorfulness(img, 'grey', n=6)
  avg.saturation <- avgSaturation(img.hsv)
  avg.hue <- avgHue(img.hsv)
  avg.central.hue <- avgCentralHue(img.hsv)
  avg.central.saturation <- avgCentralSaturation(img.hsv)
  avg.central.intensity <- avgCentralIntensity(img.hsv)
  tx <- texture(img.hsv)
  size.feat <- sizeFeature(img)
  aspect.ratio <- aspectRatioFeature(img)
  datta.seg <- regionCompositionFeatures(img)
  dof.feat <- dof(img.hsv)
  ddc <- dispersionDominantColor(img.hsv)

  # Shorthand convenience variables:
  dsap <- datta.seg$avg.patch.hsv
  dsrp <- datta.seg$rel.patch.sizes
  dssp <- datta.seg$segment.positions
  dssd <- datta.seg$segment.distances

  return(c(img.id, blur,
           mdwe.score[['score']], mdwe.score[['gaussian']], mdwe.score[['jpeg2k']],
           exposure, rms.contrast, interval.contrast,
           avg.intensity, colorfulness.score, colorfulness.grey.score, avg.saturation, avg.hue,
           avg.central.hue, avg.central.saturation, avg.central.intensity,
           tx$hue.1, tx$hue.2, tx$hue.3,
           tx$sat.1, tx$sat.2, tx$sat.3,
           tx$val.1, tx$val.2, tx$val.3,
           tx$hue.sum, tx$sat.sum, tx$val.sum,
           size.feat, aspect.ratio[1], aspect.ratio[2],
           datta.seg$num.large.patches, datta.seg$num.clusters,
           dsap[1], dsap[2], dsap[3], dsap[4], dsap[5],       # average Hues for each of largest five segments
           dsap[6], dsap[7], dsap[8], dsap[9], dsap[10],      # average Saturations for each of largest five segments
           dsap[11], dsap[12], dsap[13], dsap[14], dsap[15],  # average Values for each of largest five segments
           dsrp[1], dsrp[2], dsrp[3], dsrp[4], dsrp[5],       # relative sizes of largest five segments
           dssp[1], dssp[2], dssp[3], dssp[4], dssp[5],       # position codes of largest five segments (11,12,13,21,...,32,33)
           dssd[1], dssd[2], dssd[3], dssd[4], dssd[5],       # distance from center (largest five segments)
           dof.feat$hue.dof, dof.feat$sat.dof, dof.feat$val.dof,
           datta.seg$shape.convexity,
           ddc$mu, ddc$kappa, ddc$pi, ddc$ds, ddc$custom.ds
  ))
}


training <- function(photo.dir, target.path) {
  photos <- list.files(photo.dir, full.names=TRUE)[1:3]
  headers <- c('Photo', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k',
               'Exposure', 'RmsContrast', 'IntervalContrast',
               'AverageIntensityD01', 'ColorfulnessD02', 'ColorfulnessGreyD02E',
               'AverageSaturationD03', 'AverageHueD04',
               'AverageCentralHueD05', 'AverageCentralSaturationD06', 'AverageCentralIntensityD07',
               'TextureHueLevel1D10', 'TextureHueLevel2D11', 'TextureHueLevel3D12', 
               'TextureSatLevel1D13', 'TextureSatLevel2D14', 'TextureSatLevel3D15', 
               'TextureValLevel1D16', 'TextureValLevel2D17', 'TextureValLevel3D18', 
               'TextureHueSumD19', 'TextureSatSumD20', 'TextureValSumD21', 
               'SizeFeatureD22', 'AspectRatioD23', 'AspectRatioMaxD23E', 
               'NumLargePatchesD24', 'NumClustersD25',
               'AvgHuePatch1D26', 'AvgHuePatch2D27', 'AvgHuePatch3D28', 'AvgHuePatch4D29', 'AvgHuePatch5D30', 
               'AvgSatPatch1D31', 'AvgSatPatch2D32', 'AvgSatPatch3D33', 'AvgSatPatch4D34', 'AvgSatPatch5D35', 
               'AvgValPatch1D36', 'AvgValPatch2D37', 'AvgValPatch3D38', 'AvgValPatch4D39', 'AvgValPatch5D40', 
               'RelSizePatch1D41', 'RelSizePatch2D42', 'RelSizePatch3D43', 'RelSizePatch4D44', 'RelSizePatch5D45', 
               'PositionPatch1D48', 'PositionPatch2D49', 'PositionPatch3D50', 'PositionPatch4D51', 'PositionPatch5D52', 
               'DistanceFromCenterPatch1D48E', 'DistanceFromCenterPatch2D49E',
               'DistanceFromCenterPatch3D50E', 'DistanceFromCenterPatch4D51E', 'DistanceFromCenterPatch5D52E', 
               'DepthOfFieldHueD53', 'DepthOfFieldSatD54', 'DepthOfFieldValD55', 
               'ShapeConvexityD56',
               'DdcMean', 'DdcCompactness', 'DdcDominance', 'DdcSpatial', 'DdcNormSpatial')
  m <- matrix(sapply(photos, measure), ncol=length(headers), byrow=TRUE)
  df <- as.data.frame(m)
  colnames(df) <- headers
  write.table(df, file=target.path, row.names=FALSE, quote=FALSE, sep=';')
}

training(PHOTO_DIR, TARGET_PATH)

