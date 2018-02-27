# Analyzing trainign results and correlations
#
# Now correlations around 0.59 (for B&W) and 0.46 (for colored images) after linear modeling. Should be able to do better.
# Add more measures:
# - fine-tune existing measures: exposure and jpeg2k, DDC compacness (Inf-->BW?)
# - intervalcontrast with different intervals: 95% vs ?
# - vs. standard deviation? (instead of rating?). The avg rating is not objective truth; the distribution of ratings is.
# - other scientific papers
# - remove unnecessary measures that have power for predictions
#
# Esa Junttila, 2017-08-26 (originally)

library(randomForest)

dpchallenge <- read.table('dpchallenge_dataset.txt', header=TRUE)
measures <- read.table('measures_datta_16504.csv', header=TRUE, sep=';', stringsAsFactors=FALSE)

# All photos
df <- merge(dpchallenge, measures, by='Photo')

# Fix problems for data analysis:
# * Use zero hue (=red) for white, black, and non-existent segments.
# * Use zero saturations and value (=non-saturated black) for non-existent segments.
# * Use zero relative size for non-existent patches
patch.hsv.size.columns <- c('AvgHuePatch1D26', 'AvgHuePatch2D27', 'AvgHuePatch3D28', 'AvgHuePatch4D29', 'AvgHuePatch5D30',
                 'AvgSatPatch1D31', 'AvgSatPatch2D32', 'AvgSatPatch3D33', 'AvgSatPatch4D34', 'AvgSatPatch5D35',
                 'AvgValPatch1D36', 'AvgValPatch2D37', 'AvgValPatch3D38', 'AvgValPatch4D39', 'AvgValPatch5D40',
                 'RelSizePatch1D41', 'RelSizePatch2D42', 'RelSizePatch3D43', 'RelSizePatch4D44', 'RelSizePatch5D45')
for (pc in patch.hsv.size.columns) df[is.na(df[[pc]]), pc] <- 0

# Position code 22 (in the center) for non-existent segments:
patch.pos.columns <- c('PositionPatch1D48', 'PositionPatch2D49', 'PositionPatch3D50', 'PositionPatch4D51', 'PositionPatch5D52')
for (pc in patch.pos.columns) df[is.na(df[[pc]]), pc] <- 22

# Distance zero from the center for non-existent segments:
patch.dist.columns <- c('DistanceFromCenterPatch1D48E', 'DistanceFromCenterPatch2D49E', 'DistanceFromCenterPatch3D50E', 'DistanceFromCenterPatch4D51E', 'DistanceFromCenterPatch5D52E')
for (pc in patch.dist.columns) df[is.na(df[[pc]]), pc] <- 0


print('All photos correlations:')
cor.all <- cor(df[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast',
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
                 'ShapeConvexityD56'
)])
#print(cor.all)

#plot(df$IntervalContrast, df$Rating)
#hist(df$Blur, breaks=20)
model <- lm(df$Rating ~ df$Blur + df$MdweBlur + df$MdweGaussian + df$MdweJpeg2k + df$Exposure + df$RmsContrast + df$IntervalContrast +
            df$AverageIntensityD01 + df$ColorfulnessD02 + df$ColorfulnessGreyD02E +
            df$AverageCentralIntensityD07 +
            df$TextureValLevel1D16 + df$TextureValLevel2D17 + df$TextureValLevel3D18 + 
            df$TextureValSumD21 +
            df$SizeFeatureD22 + df$AspectRatioD23 + df$AspectRatioMaxD23E + 
            df$NumLargePatchesD24 + df$NumClustersD25 +
            df$AvgValPatch1D36 + df$AvgValPatch2D37 + df$AvgValPatch3D38 + df$AvgValPatch4D39 + df$AvgValPatch5D40 +
            df$RelSizePatch1D41 + df$RelSizePatch2D42 + df$RelSizePatch3D43 + df$RelSizePatch4D44 + df$RelSizePatch5D45 +
            df$PositionPatch1D48 + df$PositionPatch2D49 + df$PositionPatch3D50 + df$PositionPatch4D51 + df$PositionPatch5D52 +
            df$DistanceFromCenterPatch1D48E + df$DistanceFromCenterPatch2D49E + 
            df$DistanceFromCenterPatch3D50E + df$DistanceFromCenterPatch4D51E + df$DistanceFromCenterPatch5D52E +
            df$DepthOfFieldValD55 + 
            df$ShapeConvexityD56)
modeled.ratings <- predict(model)
print(cor(df$Rating, modeled.ratings))
sapply(colnames(df), function(x) {
  plot(df[[x]], df$Rating, xlab=x)
  mask <- is.finite(df[[x]])
  abline(lm(df$Rating[mask] ~ df[[x]][mask]), col='red')
})


# Color photos
##############
is.bw <- is.na(df$DdcCompactness) | df$DdcCompactness < 1e-5 | is.infinite(df$DdcCompactness) | df$AverageCentralSaturationD06 < 1e-5
colored <- df[!is.bw, ]
print('Colored correlations:')
cor.color <- cor(colored[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast',
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
                      'DdcMean', 'DdcCompactness', 'DdcDominance', 'DdcSpatial', 'DdcNormSpatial')])
#print(cor.color)
colored.model <- lm(colored$Rating ~ colored$Blur + colored$MdweBlur + colored$MdweGaussian + colored$MdweJpeg2k + colored$Exposure + colored$RmsContrast + colored$IntervalContrast +
                      colored$AverageIntensityD01 + colored$ColorfulnessD02 + colored$ColorfulnessGreyD02E +
                      colored$AverageSaturationD03 + colored$AverageHueD04 +
                      colored$AverageCentralHueD05 + colored$AverageCentralSaturationD06 + colored$AverageCentralIntensityD07 +
                      colored$TextureHueLevel1D10 + colored$TextureHueLevel2D11 + colored$TextureHueLevel3D12 + 
                      colored$TextureSatLevel1D13 + colored$TextureSatLevel2D14 + colored$TextureSatLevel3D15 + 
                      colored$TextureValLevel1D16 + colored$TextureValLevel2D17 + colored$TextureValLevel3D18 + 
                      colored$TextureHueSumD19 + colored$TextureSatSumD20 + colored$TextureValSumD21 +
                      colored$SizeFeatureD22 + colored$AspectRatioD23 + colored$AspectRatioMaxD23E + 
                      colored$NumLargePatchesD24 + colored$NumClustersD25 +
                      colored$AvgHuePatch1D26 + colored$AvgHuePatch2D27 + colored$AvgHuePatch3D28 + colored$AvgHuePatch4D29 + colored$AvgHuePatch5D30 +
                      colored$AvgSatPatch1D31 + colored$AvgSatPatch2D32 + colored$AvgSatPatch3D33 + colored$AvgSatPatch4D34 + colored$AvgSatPatch5D35 +
                      colored$AvgValPatch1D36 + colored$AvgValPatch2D37 + colored$AvgValPatch3D38 + colored$AvgValPatch4D39 + colored$AvgValPatch5D40 +
                      colored$RelSizePatch1D41 + colored$RelSizePatch2D42 + colored$RelSizePatch3D43 + colored$RelSizePatch4D44 + colored$RelSizePatch5D45 +
                      colored$PositionPatch1D48 + colored$PositionPatch2D49 + colored$PositionPatch3D50 + colored$PositionPatch4D51 + colored$PositionPatch5D52 +
                      colored$DistanceFromCenterPatch1D48E + colored$DistanceFromCenterPatch2D49E + 
                      colored$DistanceFromCenterPatch3D50E + colored$DistanceFromCenterPatch4D51E + colored$DistanceFromCenterPatch5D52E +
                      colored$DepthOfFieldHueD53 + colored$DepthOfFieldSatD54 + colored$DepthOfFieldValD55 + 
                      colored$ShapeConvexityD56 +
                      colored$DdcMean + colored$DdcCompactness + colored$DdcDominance + colored$DdcSpatial + colored$DdcNormSpatial)
print(cor(colored$Rating, predict(colored.model)))
#write.table(cor(colored), file='correlations_1000_colored.csv')


# Black & White photos
######################
bw <- df[is.bw, ]
print('BW correlations:')
cor.bw <- cor(bw[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast',
                 'AverageIntensityD01', 'ColorfulnessD02', 'ColorfulnessGreyD02E',
                 'AverageCentralIntensityD07',
                 'TextureValLevel1D16', 'TextureValLevel2D17', 'TextureValLevel3D18', 
                 'TextureValSumD21', 
                 'SizeFeatureD22', 'AspectRatioD23', 'AspectRatioMaxD23E', 
                 'NumLargePatchesD24', 'NumClustersD25',
                 'AvgValPatch1D36', 'AvgValPatch2D37', 'AvgValPatch3D38', 'AvgValPatch4D39', 'AvgValPatch5D40', 
                 'RelSizePatch1D41', 'RelSizePatch2D42', 'RelSizePatch3D43', 'RelSizePatch4D44', 'RelSizePatch5D45', 
                 'PositionPatch1D48', 'PositionPatch2D49', 'PositionPatch3D50', 'PositionPatch4D51', 'PositionPatch5D52', 
                 'DistanceFromCenterPatch1D48E', 'DistanceFromCenterPatch2D49E',
                 'DistanceFromCenterPatch3D50E', 'DistanceFromCenterPatch4D51E', 'DistanceFromCenterPatch5D52E', 
                 'DepthOfFieldValD55', 
                 'ShapeConvexityD56'
                 )])
#print(cor.bw)
bw.model <- lm(bw$Rating ~ bw$Blur + bw$MdweBlur + bw$MdweGaussian + bw$MdweJpeg2k + bw$Exposure + bw$RmsContrast + bw$IntervalContrast +
                 bw$AverageIntensityD01 + bw$ColorfulnessD02 + bw$ColorfulnessGreyD02E +
                 bw$AverageCentralIntensityD07 +
                 bw$TextureValLevel1D16 + bw$TextureValLevel2D17 + bw$TextureValLevel3D18 + 
                 bw$TextureValSumD21 +
                 bw$SizeFeatureD22 + bw$AspectRatioD23 + bw$AspectRatioMaxD23E + 
                 bw$NumLargePatchesD24 + bw$NumClustersD25 +
                 bw$AvgValPatch1D36 + bw$AvgValPatch2D37 + bw$AvgValPatch3D38 + bw$AvgValPatch4D39 + bw$AvgValPatch5D40 +
                 bw$RelSizePatch1D41 + bw$RelSizePatch2D42 + bw$RelSizePatch3D43 + bw$RelSizePatch4D44 + bw$RelSizePatch5D45 +
                 bw$PositionPatch1D48 + bw$PositionPatch2D49 + bw$PositionPatch3D50 + bw$PositionPatch4D51 + bw$PositionPatch5D52 +
                 bw$DistanceFromCenterPatch1D48E + bw$DistanceFromCenterPatch2D49E + 
                 bw$DistanceFromCenterPatch3D50E + bw$DistanceFromCenterPatch4D51E + bw$DistanceFromCenterPatch5D52E +
                 bw$DepthOfFieldValD55 + 
                 bw$ShapeConvexityD56)
print(cor(bw$Rating, predict(bw.model)))


# Learn a random-forest model for colored photos.
set.seed(1)
valid.ratings <- colored$Rating > 0  # skip photos with no ratings
dataset <- colored[valid.ratings, ]
n <- nrow(dataset)
train.set <- sample(1:n, 0.6*n)
test.set <- setdiff(1:n, train.set)
variables <- subset(dataset, select=-c(Rating, Num, Num1, Num2, Num3, Num4, Num5, Num6, Num7, Num8, Num9, Num10, Photo, Index))
rf <- randomForest(variables[train.set, ], dataset$Rating[train.set])
print(rf)
plot(rf)
pred.ratings <- predict(rf, variables[test.set, ])
real.ratings <- dataset$Rating[test.set]
plot(real.ratings, pred.ratings, pch='.', col='darkblue', title='Average ratings vs. predictions')
abline(lm(real.ratings ~ pred.ratings), col="darkred")
rf.corr <- cor(real.ratings, pred.ratings)
print(paste('Predicted random-forest corr for colored:', rf.corr))

refDeviations <- function(df) {
  mu <- df$Rating
  N <- df$Num
  cols <- paste0('Num', 1:10)
  counts <- sapply(cols, function(column) df[ , column])
  diffs <- sapply(1:10, function(x) x - mu)
  stdev <- sqrt(rowSums(counts * diffs^2) / (N - 1))
  return(stdev)
}
ref.stdev <- refDeviations(dataset[test.set, ])
ref.ratings <- rnorm(length(real.ratings), mean=real.ratings, sd=ref.stdev)
ref.corr <- cor(ref.ratings, real.ratings)
print(paste('Sampled reference corr for colored:', ref.corr))

# How much of correlation overfitting may explain (generate random data).
# Use worst-case scenario: using test data for BOTH model-training and testing.
rnd.variables <- matrix(runif(length(test.set) * ncol(measures)), nrow=length(test.set))
rnd.target <- runif(length(test.set))
linmod <- lm(rnd.target ~ rnd.variables)
rnd.corr <- cor(rnd.target, predict(linmod, data.frame(rnd.variables)))
print(paste('Randomized correlation (max explained by overfitting):', rnd.corr))

# Example results:
#   Predicted random-forest corr for colored: 0.423091314162033
#   Sampled reference corr for colored: 0.467509877868063
#   Randomized correlation (max explained by overfitting): 0.105893844003516
# => We are quite close to a level of how a human predicts, like 78 %.
