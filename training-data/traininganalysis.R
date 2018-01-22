# Analyzing trainign results and correlations
#
# Now around 0.3 correlation after linear modeling. Should be able to do better.
# Add more measures:
# - photo saturation, average/median?
# - dimensions or photo aspect ratio
# - photo size (+ or x)
# - fine-tune existing measures: exposure and jpeg2k, DDC compacness (Inf-->BW?)
# - intervalcontrast with different intervals: 95% vs ?
# - vs. standard deviation? (instead of rating?)
# - Ritendra Datta: Studying Aesthetics in Photographic Images Using a Computational Approach
#
# Esa Junttila, 2017-08-26 (originally)

dpchallenge <- read.table('dpchallenge_dataset.txt', header=TRUE)
measures <- read.table('measures_16504.csv', header=TRUE, sep=';', stringsAsFactors=FALSE)

# All photos
df <- merge(dpchallenge, measures, by='Photo')
print('All photos correlations:')
print(cor(df[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast')]))
plot(df$IntervalContrast, df$Rating)
hist(df$Blur, breaks=20)
model <- lm(df$Rating ~ df$Blur + df$MdweBlur + df$MdweGaussian + df$MdweJpeg2k + df$Exposure + df$RmsContrast + df$IntervalContrast)
modeled.ratings <- predict(model)
print(cor(df$Rating, modeled.ratings))
sapply(colnames(df), function(x) {
  plot(df[[x]], df$Rating, xlab=x)
  mask <- is.finite(df[[x]])
  abline(lm(df$Rating[mask] ~ df[[x]][mask]), col='red')
})

# Color photos
is.bw <- is.na(df$DdcCompactness) | df$DdcCompactness < 1e-5 | is.infinite(df$DdcCompactness)
colored <- df[!is.bw, ]
print('Colored correlations:')
print(cor(colored[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast', 'DdcMean', 'DdcCompactness', 'DdcDominance', 'DdcSpatial', 'DdcNormSpatial')]))
colored.model <- lm(colored$Rating ~ colored$Blur + colored$MdweBlur + colored$MdweGaussian + colored$MdweJpeg2k + colored$Exposure + colored$RmsContrast + colored$IntervalContrast + colored$DdcMean + colored$DdcCompactness + colored$DdcDominance + colored$DdcSpatial + colored$DdcNormSpatial)
print(cor(colored$Rating, predict(colored.model)))

# Black & White photos
bw <- df[is.bw, ]
print('BW correlations:')
print(cor(bw[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast')]))
bw.model <- lm(bw$Rating ~ bw$Blur + bw$MdweBlur + bw$MdweGaussian + bw$MdweJpeg2k + bw$Exposure + bw$RmsContrast + bw$IntervalContrast)
print(cor(bw$Rating, predict(bw.model)))

# Large DDC compactness:
# Photo Index DdcCompactness DdcDominance DdcSpatial DdcNormSpatial
# 305   10651 Inf      5.402542e-01 256.991718    0.846980861
# 1776  13031 Inf      1.000000e+00 123.869782    0.543859279
# 2003  13494 Inf      9.154412e-01 178.556321    0.625816707
# 2182  13796 2902.739 5.443135e-05 240.770650    0.819341477
# 2419  14929 1953.258 2.260852e-04 255.161436    0.772387298
# 3367  18204 5315.953 3.561198e-02 233.937880    0.797297346
# 5082  23811 1309.016 8.423497e-03 278.976903    0.887392972
# 7313  38474 1126.500 9.997780e-01 392.860437    1.000000000
# 8437  43380 Inf      1.000000e+00   6.389113    0.022834884
# 10496 53787 1057.933 2.500031e-05  20.897852    0.077369038
# 11260 57288 Inf      1.000000e+00 135.224579    0.444807833
# 11350 57928 Inf      1.000000e+00   6.062565    0.020950792
# 11394 58131 Inf      1.000000e+00 195.242599    0.691562631
# 11732 59663 2435.836 9.991203e-01 147.846508    1.000000000
# 11782 59912 1176.419 3.038549e-02 196.329697    0.867218015
# 12597 62443 1092.586 1.440383e-05   1.659439    0.005981826
# 12722 62900 Inf      1.000000e+00   4.148135    0.013665845
# 13067 64130 1262.327 1.351789e-03 204.791085    0.695725103
# 14677 70749 Inf      7.583333e-01 167.415594    0.555052858
# 15173 73644 Inf      2.706410e-01 223.055849    0.761133387
# 15890 76078 3045.047 3.753870e-03 190.062750    0.738364654
# 16401 79163 1188.873 5.897273e-02 306.207215    1.000000000
# 16500 79954 Inf      8.666667e-01  23.261137    0.082806309

