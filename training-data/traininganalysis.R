# Analyzing trainign results and correlations
# Esa Junttila, 2017-08-26 (originally)

dpchallenge <- read.table('dpchallenge_dataset.txt', header=TRUE)
measures <- read.table('measures_first_100.csv', header=TRUE, sep=';', stringsAsFactors=FALSE)

df <- merge(dpchallenge, measures, by='Photo')
print(cor(df[, c('Rating', 'Blur', 'MdweBlur', 'MdweGaussian', 'MdweJpeg2k', 'Exposure', 'RmsContrast', 'IntervalContrast')]))
plot(df$IntervalContrast, df$Rating)
hist(df$Blur, breaks=20)
