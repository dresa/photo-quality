##
## Test scripts that were used for testing at some point.
## Now outdated, but stored for reference.
##
## Esa Junttila, 2016-03-23


# Measuring 2D Gaussian blur:
#img <- readImage('../examples/small_grid.png')
#s <- Sys.time()
#replicate(1000, gaussianBlurred(img, 1.0))
#print(Sys.time() - s)

#print('Adjusted blur:')
#print(blurAnnoyanceQuality(readImage('../examples/small_grid.png'), f.len=3))  # conv2 & conv2dFFT: 2.002251 quality, normalized 0.2505627.

#s <- Sys.time()
#replicate(50, adjustBlur(readImage('../examples/sharp_or_blur.png'), 1.0, 5))
#print(Sys.time() - s)
# With a very small filter matrix:
#   conv2:     1.245071 secs
#   conv2dFFT: 1.492085 secs
# I suppose it takes a larger filter matrix to pay off using FFT.

#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(c(1/3), nrow=3, ncol=1)))
#     [,1] [,2] [,3] [,4] [,5] [,6]
#[1,]    7    8    9   10   11   12
#[2,]   13   14   15   16   17   18
#[3,]   19   20   21   22   23   24
#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(c(1/3), nrow=1, ncol=3)))
#     [,1] [,2] [,3] [,4]
#[1,]    2    3    4    5
#[2,]    8    9   10   11
#[3,]   14   15   16   17
#[4,]   20   21   22   23
#[5,]   26   27   28   29

#n<-9; vals <- (1:n)/(n*(n+1)/2)
#print(conv2(matrix(1:30, nrow=5, byrow=TRUE), matrix(vals, nrow=3, ncol=3)))
#     [,1] [,2] [,3] [,4]
#[1,]  6.8  7.8  8.8  9.8
#[2,] 12.8 13.8 14.8 15.8
#[3,] 18.8 19.8 20.8 21.8

# Marziliano data (estimated from plots):
#   Gaussian comparison (x: objective, y: subjective):
#   gx <- c(4.5, 4.5, 4.6, 4.7, 4.7, 4.8, 4.9, 5.1, 5.3, 5.4, 5.9, 6.0, 6.1, 6.4, 7.2, 7.3, 7.6, 8.0, 8.4, 8.6, 9.2, 9.4, 9.6, 9.8, 10.9, 11.2, 11.3, 11.4, 13.0, 13.4)
#   gy <- c(0.1, 0.2, 0.6, 0.8, 0.6, 0.3, 0.4, 0.3, 0.7, 1.2, 2.9, 2.8, 2.3, 2.5, 3.8, 2.3, 5.1, 4.9, 4.7, 6.0, 7.3, 5.5, 6.5, 7.3, 8.7, 7.0, 8.1, 7.6, 8.9, 8.0)
#   JPEG comparison (x: objective, y: subjective):
#   jx <- c(4.5, 4.5, 4.6, 4.8, 5.1, 5.2, 5.3, 5.3, 5.7, 5.7, 5.8, 5.9, 5.9, 6.2, 6.2, 6.3, 6.6, 6.7, 6.8, 6.8, 7.0, 7.0, 7.1, 7.1, 7.2, 7.7, 7.7, 7.8, 8.4, 9.2)
#   jy <- c(0.3, 0.6, 0.1, 0.8, 2.5, 0.9, 4.5, 3.5, 1.6, 4.8, 3.6, 7.1, 5.8, 6.0, 3.3, 8.5, 5.1, 9.2, 7.1, 6.6, 7.2, 6.1, 7.4, 8.7, 8.8, 8.1, 8.1, 8.0, 8.0, 8.9)
#  My own functions: objective x -> subjective y (limited to 1--10 disturbance, higher is worse)
#  Gaussian: y = x^1.04 - 4.5
#  JPEG: y = 17.5*(x-3)^(0.3)-19.5

# Von Mises test data:
# http://www.stat.sfu.ca/content/dam/sfu/stat/alumnitheses/MiscellaniousTheses/Bentley-2006.pdf
#   degrees <- c(0,0,0,15,45,68,100,110,113,135,135,140,140,155,165,165,169,180,180,180,180,180,180,180,189,206,209,210,214,215,225,226,230,235,245,250,255,255,260,260,260,260,270,270)
#   rads <- degrees / 360 * (2*pi)
#   mu <- atan2(sum(sin(rads)), sum(cos(rads))) %% (2*pi)
#   kappa <- est.kappa(rads)
# Compare:
#   hist(degrees)
#   hist(as.numeric(rvonmises(10000, mu, kappa)) / (2*pi) * 360)
# Expected: maximum likelihood parameter estimate:
#   mu; stderr(mu); kappa; stderr(kappa)
#   199.4 degrees; 12.2 degrees; 1.07; 0.26
#
# Another test data:
#   degrees <- c(1.9, 12.4, 28.1, 28.9, 41.5, 46.0, 55.5, 56.6, 72.6, 75.5,
#                86.1, 109.6, 111.0, 115.3, 123.3, 127.6, 139.6, 140.8, 142.0, 147.5,
#                147.7, 149.8, 150.3, 154.1, 160.0, 161.9, 162.1, 162.4, 162.7, 163.1,
#                163.7, 168.1, 170.2, 170.4, 171.9, 172.2, 172.5, 175.4, 175.6, 175.7,
#                176.5, 177.1, 177.7, 179.0, 179.4, 179.7, 180.6, 180.7, 180.7, 181.1,
#                181.7, 182.0, 183.8, 184.1, 185.3, 188.5, 188.8, 189.0, 189.8, 192.6,
#                193.9, 194.9, 195.5, 195.7, 195.9, 196.0, 196.2, 196.4, 196.6, 198.0,
#                199.2, 199.8, 202.4, 202.9, 204.8, 206.7, 207.6, 210.5, 210.9, 212.4,
#                212.5, 213.1, 214.8, 218.0, 219.6, 220.6, 220.7, 224.5, 228.2, 232.4,
#                254.5, 255.0, 266.0, 277.4, 282.9, 289.4, 295.4, 301.1, 326.4, 354.9)
#   rads <- degrees / 360 * (2*pi)
#   mu <- atan2(sum(sin(rads)), sum(cos(rads))) %% (2*pi)
#   kappa <- est.kappa(rads)
# Expected: maximum likelihood parameter estimate:
#   mu; stderr(mu); kappa; stderr(kappa)
#   183.3 degrees; 5.9 degrees; 1.55; 0.21

