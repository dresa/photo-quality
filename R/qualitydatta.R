# Photo quality measures by Datta et al
#
# Implementation based on the description in the following article:
#   Studying Aesthetics in Photographic Images Using a Computational Approach
#   Datta, Ritendra; Joshi, Dhiraj; Li, Jia; Wang, James Z.
#   Computer Vision – ECCV 2006: 9th European Conference on Computer Vision, Graz, Austria, May 7-13, 2006, Proceedings, Part III
#   https://doi.org/10.1007/11744078_23
#
# Esa Junttila
# 2017-08-27 (originally)

library(emdist)  # Earth Mover's Distance

source('photo.R')
source('colorspace.R')
source('imageio.R')

# Datta 1: Average pixel intensity
##########
avgIntensity <- function(img.hsv) { mean(extractHSVChannel(img.hsv, VALUE)) }

# Datta 2: Colorfulness
# Not tested or verified yet.
# Lower EMD distances mean that colors are evenly spread, which means higher colorfulness.
# For example EMD for Black & White images may be around 82,
# a narrow set of pure colors may get EMD like 48,
# whereas an image where a wide variety colors and tones are represented has distances like.
# Not that scaling of distance is based on LUV space dimensions.
##########
colorfulness <- function(img) {
  n <- 4  # number of buckets per channel, in total n^3 buckets
  size <- prod(dim(img)[1:2])
  cuts <- cut(img, breaks = seq(0, 1, length.out=n+1), labels = 0:(n-1), include.lowest=TRUE)
  b <- array(as.integer(cuts), dim=dim(img)) - 1
  red <- 1
  green <- 2
  blue <- 3
  buckets <- n * (n * b[,,blue] + b[,,green]) + b[,,red] + 1  # bucket ID [1;64] for each _pixel_
  n.buckets <- n^3
  # Example of bucket IDs: ID's [0,n-1] have Blue and Green at lowest level, while Red varies from low to high.
  D2 <- table(buckets)
  # Bucket centers (RGB values)
  midpoints <- seq(1/(2*n), 1, 1/n)  # average RGB channel values for buckets
  domain <- expand.grid(midpoints, midpoints, midpoints)
  channels <- c('Red', 'Green', 'Blue')
  colnames(domain) <- channels
  #print(domain)
  bucket.rgb <- array(
    c(sapply(channels, function(x) unlist(domain[x]))),
    dim=c(n.buckets, 1, 3),
    dimnames=list(NULL, NULL, channels)
  )
  #print(bucket.rgb)
  bucket.luv <- toLUV(bucket.rgb)  # there is no black color (or NAs) since we use bucket midpoints
  #print(bucket.luv)
  
  # Create an EMD problem instance (Earth Mover's Distance): location and weights
  locations <- matrix(bucket.luv, ncol=3)
  #print(locations)
  #print(D2)
  from.w <- sapply(1:n^3, function(x) D2[toString(x)]) / size
  from.w[is.na(from.w)] <- 0
  names(from.w) <- 1:n.buckets
  #print(from.w)
  to.w <- replicate(n.buckets, 1 / n.buckets)
  names(to.w) <- 1:n.buckets
  #print(to.w)

  # Solve EMD instance
  #print(dim(locations))
  #print(length(from.w))
  #print(length(to.w))
  #print(locations)
  #print(from.w)
  #print(to.w)
  
  e <- emdw(locations, from.w,  locations, to.w) #, flows=TRUE)
  print(paste('EMD distance (not tested):', e))
  
  
  #rgb.centers=as.list(sapply(, function(x) ))
  #names(rgb.centers) <- sapply(0:(n.buckets-1), toString)
  
  # Convert to LUV color space
  # http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html
  # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  # ?convertColor
  # http://pages.cs.wisc.edu/~dyer/cs766/hw/hw2/code/rgb2luv.m
  # https://www.easyrgb.com/en/math.php
  # http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Luv.html
  # any good?

  
  # Earth mover's distance
  # http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/RUBNER/emd.htm
  #
  # Example:
  # Weight:   20 10  0  -->  10  2 18   
  # Location:  1  2  3        1  2  3   
  # Optimal moves: move 10 from 1 to 2, and move 18 from 2 to 3. Or move 10 from 1 to 3, and move 8 from 2 to 3.
  # The cost of moves is 28. We have 30 units in total, so distance is 28/30=0.9333333.
  #
  # > emdw(matrix(c(1,2,3)), matrix(c(20,10,0)),  matrix(c(1,2,3)), matrix(c(10,2,18)))
  # [1] 0.9333333
  # Specifically, EMD assigns the following flows (from 'flows' argument):
  # 0->0, 10 units; 0->1, 2 units; 0->2, 8 units; 1->2, 10 units; in total 30, so 28/30=0.9333333
  #
  # Replicated 2D example:
  # > A <- matrix(1:6 / sum(1:6), 2)
  # > B <- matrix(c(0, 0, 0, 0, 0, 1), 2)
  # > emd2d(A, B)
  # [1] 0.9275576
  # > A
  # [,1]      [,2]      [,3]
  # [1,] 0.04761905 0.1428571 0.2380952
  # [2,] 0.09523810 0.1904762 0.2857143
  # > B
  # [,1] [,2] [,3]
  # [1,]    0    0    0
  # [2,]    0    0    1
  # Replicated as Euclidean distance on the grid:
  # > f <- function(i,j) A[i,j]*sqrt((2-i)^2 + (3-j)^2)
  # > sum(f(1,1), f(1,2), f(1,3), f(2,1), f(2,2), f(2,3))
  # [1] 0.9275576
  #
  # Third example:
  # Weight:      w 3 20  5 19  3  -->  5 12 12 16  5
  # 2D Location: x 2  3  5  8  9       2  3  5  8  9   
  #              y 1  2  4  7  9       1  2  4  7  9
  # Locations:
  #   L <- t(matrix(c(2,1, 3,2, 5,4, 8,7, 9,9), 2))  # locations by rows
  # Signatures (weights):
  #   from.W <- c(3, 20, 5, 19, 3)
  #   to.W <- c(5, 12, 12, 16, 5)
  # > emdw(L, from.W,  L, to.W, flows=TRUE)
  # [1] 0.5702753
  # attr(,"flows")
  # attr(,"flows")[[1]]
  # [1] 0 4 1 1 3 3 2 1 3  # from location
  # attr(,"flows")[[2]]
  # [1] 0 4 1 0 3 4 2 2 2  # to location
  # attr(,"flows")[[3]]
  # [1]  3  3 12  2 16  2  5  6  1  # how many units to move?

  
  # Hungarian algorithm?

  'Work in progress!'
}



img <- readImage('../examples/small_grid.png')     # 48? (EMD colorfulness distance), not verified
#img <- readImage('../examples/sharp_or_blur.png')  # 83?
#img <- readImage('../examples/penguin.jpg')        # 61?
#img <- readImage('../examples/no_shift.png')       # 57?
#img <- readImage('../examples/many_colors.png')    # 18?
#img <- readImage('../examples/niemi.png')          # 49?
#img <- readImage('../examples/almost_black.png')   # 83?
#img <- readImage('../examples/dark_city.png')      # 65?

#source('viewer.R')
#view(img)
#print(toXYZ(img))

print(colorfulness(img))


#img <- readImage('../examples/K5_10994.JPG')
#n <- 10
#start.time <- Sys.time()
#replicate(n, dim(toXYZ(img)))
#end.time <- Sys.time()
#print((end.time - start.time) / n)

