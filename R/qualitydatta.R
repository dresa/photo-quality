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
##########
centers <- 'FILL!' # LUV values of bucket centers, precomputed

# Convert RGB image to an LUV image. RGB values are assumed to be between 0 and 1.
toLUV <- function(img.rgb) {
}
# Function:
colorfulness <- function(img) {
  n <- 4  # number of buckets per channel
  cuts <- cut(img, breaks = seq(0, 1, length.out=n+1), labels = 0:(n-1), include.lowest=TRUE)
  b <- array(as.integer(cuts), dim=dim(img)) - 1
  buckets <- n * (n * b[,,1] + b[,,2]) + b[,,3]
  D2 <- table(buckets)

  # Convert to LUV color space
  # http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html
  # http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
  # ?convertColor
  # http://pages.cs.wisc.edu/~dyer/cs766/hw/hw2/code/rgb2luv.m
  # https://www.easyrgb.com/en/math.php
  # any good?
  toXYZ(img)

  # Earth mover's distance
  # http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/RUBNER/emd.htm

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



#img <- readImage('../examples/small_grid.png')
#print(colorfulness(img))

#print(toXYZ(img))
#print(toXYZ(img, 3) - toXYZ(img, 4))
#print(toXYZ(img, 1) - toXYZ(img, 4))

#stop()

#img <- readImage('../examples/K5_10994.JPG')
#n <- 10
#for (m in 1:4) {
#  start.time <- Sys.time()
#  replicate(n, dim(toXYZ(img, m)))
#  end.time <- Sys.time()
#  time.m <- (end.time - start.time) / n
#  print(paste('Method', m, 'time:', time.m))
#}

