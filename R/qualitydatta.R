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

# Datta 1: Average pixel intensity
avgIntensity <- function(img.hsv) { mean(extractHSVChannel(img.hsv, VALUE)) }

# Datta 2: Colorfulness
centers <- # LUV values of bucket centers, precomputed
# Function:
colorfulness <- function(img) {
  n <- 4  # number of buckets per channel
  cuts <- cut(img, breaks = seq(0, 1, length.out=n+1), labels = 0:(n-1), include.lowest=TRUE)
  b <- array(as.integer(cuts), dim=dim(img)) - 1
  buckets <- n * (n * b[,,1] + b[,,2]) + b[,,3]
  D2 <- table(buckets)

  # Convert to LUV color space
  # http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html
  # http://127.0.0.1:28721/library/grDevices/html/convertColor.html
  # http://pages.cs.wisc.edu/~dyer/cs766/hw/hw2/code/rgb2luv.m
  # https://www.easyrgb.com/en/math.php
  # any good?
  
  # Earth mover's algorithm

  # Hungarian algorithm?

  'Work in progress!'
}

