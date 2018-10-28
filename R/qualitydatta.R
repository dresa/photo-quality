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

library(emdist)    # Earth Mover's Distance
library(waveslim)  # Discrete Wavelet Transform 2D
library(mmand)     # Segmented values
#library(NbClust)   # Optimal number of clusters

source('R/photo.R')
source('R/colorspace.R')
source('R/imageio.R')
source('R/common.R')  # k-means++
source('R/viewer.R')
source('R/mdl.R')
source('R/convexity.R')



# Datta 2: Colorfulness
# Lower EMD distances mean that colors are evenly spread, which means higher
# colorfulness. For example EMD for Black & White images may be around 82,
# a narrow set of pure colors may get EMD like 48, whereas images with a wide
# variety of colors have distances like 20. Note that the distance and its
# range is based on the CIE LUV color space.
##########

# Map the colors in an RGB image into a small number (n^3) of buckets.
allocateColorsToBuckets <- function(img, n) {
  DEBUG_BUCKETS <- FALSE
  # Breaks that separate RGB color-value intervals 1..n
  breaks <- seq(0, 1, length.out=n+1)
  # Map all R, G, and B channels values [0;1] into inverval IDs {0,1,..,n}
  b <- array(findInterval(img, breaks, rightmost.closed=TRUE) - 1, dim=dim(img))
  # Combine R, G, and B interval codes to form pixel-wise bucket codes, ranged {1,2,...,n^3}
  red <- 1
  green <- 2
  blue <- 3
  genBucketCodes <- function() n * (n * b[,,blue] + b[,,green]) + b[,,red] + 1
  # Example: bucket ID's [0;n-1] have low Blue and Green, while Red varies from low to high.
  # Frequency of bucket-colors within the RGB image:
  bucket.freqs <- tabulate(genBucketCodes(), nbins=n^3)
  
  # Debugging: show distribution of color-buckets:
  if (DEBUG_BUCKETS) {
    signs <- function(ch, x) paste(rep(ch, x), collapse='')
    minus.signs <- sapply(floor(n/2):1, function(x) signs('-', x))
    neutral.signs <- rep('', if (n%%2) 1 else 0)
    plus.signs <- sapply(1:floor(n/2), function(x) signs('+', x))
    dn <- function(s) paste0(s, c(minus.signs, neutral.signs, plus.signs))
    print(array(bucket.freqs, dim=c(n,n,n), dimnames=list(dn('r'), dn('g'), dn('b'))))
    print(array(1:n^3, dim=c(n,n,n)))
    p <- bucket.freqs/sum(bucket.freqs)
    print(paste('Bucket entropy:', -sum(p*log2(p), na.rm=TRUE)))
  }
  
  return(bucket.freqs)  # variable D2 in original article
  # Note that 'findInterval' is 8x faster than 'cut' in this case. Also, 'tabulate' is 8x faster
  # than 'table' function (because of implicit int->str).
}

# Generate the RGB colors of bucket centerpoints (fixed for every image).
bucketColorsRGB <- function(n) {
  # Bucket center colors: central RGB values converted into LUV values
  # There is no black color (or NAs) since we use bucket midpoints, so no problem with LUV.
  n.buckets <- n^3
  midpoints <- seq(1/(2*n), 1, 1/n)  # average RGB channel values for buckets
  domain <- expand.grid(midpoints, midpoints, midpoints)
  channels <- c('Red', 'Green', 'Blue')
  colnames(domain) <- channels
  bucket.rgb <- array(
    c(sapply(channels, function(x) unlist(domain[x]))),
    dim=c(n.buckets, 1, 3),
    dimnames=list(NULL, NULL, channels)
  )
  #print(bucket.rgb)  # debug
  return(bucket.rgb)
}

BUCKET_LUV_CACHE <- list()

# Generate the LUV colors of bucket centerpoints (fixed for every image).
bucketColorsLUV <- function(n) {
  if (is.null(unlist(BUCKET_LUV_CACHE[n]))) {
    BUCKET_LUV_CACHE[[n]] <<- toLUV(bucketColorsRGB(n))
  }
  return(BUCKET_LUV_CACHE[[n]])
}

# Create an instance of a Earth Mover's Distance problem (EMD),
# where we model color coordinates in LUV space as locations,
# the frequency of bucket colors with the image as source distribution,
# and uniform distribution of bucket colors as target distribution (signature).
createEMDProblem <- function(bucket.luv, bucket.freqs, n, target='all') {
  n.buckets <- n^3
  locations <- matrix(bucket.luv, ncol=3, dimnames=list(NULL, dimnames(bucket.luv)[[3]]))
  # Weights, distributions, and signatures all mean the same thing.
  from.w <- bucket.freqs / sum(bucket.freqs)  # SOURCE distribution
  from.w[is.na(from.w)] <- 0    # missing buckets have zero frequency
  names(from.w) <- 1:n.buckets  # just for debugging
  to.w <- switch(tolower(target),
                 uniform={setNames(replicate(n.buckets, 1/n.buckets), 1:n.buckets)},
                 grey={
                   g <- integer(n.buckets)  # zeros
                   # use grey buckets; with n=4, buckets are 1, 22, 43, 64, from black to white
                   GREY_BUCKETS <- seq(1, n.buckets, (n.buckets-1)/(n-1))
                   g[GREY_BUCKETS] <- 1/length(GREY_BUCKETS)
                   setNames(g, 1:n.buckets)}
  )
  #if (lower(target) == 'all') to.w <- replicate(n.buckets, 1 / n.buckets)  # evenly distributed target
  #else if (lower(target) == 'grey') to.w <- replicate(n.buckets, 1 / n.buckets)  # evenly distributed target
  return(list(locations=locations, from.weights=from.w, to.weights=to.w))
}

# Debugging: using the EMD solution flows to replicate the result
# based on Euclidean distance. Slow, but independent.
replicateDistance <- function(emd.solution, locations, n) {
  n.buckets <- n^3
  flows <- attributes(emd.solution)[['flows']]
  froms <- flows[[1]] + 1  # from 0..(n^3-1) domain to 1:n^3
  tos <- flows[[2]] + 1
  amounts <- flows[[3]]
  euclidean <- function(x) sqrt(sum((locations[x[1], ] - locations[x[2], ])^2))
  distances <- matrix(apply(expand.grid(1:n.buckets, 1:n.buckets), 1, euclidean), ncol=n.buckets)
  d <- setNames(apply(distances, 1, sum)/n.buckets, 1:n.buckets)
  print(paste('Farthest bucket', which.max(d), 'has uniform distance', max(d)))
  print(paste('Closest bucket', which.min(d), 'has uniform distance', min(d)))
  f <- function(from, to, amount) distances[from, to] * amount
  costs <- mapply(f, froms, tos, amounts)
  print(paste('Flow amount', amounts, 'from', froms, 'to', tos, 'with cost', costs), collapse='\n')
  total.cost <- sum(costs)
  normalized.cost <- total.cost / sum(amounts)
  return(normalized.cost)
}

# Datta 2 measure: colorfulness for an RGB image, as measured in LUV space
# n is the number of buckets per channel, in total n^3 buckets
###########
#' Colorfulness of an image (Datta measure 2)
#' 
#' Compute \emph{colorfulness} of an image: the distance to a reference
#' image that has uniformly distributed colors in LUV color space.
#' 
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' 
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3)
#' @param mode name of the reference image composition, distribution of colors
#' @param n number of color buckets; uses \emph{n^3} buckets to
#'   count frequencies
#' @return numeric \emph{colorfulness} measure; the maximum value depends
#'   on the number of buckets \emph{n^3}.
#' @examples
#' set.seed(1)
#' img.rgb <- array(runif(4*5*3, 0.4, 1.0), dim=c(4,5,3))
#' abs(colorfulness(img.rgb) - 41.06221) < 1e-5
#' set.seed(1)
#' img.unif <- array(runif(4*5*3), dim=c(4,5,3))
#' abs(colorfulness(img.unif) - 29.65533) < 1e-5
#' set.seed(1)
#' img.purple <- array(c(rep(1, 4*5), runif(4*5), runif(4*5, 0.5, 1)), dim=c(4,5,3))
#' abs(colorfulness(img.purple) - 56.10874) < 1e-5
#' 
#' @export
colorfulness <- function(img, mode='uniform', n=4) {
  # Preliminary constants
  EMD_DEBUG <- FALSE
  
  # Label image pixels with a number of color-buckets and create an EMD problem instance:
  bucket.freqs <- allocateColorsToBuckets(img, n)
  bucket.luv <- bucketColorsLUV(n)
  p <- createEMDProblem(bucket.luv, bucket.freqs, n, mode)
  # Solve EMD instance. We can choose 'from' and 'to' either way (Euclidean is symmetric).
  e <- emdw(p$locations, p$from.weights, p$locations, p$to.weights, dist='euclidean', flows=TRUE)
  
  # Debugging: given the flows, replicate the EMD distance?
  if (EMD_DEBUG) {
    print(bucket.luv)
    print(p)
    rcost <- replicateDistance(e, p$locations, n)
    print(paste('DEBUG: Normalized EMD cost (replicated):', rcost))
  }
  
  return(as.numeric(e))  # drop flow information
}


#' @rdname avgHsvPixels
#' @name avgHsvPixels
#' @aliases avgHue
#' @aliases avgSaturation
#' @aliases avgIntensity
#' @aliases avgCentralHue
#' @aliases avgCentralSaturation
#' @aliases avgCentralIntensity
#' @title
#' Pixel HSV Averages (Datta measures 1 & 3--7)
#' @description
#' Compute average HSV (hue, saturation, value) statistics for pixels
#' within a full image or in a central area of an image, as defined
#' by \emph{Datta measures 1 & 3--7}.
#' @details
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' 
#' In and \code{avgHue} and \code{avgCentralHue}, the pixels that do not
#' have a hue, such as grayscale pixels, are excluded. The functions use the
#' \emph{arithmetic} definition of "mean" rather than a \emph{circular} one.
#' @param img.hsv photo as an HSV image array (\emph{m} x \emph{n} x 3)
#' @examples
#' set.seed(1)
#' img.rgb <- array(runif(4*5*3, 0.4, 1.0), dim=c(4,5,3))
#' img.hsv <- toHSV(img.rgb)
#' 
NULL


# Datta 4: Average pixel hue (circular)
##########
#' @description
#' \code{avgHue} computes the average pixel \strong{hue}
#' in an image, as in \emph{Datta measure 4}.
#' @return
#'   \code{avgHue} returns the mean hue of pixels, within \code{[0;360]}
#'   degrees or \code{[0;2*pi]} radians.
#' @examples
#' abs(avgHue(img.hsv) - 201.3679604) < 1e-6
#' abs(avgHue(toHSV(img.rgb, radians=TRUE)) - 201.3679604/360*2*pi) < 1e-6
#' 
#' @rdname avgHsvPixels
#' @export
avgHue <- function(img.hsv) { mean(extractHSVChannel(img.hsv, HUE), na.rm=TRUE) }


# Datta 3: Average pixel saturation
##########
#' @description
#' \code{avgSaturation} computes the average pixel \strong{saturation}
#' in an image, as in \emph{Datta measure 3}.
#' @return \code{avgSaturation} returns the mean saturation of pixels, ranging from 0 to 1.
#' @examples
#' abs(avgSaturation(img.hsv) - 0.3379017) < 1e-6
#' 
#' @rdname avgHsvPixels
#' @export
avgSaturation <- function(img.hsv) { mean(extractHSVChannel(img.hsv, SATURATION)) }


# Datta 1: Average pixel intensity
##########
#' @description
#' \code{avgIntensity} returns the average pixel \strong{intensity}
#' in an image, as in \emph{Datta measure 1}.
#' @return
#'   \code{avgIntensity} returns the mean intensity (\emph{value}) of pixels, ranging from 0 to 1.
#' @examples
#' abs(avgIntensity(img.hsv) - 0.852) < 1e-3
#' 
#' @rdname avgHsvPixels
#' @export
avgIntensity <- function(img.hsv) { mean(extractHSVChannel(img.hsv, VALUE)) }


avgCentral <- function(img.hsv.channel) {
  ch <- img.hsv.channel
  mid.cols <- seq(floor(ncol(ch)/3), ceiling(2*ncol(ch)/3))
  mid.rows <- seq(floor(nrow(ch)/3), ceiling(2*nrow(ch)/3))
  return(mean(ch[mid.rows, mid.cols], na.rm=TRUE))
}

# Datta 5: Average pixel hue (circular) in the middle part of the image [one-third,two-thirds]
##########
#' @description
#' \code{avgCentralHue} computes the average pixel \strong{hue}
#' in the \strong{central} part an image, as in \emph{Datta measure 5}.
#' @return
#'   \code{avgCentralHue} returns the mean hue of pixels in central part,
#'   within \code{[0;360]} degrees or \code{[0;2*pi]} radians.
#' @examples
#' abs(avgCentralHue(img.hsv) - 208.99541170) < 1e-6
#' abs(avgCentralHue(toHSV(img.rgb, radians=TRUE)) - 208.99541170/360*2*pi) < 1e-6
#' 
#' @rdname avgHsvPixels
#' @export
avgCentralHue <- function(img.hsv) { return(avgCentral(extractHSVChannel(img.hsv, HUE))) }


# Datta 6: Average pixel saturation in the middle part of the image
##########
#' @description
#' \code{avgCentralSaturation} computes the average pixel \strong{saturation}
#' in the \strong{central} part an image, as in \emph{Datta measure 6}.
#' @return
#'   \code{avgCentralSaturation} returns the mean saturation of pixels
#'   in central part, ranging from 0 to 1.
#' @examples
#' abs(avgCentralSaturation(img.hsv) - 0.3177584) < 1e-6
#' 
#' @rdname avgHsvPixels
#' @export
avgCentralSaturation <- function(img.hsv) { return(avgCentral(extractHSVChannel(img.hsv, SATURATION))) }


# Datta 7: Average pixel intensity in the middle part of the image
##########
#' @description
#' \code{avgCentralIntensity} computes the average pixel \strong{intensity}
#'  (value) in the \strong{central} part an image, as in \emph{Datta measure 7}.
#' @return
#'   \code{avgCentralIntensity} returns the mean intensity (value) of pixels
#'   in central part, ranging from 0 to 1.
#' @examples
#' abs(avgCentralIntensity(img.hsv) - 0.8340243) < 1e-6
#' 
#' @rdname avgHsvPixels
#' @export
avgCentralIntensity <- function(img.hsv) { return(avgCentral(extractHSVChannel(img.hsv, VALUE))) }


# Datta 8 (top-20 familiarity) skipped due to lack of data!
##########

# Datta 9 (top-100 familiarity) skipped due to lack of data!
##########


# Datta 10--21 Texture (graininess)
##########

# Note: function "waveletTexure" used in Datta 53--55 as well.
waveletTexture <- function(img.channel) {
  ASSUMED_LEVELS <- 3
  levels <- min(c(ASSUMED_LEVELS, floor(log2(dim(img.channel)))))
  wavelet.filter <- 'd4'  # Daubechies wavelet 'd4'
  # Perform 2D Discrete Wavelet Transform
  d <- dwt.2d(img.channel, wf=wavelet.filter, J=levels)  # from 'waveslim' package
  # reconstructed <- idwt.2d(d)  # up to image dimensions that are power of two
  return(d)
}

waveletDattaFrequencies <- function(img.channel) {
  d <- waveletTexture(img.channel)
  num.levels <- attr(d, 'J')
  # Compute Datta measures on different frequencies (levels):
  freq.1 <- if (num.levels>=1) sum(abs(d$HH1) + abs(d$HL1) + abs(d$LH1), na.rm=TRUE) / (3 * length(d$HH1)) else NA
  freq.2 <- if (num.levels>=2) sum(abs(d$HH2) + abs(d$HL2) + abs(d$LH2), na.rm=TRUE) / (3 * length(d$HH2)) else NA
  freq.3 <- if (num.levels>=3) sum(abs(d$HH3) + abs(d$HL3) + abs(d$LH3), na.rm=TRUE) / (3 * length(d$HH3)) else NA
  return(c(freq.1, freq.2, freq.3))
}  # TODO: should check if "sum of abs" is the correct way of doing this measure.


#' Texture of an image (Datta measures 10--21)
#' 
#' Compute the \emph{texture} features of an image, as in Datta measures 10--21.
#' 
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' 
#' @param img.hsv photo as an HSV image array (\emph{m} x \emph{n} x 3)
#' @return list of twelve numeric \emph{texture} features:
#'   hues (at three levels 1, 2, and 3), saturation (at three levels),
#'   value (at three levels), sum of hues, sum of saturations,
#'   and sum of values. Hues have been rescaled from degrees onto \code{[0;1]}.
#'   By default three texture levels are computed, but if the image has
#'   fewer than \emph{2^3=8} rows or columns, some texture levels will
#'   be \code{NA}.
#' @examples
#' set.seed(1)
#' img.rgb <- array(runif(8*10*3, 0.4, 1.0), dim=c(8,10,3))
#' expected <- c(0.2303976, 0.255253, 0.1964088, 0.1054641, 0.1160977, 0.1132417, 0.08682113, 0.08736958, 0.1347886, 0.6820594, 0.3348034, 0.3089793)
#' all(abs(unlist(texture(toHSV(img.rgb))) - expected) < 1e-5)
#' 
#' @export
texture <- function(img.hsv) {
  img.hsv[is.na(img.hsv)] <- 0  # local changes; DWT works only non-NA values
  res <- lapply(HSV, function(channel) waveletDattaFrequencies(extractHSVChannel(img.hsv, channel)))
  return(list(
    hue.1=res[[HUE$idx]][1] / 360,  # rescale HUE from 0--360 degrees to 0--1
    hue.2=res[[HUE$idx]][2] / 360,
    hue.3=res[[HUE$idx]][3] / 360,
    sat.1=res[[SATURATION$idx]][1],
    sat.2=res[[SATURATION$idx]][2],
    sat.3=res[[SATURATION$idx]][3],
    val.1=res[[VALUE$idx]][1],
    val.2=res[[VALUE$idx]][2],
    val.3=res[[VALUE$idx]][3],
    hue.sum=sum(res[[HUE$idx]], na.rm=TRUE) / 360,  # rescale HUE from 0--360 degrees to 0--1
    sat.sum=sum(res[[SATURATION$idx]], na.rm=TRUE),
    val.sum=sum(res[[VALUE$idx]], na.rm=TRUE)
  ))  # should we use absolute value??
}


#' @rdname dimFeatures
#' @name dimFeatures
#' @aliases sizeFeature
#' @aliases aspectRatioFeature
#' @title
#' Dimensional Photo Features (Datta measures 22 & 23)
#' @description
#' Identify dimensional features of an image, as defined
#' by \emph{Datta measures 22 & 23}.
#' @details
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' @param img photo as an image array (\emph{m} x \emph{n} x ?)
#' @examples
#' img <- array(runif(4*5*3, 0.4, 1.0), dim=c(4,5,3))
#' 
NULL

# Datta 22: size
###########
#' @description
#'   \code{sizeFeature} measures image size: rows \strong{plus} columns,
#'   as in \emph{Datta measure 22}.
#' @return
#'   \code{sizeFeature} returns the image size that is \emph{m} \strong{+} \emph{n}.
#' @examples
#' sizeFeature(img) == 9
#' 
#' @rdname dimFeatures
#' @export
sizeFeature <- function(img) {
  d <- dim(img)
  return(d[1] + d[2])
}

# Datta 23: aspect ratio
###########
#' @description
#'   \code{aspectRatioFeature} measures the aspect ratio of the image,
#'   as in \emph{Datta measure 23}.
#' @return
#'   \code{aspectRatioFeature} returns the aspect ratio of image dimensions.
#'   It is a vector of two elements:
#'   \itemize{
#'     \item \emph{m/n}, which is the number of rows divided by that of columns
#'     \item \emph{m/n} or its reciprocal \emph{n/m}, whichever is larger
#'       (ratio is at least 1)
#'   }.
#' @examples
#' all(aspectRatioFeature(img) == c(4/5, 5/4))
#' 
#' @rdname dimFeatures
#' @export
aspectRatioFeature <- function(img) {
  d <- dim(img)
  return(c(d[1]/d[2], max(d[1]/d[2], d[2]/d[1])))
}


# Datta 24--52: Region Composition
###########

photoSegmentationSample <- function(img.luv) {
  # Create a dataset, might downsample if needed
  d <- dim(img.luv)
  nr <- d[1]
  nc <- d[2]
  pixels <- nr * nc
  SIZE_LIMIT <- 10000  # in pixels, for performance reasons
  needs.thinning <- SIZE_LIMIT < pixels
  if (needs.thinning) {  # downsampling needed
    thin.factor <- sqrt(pixels / SIZE_LIMIT)
    ds.rows <- round(seq(1, nr, length.out=nr/thin.factor))
    ds.cols <- round(seq(1, nc, length.out=nc/thin.factor))
    x <- matrix(img.luv[ds.rows, ds.cols, ], ncol=3)
    d.new <- c(length(ds.rows), length(ds.cols))
  } else {
    x <- matrix(img.luv, ncol=3)
    d.new <- c(nr, nc)
  }
  colnames(x) <- c('L', 'U', 'V')
  return(list(points=x, dim=d.new))
}

photoSegmentationClustering <- function(x, k, iter.max, restarts) {
  # Clustering in LUV space
  # If k=0, then simply return the column means;
  # not a clustering but to bring robustness in viewing images.
  if (k>=1) { clustering <- kmeanspp(x, k, iter.max=iter.max, restarts=restarts) }
  else { clustering <- list(centers=matrix(colMeans(x), ncol=3), k=0) }
  return(clustering)
}

photoSegmentationMappings <- function(pixels, centers, nr) {
  # For each center, find closest pixels iteratively.
  # Number of centers should be small, < 100, so for loop should be fast enough.
  pxs <- t(pixels)
  min.distance <- rep(Inf, nrow(pixels))  # min dist to a center this far
  min.mapping <- rep(NA, nrow(pixels))  # closest center this far
  for (c.idx in 1:nrow(centers)) {  # check one center at a time
    d <- colSums((pxs - centers[c.idx, ])^2)  # distances to this center
    upd <- d < min.distance  # new closest center found for pixels?
    min.distance[upd] <- d[upd]
    min.mapping[upd] <- c.idx
  }
  # return center mappings in original matrix dimensions
  return(matrix(unlist(min.mapping), nrow=nr))
}

# Slower (by x20), but elegant version. A call for each pixel is too much!
#photoSegmentationMappings <- function(pixels, centers, nr) {
#  C <- t(centers)
#  chooseCluster <- function(p) which.min(colSums((C - p)^2))
#  mapping <- apply(pixels, 1, chooseCluster)
#  return(matrix(unlist(mapping), nrow=nr))
#}

photoSegmentationReconstruct <- function(centers, map, num.rows) {
  return(t(sapply(map, function(x) centers[x, ])))
}

#photoSegmentationComponents <- function(clusters.2d) {
#  Recursive version, too deep
#  segDFS <- function(i, j) {  # recursive Depth-First Search for identifying connected segments
#    segments.2d[i, j] <<- seg.id
#    if (i<nr && segments.2d[i+1,j  ] == NO_SEG && clusters.2d[i+1,j  ] == cluster.id) segDFS(i+1, j  )
#    if (j<nc && segments.2d[i  ,j+1] == NO_SEG && clusters.2d[i  ,j+1] == cluster.id) segDFS(i  , j+1)
#    if (i>1  && segments.2d[i-1,j  ] == NO_SEG && clusters.2d[i-1,j  ] == cluster.id) segDFS(i-1, j  )
#    if (j>1  && segments.2d[i  ,j-1] == NO_SEG && clusters.2d[i  ,j-1] == cluster.id) segDFS(i  , j-1)
#  }
#  NO_SEG <- 0
#  nr <- nrow(clusters.2d)
#  nc <- ncol(clusters.2d)
#  segments.2d <- matrix(NO_SEG, nrow=nr, ncol=nc)
#  seg.id <- 0
#  if (nr >= 1 & nc >= 1) {
#    # quite inefficient, I know...should replace by something?? TODO
#    for (i in 1:nr) {
#      for (j in 1:nc) {
#        if (segments.2d[i, j] == NO_SEG) {
#          seg.id <- seg.id + 1  # new id for another connected component
#          cluster.id <- clusters.2d[i, j]
#          segDFS(i, j)
#        }
#      }
#    }
#  }
#  return(segments.2d)
#}


#photoSegmentationComponents <- function(clusters.2d) {
#  # Very slow whileloop + stack segmentation. Should find a reasonable method instead.
#  NO_SEG <- 0
#  nr <- nrow(clusters.2d)
#  nc <- ncol(clusters.2d)
#  segments.2d <- matrix(NO_SEG, nrow=nr, ncol=nc)
#  seg.id <- 0
#  if (nr >= 1 & nc >= 1) {
#    # Initialize a stack: all pixels
#    stack <- mapply(c, rev(rep(1:nr, nc)), rev(rep(1:nc, each=nr)), SIMPLIFY=FALSE)
#    # Introduce stack methods: pop and push
#    st.pop <- function() {
#      n <- length(stack)
#      rc <- stack[[n]]
#      stack <- stack[-n]
#      stack <<- stack
#      return(rc)
#    }
#    st.push <- function(rc) {
#      n <- length(stack)
#      stack[n + 1] <- list(rc)
#      stack <<- stack
#    }
#    # Perform a DFS (depth-first search) with a while loop and stack
#    pixels.left <- length(stack)
#    while (length(stack)) {
#      rc <- st.pop()
#      i <- rc[1]
#      j <- rc[2]
#      if (segments.2d[i, j] == NO_SEG) {
#        if (length(stack) < pixels.left) {
#          pixels.left <- length(stack)
#          seg.id <- seg.id + 1
#          cluster.id <- clusters.2d[i, j]
#          print(length(stack))
#        }
#        segments.2d[i, j] <- seg.id
#        if (i<nr && segments.2d[i+1,j  ] == NO_SEG && clusters.2d[i+1,j  ] == cluster.id) st.push(c(i+1, j  ))
#        if (j<nc && segments.2d[i  ,j+1] == NO_SEG && clusters.2d[i  ,j+1] == cluster.id) st.push(c(i  , j+1))
#        if (i>1  && segments.2d[i-1,j  ] == NO_SEG && clusters.2d[i-1,j  ] == cluster.id) st.push(c(i-1, j  ))
#        if (j>1  && segments.2d[i  ,j-1] == NO_SEG && clusters.2d[i  ,j-1] == cluster.id) st.push(c(i  , j-1))
#      }
#    }
#  }
#  return(segments.2d)
#}



# Third try, using "mmand" package in a clumsy way, separately for each cluster.

photoSegmentationComponents <- function(clusters.2d) {
  clusters <- sort(unique(c(clusters.2d)))
  nr <- nrow(clusters.2d)
  nc <- ncol(clusters.2d)
  # Separate bottom and top with an extra row of NAs ([nr,j] and [1,j+1] are neighbors)
  component.base <- rbind(clusters.2d, NA)
  segments <- matrix(nrow=nr+1, ncol=nc)  # initialized with NA
  neighbor.kernel <- matrix(c(0,1,0, 1,1,1, 0,1,0), byrow=TRUE, ncol=3)
  seg.cumu.num <- 0
  # Accumulate segments by processing each cluster (each may give rise to multiple segments)
  for (c.id in clusters) {
    x <- component.base
    x[x != c.id] <- NA
    # Compute connected components, separated by NA, with the 'mmand' package.
    s <- components(x, neighbor.kernel)
    new.segs <- s[!is.na(s)]
    segments[!is.na(s)] <- new.segs + seg.cumu.num
    seg.cumu.num <- seg.cumu.num + max(new.segs)
  }
  return(segments[1:nr, ])  # drop the extra row of NAs
}


photoSegmentationShow <- function(img.luv, centers, clustered.img, conn.components) {
  nr <- dim(img.luv)[1]
  nc <- dim(img.luv)[2]
  img.rgb <- array(convertColor(matrix(img.luv, ncol=3), from='Luv', 'sRGB'), dim=dim(img.luv))
  view(img.rgb, title='Original image')
  reconstr <- photoSegmentationReconstruct(centers, clustered.img, nr)
  reconstr.rgb <- array(convertColor(reconstr, from='Luv', 'sRGB'), dim=dim(img.luv))
  nclus <- nrow(centers)
  view(reconstr.rgb, title=paste0('Reconstructed CIE LUV (', nclus, ' nearest clusters)'))
  m <- max(conn.components)
  x <- conn.components
  rnd.title <- paste('Random-colored', m, 'connected segments from', nclus, 'clusters')
  view(array(c(x/m, ((3*x)%%m)/m, ((7*x)%%m)/m), dim=dim(img.luv)), title=rnd.title)
  seg.sizes <- rev(sort(table(conn.components)))
  NUM_SEG <- 10
  largest.ids <- as.integer(names(head(seg.sizes, NUM_SEG)))
  largest.comps <- conn.components
  largest.comps[!(largest.comps %in% largest.ids)] <- 0  # non-segments are black
  lc.val <- largest.comps
  u <- sort(unique(c(largest.comps)))
  for (i in 1:length(u)) {
    lc.val[lc.val == u[i]] <- (i-1)/(length(u)-1)
  }
  view(matrix(lc.val, nrow=nr, ncol=nc), title=paste0('Largest connected components (', NUM_SEG, ')'))
}

# Choose optimal number of clusters based on MDL (Minimum Description Length) principle.
# This function also defines the segments in the image.
chooseNumberClusters <- function(img.luv.points, num.rows, num.cols, k.min=0, k.max=30, iter.max=20, restarts=3) {
  s <- img.luv.points  # image with points as rows and LUV values as columns
  nr <- num.rows
  nc <- num.cols
  stopifnot(nr * nc == nrow(s))

  ## Using several methods to choose k is too slow. Infeasible even for moderate-size images.
  #num.clust <- NbClust(data=s, distance='euclidean', min.nc=2, max.nc=25, method='kmeans', index='alllong')
  #print(num.clust)

  computeClusteringMDL <- function(k) {
    # Sample L, U, V coordinates from uniform distributions, measured at integer granularity.
    # Computed L component values are in the range [0 to 100].
    # Computed U component values are in the range [-124 to 220].
    # Computed V component values are in the range [-140 to 116].
    # Encode model, with k clusters:
    mdl.dim <- eliasDelta(nr) + eliasDelta(nc)
    mdl.k <- log2(nr*nc)
    # Cluster centers and stdev from uniform distribution, with unit precision
    mdl.centers <- k * (log2(101) + log2(345) + log2(257))
    mdl.std <- k * log2(100)
    delta <- 0.5
    #print(paste('k=', k))
    if (k==0) {
      mdl.map <- 0
      mdl.diffs <- nr * nc * (log2(101) + log2(345) + log2(257))
    } else if (k>=1) {
      res.k <- photoSegmentationClustering(s, k, iter.max, restarts)
      s.diff <- s - res.k$centers[res.k$cluster, ]
      cluster.mat <- matrix(res.k$cluster, nrow=nr)
      mdl.map.old <- nr * nc * log2(k)  # for each pixel, record its cluster
      mdl.map.details <- mdlRle2d(cluster.mat, k)
      mdl.map <- mdl.map.details$values + mdl.map.details$widths
      mdl.diffs <- 0
      for (i in 1:k) {
        diff.part <- s.diff[res.k$cluster == i, ]
        # Encode data:
        part.sd <- max(1, round(sd(diff.part)))
        xs <- round(diff.part)
        part.probs <- (dnorm(xs-delta, sd=part.sd) + dnorm(xs, sd=part.sd) + dnorm(xs+delta, sd=part.sd)) / 3
        mdl.diffs <- mdl.diffs + sum(log2(1/(part.probs)))
      }
    } else { stop(paste('Illegal number of clusters:', k)) }
    #print(paste(c('MDLs', mdl.k, mdl.centers, mdl.std, mdl.map, mdl.diffs), collapse=' ', sep=', '))
    mdl <- mdl.k + mdl.centers + mdl.std + mdl.map + mdl.diffs
    return(mdl)
  }

  mdl <- sapply(k.min:k.max, computeClusteringMDL)
  names(mdl) <- k.min:k.max
  mdl.k <- as.integer(names(mdl)[which.min(mdl)])
  clustering <- photoSegmentationClustering(s, mdl.k, iter.max, restarts)
  return(list(clustering=clustering, mdl.values=mdl, optimal.k=clustering$k, initial.k=mdl.k))
}


# Segmentation tool for images, automatically chooses k (optimal number of clusters)
photoSegmentation <- function(img.rgb, k.min=0) {
  # Preprocess:
  # * convert to LUV color space (for meaningful clustering)
  # * and downsample resolution (for performance)
  img.luv <- toLUV(img.rgb)
  img.luv[is.na(img.luv)] <- 0  # for black 'u' and 'v' are missing (limit to inf is zero)
  sample <- photoSegmentationSample(img.luv)

  k.clusters <- chooseNumberClusters(sample$points, sample$dim[1], sample$dim[2], k.min=k.min)

  # Now we have chosen 'k' and k cluster centers
  # Now use the full image (instead of a downsampled image)
  nr.img <- nrow(img.rgb)
  nc.img <- ncol(img.rgb)
  k.centers <- k.clusters$clustering$centers
  clustered.img <- photoSegmentationMappings(matrix(img.luv, ncol=3), k.centers, nr.img)
  conn.components <- photoSegmentationComponents(clustered.img)
  #print(conn.components)

  # When debugging, enable the folloing line. Debug by viewing images:
  #photoSegmentationShow(img.luv, k.centers, clustered.img, conn.components)  # DEBUG: show intermediate images

  return(
    list(
      img.rgb=img.rgb,
      img.luv=img.luv,
      clustering = k.clusters$clustering,
      mdl.values=k.clusters$mdl.values,
      optimal.k=k.clusters$optimal.k,
      conn.components=conn.components
    )
  )
}


# Return the bounding rectangular for each id. That is, the minimum
# rectangular that contains all given id values within matrix m.
# Assuming the IDs are 1, 2, 3, ..., max.id.
getBoundingBoxes <- function(m, max.id) {
  first.row <- integer(length=max.id) * NA
  last.row <- integer(length=max.id) * NA
  first.col <- integer(length=max.id) * NA
  last.col <- integer(length=max.id) * NA
  for (r in 1:nrow(m)) {
    found.ids <- unique(m[r, ])
    upd.mask <- is.na(first.row[found.ids])
    first.row[found.ids[upd.mask]] <- r
    last.row[found.ids] <- r
  }
  for (col in 1:ncol(m)) {
    found.ids <- unique(m[ , col])
    upd.mask <- is.na(first.col[found.ids])
    first.col[found.ids[upd.mask]] <- col
    last.col[found.ids] <- col
  }
  df <- data.frame(
    FirstRow=first.row,
    LastRow=last.row,
    FirstColumn=first.col,
    LastColumn=last.col
  )
  return(df)
}


# Visualize convex shapes found in the photo
plotConvexShapes <- function(dims, convex.shapes) {
  nr <- dims[1]
  nc <- dims[2]
  default.par <- par(mar=c(1,1,2,1))
  plot(c(), xlim=c(1,nc), ylim=c(1,nr), main='Convex shapes in image', xlab=NULL, ylab=NULL, xaxt='n', yaxt='n')
  if (length(convex.shapes) >= 1) {
    for (idx in 1:length(convex.shapes)) {
      shape <- convex.shapes[[idx]]
      s.x <- c(shape$hull.x, shape$hull.x[1])
      s.y <- nr + 1 - c(shape$hull.y, shape$hull.y[1])
      lines(s.x, s.y, col=shape$hull.color, type='l')
      transparent.color <- do.call(rgb, as.list(c(col2rgb(shape$hull.color),170)/255))
      polygon(s.x, s.y, border=shape$hull.color, col=transparent.color)
    }
  }
  par(default.par)  # restore graphical parameters
}


# Given a subimage of a connected shape and associated RGB values,
# return the conver hull information, possible adjusted by anchor point
# that is the top-left point in the subimage.
extractConvexShape <- function(shape.rectangular, subimg.rgb, anchor.row=1, anchor.col=1) {
  sr <- nrow(shape.rectangular)
  sc <- ncol(shape.rectangular)
  x.mask <- matrix(rep(1:sc, each=sr), nrow=sr)
  y.mask <- matrix(rep(-1:-sr, sc), nrow=sr)
  # convert pixels into points in Euclidean space (use pixel centers)
  x <- x.mask[shape.rectangular]
  y <- y.mask[shape.rectangular]
  ch <- findConvexHull2D(x, y)
  inside.ch <- insideConvexHull2D(x.mask, y.mask, x[ch], y[ch])
  shape.size <- sum(shape.rectangular)
  hull.avg.col <- unlist(lapply(
    list(RED,GREEN,BLUE),
    function(channel) mean(extractRGBChannel(subimg.rgb, channel)[shape.rectangular])
  ))
  return(list(
    hull.x= x[ch] + anchor.col - 1,
    hull.y= -y[ch] + anchor.row - 1,
    hull.color=do.call(rgb, as.list(hull.avg.col)),
    hull.compsize=shape.size,
    hull.coverage=shape.size/sum(inside.ch)
  ))
}

# Datta measure 56, "shape convexity"
###########
shapeConvexity <- function(conn.components, img.rgb) {
  nr <- nrow(conn.components)
  nc <- ncol(conn.components)
  tab <- tabulate(conn.components)
  names(tab) <- 1:length(tab)
  ordered.components <- rev(sort(tab))
  threshold <- nr*nc/200
  comps <- ordered.components[ordered.components >= threshold]
  bound.boxes <- getBoundingBoxes(conn.components, length(tab))
  convex.shapes <- list()
  if (length(comps) >= 1) {
    for (idx in 1:length(comps)) {
      id <- as.integer(names(comps[idx]))
      box <- bound.boxes[id, ]
      bounding.rows <- box$FirstRow:box$LastRow
      bounding.cols <- box$FirstColumn:box$LastColumn
      rectangular <- conn.components[bounding.rows, bounding.cols, drop=FALSE]
      shape.rectangular <- rectangular == id
      subimg.rgb <- img.rgb[box$FirstRow:box$LastRow, box$FirstColumn:box$LastColumn, , drop=FALSE]
      convex.shapes[[idx]] <- extractConvexShape(shape.rectangular, subimg.rgb, box$FirstRow, box$FirstColumn)
    }
  }
  DO_VIEW_CONVEX <- FALSE  # TRUE
  if (DO_VIEW_CONVEX) { plotConvexShapes(dim(img.rgb), convex.shapes) }
  CONVEX_THRESHOLD <- 0.8
  incl.mask <- unlist(lapply(convex.shapes, function(x) x$hull.coverage >= CONVEX_THRESHOLD))
  shape.convexity.feature <- if (length(comps)) sum(as.integer(comps[incl.mask])) / (nr*nc) else 0
  return(shape.convexity.feature)
}


segmentBlockLocations <- function(conn.components, largest.ids) {
  seg.idx <- lapply(largest.ids, function(x) which(conn.components == x))
  names(seg.idx) <- largest.ids
  BLOCKS <- 3
  nr <- nrow(conn.components)
  nc <- ncol(conn.components)
  seg.avg.row <- lapply(seg.idx, function(idx) mean(((idx-1) %% nr) + 1))
  seg.avg.col <- lapply(seg.idx, function(idx) mean(((idx-1) %/% nr) + 1))
  block.row <- findInterval(seg.avg.row, seq(1, nr, length.out=BLOCKS+1), rightmost.closed=TRUE)
  block.col <- findInterval(seg.avg.col, seq(1, nc, length.out=BLOCKS+1), rightmost.closed=TRUE)
  block.code <- 10 * block.row + block.col  # Datta 48--52

  # My own location variant: distance from center, max 1.0
  sr <- unlist(seg.avg.row)
  sc <- unlist(seg.avg.col)
  mid.r <- (nr + 1) / 2
  mid.c <- (nc + 1) / 2
  max.dist <- sqrt((mid.r-1)^2 + (mid.c-1)^2)
  avg.dist <- sqrt((sr-mid.r)^2 + (sc-mid.c)^2)
  block.deviation <- avg.dist / max.dist  # Esa's proxy to Datta 48--52

  return(list(codes=block.code, distances=block.deviation))
}


# Datta measures 24--52 based on image segmentation (except 46 & 47)
## Datta measure 56, "shape convexity", included

#' Image Region Composition Features (Datta measures 24--45 & 48--52 & 56)
#' 
#' Compute regional \emph{composition} features of an image, as in
#' \emph{Datta measures} 24--45 & 48--52 & 56.
#' 
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' 
#' @param img.rgb photo as an RGB image array (\emph{m} x \emph{n} x 3)
#' @param num.segments maximum number of image segments to identify and analyze.
#' @return numeric \emph{compositional} region features as a list:
#'   \itemize{
#'     \item \code{"num.large.patches"}: number of segments that contain
#'       more than 1\% of image's pixels (Datta 24)
#'     \item \code{"num.clusters"}: optimal number of segments found (Datta 25)
#'     \item \code{"avg.patch.hsv"}: average values on HSV dimensions of
#'       largest segments, or \code{NA} if non-existent (Datta 26--40)
#'     \item \code{"rel.patch.sizes"}: proportional size of each of the largest
#'       segments w.r.t. full image size (Datta 41--45)
#'     \item \code{"segment.positions"}: coded locations of largest segments,
#'       determined by which cell in a 3x3 grid the segments' mass-centers hit
#'       (Datta 48--52): 11, 12, 13, 21, 22, 23, 31, 32, or 33.
#'     \item \code{"segment.distances"}: Esa Junttila's proxy for Datta 48--52:
#'       distances from image midpoint to segments' mass centers,
#'       normalized onto \code{[0;1]}
#'     \item \code{"shape.convexity"}: proportion of image covered by
#'     almost-convex segments, having sizes at least 0.5\%
#'     and that compose at least 80\% of their convex hull shapes (Datta 56)
#'   }
#' 
#' @export
regionCompositionFeatures <- functionregionCompositionFeatures <- function(img.rgb, num.segments=5) {
  nr.img <- nrow(img.rgb)
  nc.img <- ncol(img.rgb)
  seg <- photoSegmentation(img.rgb)  # used 'k.min=num.segments' but number of clusters is sometimes smaller
  conn.components <- seg$conn.components
  k <- seg$optimal.k

  # Compute Datta measures:
  freqs <- tabulate(conn.components)  # "tabulate" is much faster than "table"
  names(freqs) <- 1:length(freqs)
  comp.sizes <- rev(sort(freqs))
  largest <- head(comp.sizes, num.segments)
  largest.ids <- as.integer(names(largest))
  num.threshold.comps <- sum(largest / (nr.img * nc.img) > 0.01)  # Datta feature 24
  num.clusters <- k  # Datta feature 25
  hsv.pixels <- matrix(toHSV(img.rgb), ncol=3)
  # Make sure we have enough HSV values, even if they are NA:
  largest.hsv.avg <- matrix(NA, nrow=num.segments, ncol=3)
  largest.hsv.avg[1:length(largest.ids), 1:3] <- t(
    sapply(
      largest.ids,
      function(x) colMeans(hsv.pixels[conn.components == x, , drop=FALSE], na.rm=TRUE)
  ))
  #print(largest.hsv.avg)  # Datta 26--40
  rel.sizes <- as.numeric(largest / (nr.img * nc.img)) # Datta 41--45

  ## Skipping these two measures, as their definitions seem unclear and non-intuitive.
  #avg.col.spread <-  # Datta 46
  #avg.col.complmentary <-   # Datta 47

  # Datta 48--52 and Esa's proxy measures (distances from center)
  locations <- segmentBlockLocations(conn.components, largest.ids)

  convexity <- shapeConvexity(conn.components, img.rgb)  # Datta 56

  return(list(
    num.large.patches=num.threshold.comps,  # Datta 24
    num.clusters=num.clusters,  # Datta 25
    avg.patch.hsv=largest.hsv.avg,  # Datta 26--40; hues may be NA for pure black and white
    rel.patch.sizes=rel.sizes,  # Datta 41--45
      # excluding Datta 46 & 47
    segment.positions=locations$codes,  # Datta 48--52
    segment.distances=locations$distances,  # My own proxy for Datta 48--52: distances from center
    shape.convexity=convexity  # Datta 56
  ))
}


# Datta measures 53--55: Depth-of-field, DOF
###########

# I'm not quite sure how we should aggragete the different coefficients: wLH, HL, and HH.
# Now summing them together with the rest of the values.

waveletDattaDof <- function(img.channel) {
  d <- waveletTexture(img.channel)  # Already defined in Datta 10
  num.levels <- attr(d, 'J')
  if (num.levels < 3) return(NA)
  nr3 <- nrow(d$HH3)
  nc3 <- ncol(d$HH3)
  BLOCKS <- 5
  r.idx <- round(seq(1, nr3, length.out=BLOCKS+1)) # split into BLOCKS by endpoints
  c.idx <- round(seq(1, nc3, length.out=BLOCKS+1))
  r.mid <- r.idx[2]:r.idx[BLOCKS]  # mid-part, exclude first and last block
  c.mid <- c.idx[2]:c.idx[BLOCKS]
  expectation <- length(r.mid)*length(c.mid) / (nr3 * nc3)
  coeffSum3 <- function(ri, ci) sum(abs(d$HH3[ri,ci]) + abs(d$HL3[ri,ci]) + abs(d$LH3[ri,ci]), na.rm=TRUE)

  # Debug:
  #m <- abs(d$HH3[1:nr3,1:nc3]) + abs(d$HL3[1:nr3,1:nc3]) + abs(d$LH3[1:nr3,1:nc3])
  #view((m - min(m)) / (max(m)-min(m)), title='aggregated abs features')

  # I decided to make the "dof.indicator" measure relative.
  # It is the deviation from expectation, measured by expectation-units:
  # * uniform photo gets ~0, center-dominant gets >0, edge-dominant gets <0
  # * Maximum value is 1/(3^2/5^2) - 1 = 1.778 when all detail is in center; minimum is -1 (no detail in center).
  full.normalizer <- coeffSum3(1:nr3, 1:nc3)
  if (full.normalizer == 0) dof.indicator <- 0  # in extreme cases we assume uniform image
  else dof.indicator <- (coeffSum3(r.mid, c.mid) / full.normalizer - expectation) / expectation
  #print(dof.indicator)
  return(dof.indicator)
}


#' Depth-of-field Features of an Image (Datta measures 53--55)
#' 
#' Compute depth-of-field (DOF) features in the image,
#' as in \emph{Datta measures} 53--55, with Esa Junttila's
#' relative normalization. It measures if there is any extra detail in the
#' center part, compared to the rest of the image.
#' 
#' For description of Datta measures, refer to \url{https://doi.org/10.1007/11744078_23}
#' 
#' @param img.hsv image as an HSV image array (\emph{m} x \emph{n} x 3)
#' @return List of three numeric \emph{DOF} features, one for each of
#'   HSV channels. The returned values are the deviations from expectation,
#'   measured by expectation-units:
#'   \itemize{
#'     \item uniform photo gets ~0,
#'     \item center-dominant gets >0,
#'     \item edge-dominant gets <0.
#'   }
#'   Maximum DOF value is 1/(3^2/5^2) - 1 = 1.778 when all detail is in center;
#'   minimum is -1 (no detail in center).
#' @examples
#' set.seed(1)
#' img <- array(runif(40*60*3, 0.4, 0.6), dim=c(40,60,3))  # uniform starting-point
#' img[15:25, 20:40, 1:3] <- rnorm(11*21*3, 0.5, 0.15)  # concentration at center
#' all(abs(unlist(dof(toHSV(img))) - c(-0.1018337, 0.2587495, 0.1880983)) < 1e-5)
#' 
#' @export
dof <- function(img.hsv) {
  img.hsv[is.na(img.hsv)] <- 0  # local changes; DWT works only non-NA values
  res <- lapply(HSV, function(channel) waveletDattaDof(extractHSVChannel(img.hsv, channel)))
  return(list(
    hue.dof=res[[HUE$idx]],
    sat.dof=res[[SATURATION$idx]],
    val.dof=res[[VALUE$idx]]
  ))  # not quite sure how we should compute the DOF measure
}




# References, links and tests for LUV and EMD.

#rgb.centers=as.list(sapply(, function(x) ))
#names(rgb.centers) <- sapply(0:(n.buckets-1), toString)

# Convert to LUV color space
# http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html
# http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
# ?convertColor
# http://pages.cs.wisc.edu/~dyer/cs766/hw/hw2/code/rgb2luv.m
# https://www.easyrgb.com/en/math.php
# http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Luv.html


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


# Test example:
# Colors, seven of them, by names:
#   orange        red           red     dark-purple  black
#   light-orange  orange        red     dark-purple  dark-blue
#   pale-white    light-orange  orange  red          dark-purple
# Colors by RGB and related buckets(n=4):
#   color name     RGB value     bucketcodes (0--3, -->r+g*n+b*n*n+1) from 1 to 64
#   orange       = 255,130,  0   320 -> 3+4*2+4*4*0+1 = 12
#   light-orange = 255,200,100   331 -> 3+4*3+4*4*1+1 = 32
#   pale-white   = 255,240,220   333 -> 3+4*3+4*4*3+1 = 64
#   red          = 255,  0,  0   300 -> 3+4*0+4*4*0+1 =  4
#   dark-purple  = 125,  0,125   101 -> 1+4*0+4*4*1+1 = 18
#   dark-blue   =  60, 65,110   011 -> 0+4*1+4*4*1+1 = 21
#   black        =   0,  0,  0   000 -> 0+4*0+4*4*0+1 =  1
# Pixel RGB buckets:
#   12  4  4 18  1
#   32 12  4 18 21
#   64 32 12  4 18
# Bucket frequencies:
#    1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64
#    1  0  0  4  0  0  0  0  0  0  0  3  0  0  0  0  0  3  0  0  1  0  0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1
# Bucket center RGB values:
#centralRGB <- function(x, n) {
#  b <- (x-1) %/% n^2
#  g <- ((x-1) - b*n^2) %/% n
#  r <- ((x-1) - b*n^2 - g*n)
#  red <- mean(c(0, 1/n) + r/n)
#  green <- mean(c(0, 1/n) + g/n)
#  blue <- mean(c(0, 1/n) + b/n)
#  return(c(red,green,blue))
#}
#n <- 4
#centrals <- array(NA, dim=c(3,5,3), dimnames=list(NULL, NULL, c('R','G','B')))
#px <- matrix(c(
#  12,  4,  4, 18,  1,
#  32, 12,  4, 18, 21,
#  64, 32, 12,  4, 18), nrow=3, byrow=TRUE)
#for (row in 1:nrow(px)) {
#  for (col in 1:ncol(px)) {
#    px.rgb <- centralRGB(px[row,col], n)
#    for (i in 1:3) centrals[row,col,i] <- px.rgb[i]
#  }
#}
#source('viewer.R')
#view(centrals) # checked colors visually and with a color picker tool
# Transform to LUV with http://colormine.org/convert/rgb-to-luv, for example:
# Bucket  RGB-central  LUV-central
#  1       32, 32, 32  12.2500301015228,   0.0001091072937017, −0.0020957293787555
#  4      223, 32, 32  48.0391344061145, 142.982255605226,     30.8335943777994
# 12      223,159, 32  69.8866725638048,  53.0398933951075,    67.8563682857694
# 18       96, 32, 96  24.8411338297045,  25.6808799074132,   −33.2072534707568
# 21       32, 96, 96  36.9618665690869, −24.3501763192539,    −5.25883520398827
# 32      223,223, 96  86.65626930539,     5.5432074756774,    76.8426762606292
# 64      223,223,223  88.8236315210308,   0.0007911250806805, −0.0151959050356324
#
# What we got from our code, example. There are some differences because
# of confusing color profiles, white points, and what not, so not exactly correct.
# LUV-centers
# [ 1,] 12.18863   0.001073099  9.251800e-04
# [ 4,] 48.06164 143.113620107  3.087699e+01
# [12,] 69.99378  52.850982127  6.802831e+01
# [18,] 24.73662  25.559904896 -3.303788e+01
# [21,] 36.81698 -24.246935292 -5.228755e+00
# [32,] 86.69642   5.567193840  7.705239e+01
# [64,] 88.86786   0.007824012  6.745527e-03
#img <- readImage('../examples/colorfulness-test.png')
#print(colorfulness(img))



# Unnecessary, but useful, code used for testing, viewing and profiling:
####

# Test image                                           EMD colorfulness distance and bucket entropy (max 6)
#img <- readImage('examples/uniform-buckets.png')   #  0  entropy 6 b (minimum,maximum)
#img <- readImage('examples/many_colors.png')       # 18
#img <- readImage('examples/small_grid.png')        # 48  entropy 3.14 b
#img <- readImage('examples/niemi.png')             # 49  entropy 3.67 b
#img <- readImage('examples/no_shift.png')          # 57
#img <- readImage('examples/K5_10994.JPG')          # 58
#img <- readImage('examples/puffin.jpg')            # 61
#img <- readImage('examples/colorfulness-test.png') # 64
#img <- readImage('examples/dark_city.png')         # 65
#img <- readImage('examples/almost_black.png')      # 83
#img <- readImage('examples/sharp_or_blur.png')     # 83
#img <- readImage('examples/bluehue.png')           # 86  entropy 0.60 b
#img <- readImage('examples/pure-red.png')          # 153.7  entropy 0 b (maximum,minimum)
#img <- readImage('examples/grainy.jpg')            #
#img <- readImage('examples/datta-colorfulness-high-2.png')
#img <- readImage('examples/green_grass_blue_sky.png')

# Mono-colored images have maximum distances from/to uniform distribution.
# Ranging from 75.0 (bucket 43/64: light-grey) up to 153.7 (bucket 4/64: bright red).

# Example images from the article do not behave like the authors claim:
#img <- readImage('examples/datta-colorfulness-high-1.png')  # 70  entropy 2.56 b
#img <- readImage('examples/datta-colorfulness-high-2.png')  # 63  entropy 2.50 b
#img <- readImage('examples/datta-colorfulness-low-1.png')   # 71  entropy 2.08 b
#img <- readImage('examples/datta-colorfulness-low-2.png')   # 57  entropy 2.25 b

#source('viewer.R')
#view(img)
#print(toXYZ(img))

#Rprof(filename="Rprof.out", append = FALSE, interval = 0.01,
#      memory.profiling = FALSE, gc.profiling = FALSE,
#      line.profiling=TRUE, numfiles = 100L, bufsize = 10000L)

#print(colorfulness(img, 'uniform'))
#print(colorfulness(img, 'uniform', n=6))
#print(colorfulness(img, 'grey', n=6))

#print(texture(img))

# When debugging, enable the folloing line:
#print(regionCompositionFeatures(img))


#img <- readImage('../examples/K5_10994.JPG')
#n <- 10
#start.time <- Sys.time()
#replicate(n, dim(toXYZ(img)))
#end.time <- Sys.time()
#print((end.time - start.time) / n)

