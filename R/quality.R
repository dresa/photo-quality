##
## Quality measures for photos and images.
##
## Esa Junttila, 2016-03-23


source('R/common.R')  # limit function

# The Blur Effect: Perception and Estimation with a New No-Reference Perceptual Blur Metric"
# Crete F., Dolmiere T., Ladret P., Nicolas M. - GRENOBLE - 2007
# In SPIE proceedings - SPIE Electronic Imaging Symposium Conf Human Vision and Electronic Imaging, tats-Unis d'Amrique (2007)
# https://hal.archives-ouvertes.fr/hal-00232709/document
# Error is possible (div by 0 in normalization)?
# Minimum value 0, best quality in terms of blur perception.
# Maximum value 1.0 is reported. Is it reachable? Checkerboard image: sum.gray=1.0*h*w, sum.cmp=h*w*0.5 yields 0.5 (approx)
# http://www.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric/content//blurMetric.m

blurAnnoyanceQuality <- function(img, f.len=9) {
  h <- dim(img)[1]
  w <- dim(img)[2]
  # Sharpness requires only one channel: luminance
  img.gray <- luminance(img)
  # Vertical and horizontal filters:
  f.vert <- matrix(1/f.len, nrow=f.len)
  f.horiz <- t(f.vert)
  blurred.vert <- imfilter(img.gray, f.vert, pad='replicate')[ , , 1]
  blurred.horiz <-imfilter(img.gray, f.horiz, pad='replicate')[ , , 1]
  # Differences in neighboring pixels (both within original and blurred):
  diff.gray.vert     <- abs(    img.gray[2:h, 1:w] -     img.gray[1:(h-1), 1:w])
  diff.blurred.vert  <- abs(blurred.vert[2:h, 1:w] - blurred.vert[1:(h-1), 1:w])
  diff.gray.horiz    <- abs(     img.gray[1:h, 2:w] -      img.gray[1:h, 1:(w-1)])
  diff.blurred.horiz <- abs(blurred.horiz[1:h, 2:w] - blurred.horiz[1:h, 1:(w-1)])
  # Difference between original and blurred pixels:
  cmp.vert <- pmax(0, diff.gray.vert - diff.blurred.vert)  # including all columns (not just 2:w)
  cmp.horiz <- pmax(0, diff.gray.horiz - diff.blurred.horiz)  # including all rows (not just 2:h)
  # Comparison between original (gray) and blurred: sum of differences
  sum.gray.vert <- sum(diff.gray.vert)
  sum.cmp.vert <- sum(cmp.vert)
  sum.gray.horiz <- sum(diff.gray.horiz)
  sum.cmp.horiz <- sum(cmp.horiz)
  # Normalization
  blur.vert <- (sum.gray.vert - sum.cmp.vert) / sum.gray.vert
  blur.horiz <- (sum.gray.horiz - sum.cmp.horiz) / sum.gray.horiz
  blur <- max(blur.vert, blur.horiz)  # from 0 to 1.0, smaller is better
  # Convert into a quality parameter (from 1 to 5, larger is better)
  # Constants come from the interpolation model parameters, as reported in the article.
  quality <- 3.79/(1 + exp(10.72*blur - 4.55)) + 1.13  # from 1 to 5, larger is better
  ## My own normalization
  normalized <- (quality-1)/4  # convert to 0--1 (larger is better)
  return(normalized)
}


# Helper function for the 'mdwe' method -- not a universal method.
# TODO: how to derive the magic constant? How many bins do we consider?
detectVerticalEdgeAreas <- function(img.gray) {
  f.vert <- sobelVertical()
  filtered <- imfilter(img.gray, f.vert, pad='replicate')[ , , 1]
  edge.strengths <- abs(filtered)
  bins <- 20  # granularity: just my own magic constant
  breaks <- seq(min(edge.strengths), max(edge.strengths), length=bins)
  h <- hist(c(edge.strengths), breaks=breaks, plot=FALSE)
  thres.idx <- otsuThreshold(h$counts)
  thres.val <- h$mids[thres.idx]
  edges <- edge.strengths >= thres.val
  return(edges)
}
# Returns a vector of midpoints (rounding up) for consecutive sub-sequences of TRUE.
# Input: logical vector
# For example: TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE  --> 1,5,9
# TODO: consider whether consecutive edges should allow one index without
#       edge indication, such as ...,0,0,0,1,1,1,1,0,1,1,1,0,0,...
sequenceMidpoints <- function(vec) {
  idxs <- which(vec)  # indices that are TRUE
  starts <- idxs[which(diff(c(-1, idxs)) != 1)]    # start indices of consecutive TRUEs, first TRUE included
  ends <- idxs[which(diff(c(idxs, length(vec)+2)) != 1)]  # end indices of consecutive TRUEs, last TRUE included
  mids <- as.integer(ceiling(starts + (ends - starts)/2))  # mid indices of consecutive TRUE sub-sequences
  return(mids)
}
localMinimaIndices <- function(vec) {
  n <- length(vec)
  d <- c(0, diff(vec), 0)
  return(which(d[1:n] <= 0 & d[2:(n+1)] >= 0))
}
localMaximaIndices <- function(vec) {
  n <- length(vec)
  d <- c(0, diff(vec), 0)
  return(which(d[1:n] >= 0 & d[2:(n+1)] <= 0))
}

# Pina Marziliano, Frederic Dufaux, Stefan Winkler and Touradj Ebrahimi:
# A No-Reference Perceptual Blur Metric, 2002
# TODO: derive a model to produce a quality index between 1 and 5, for example (based on data in the article)
# TODO: horizontal-edge version of the method
mdweVertical <- function(img) {
  img.gray <- luminance(img)  # use only the luminosity image
  edges.areas <- detectVerticalEdgeAreas(img.gray)  # several consecutive pixels may be marked as belonging to an edge
  total.width <- 0
  num.edges <- 0
  for (r in 1:nrow(img)) {  # for each row
    # Find edge indices:
    edge.mids <- sequenceMidpoints(edges.areas[r, ])
    # Find all local maxima and minima for the luminosity vector r:
    horiz.vec <- img.gray[r, ]  # luminosities for row pixels
    minima.idxs <- localMinimaIndices(horiz.vec)
    maxima.idxs <- localMaximaIndices(horiz.vec)
    # For each edge mid-point, find a local minimum and maximum so that mid-point is enclosed:
    min.idxs <- findInterval(edge.mids, minima.idxs, all.inside=TRUE)
    max.idxs <- findInterval(edge.mids, maxima.idxs, all.inside=TRUE)
    prev.mins <- minima.idxs[min.idxs]  # preceding local minimum
    prev.maxs <- maxima.idxs[max.idxs]
    # From previous minimum and maximum, deduce the next minimum or maximum.
    # It is either the next index or the same index (when prev idx is last or edge-mid is at a local optimum).
    not.same.mins <- (edge.mids != prev.mins)*(length(minima.idxs)+1)  # impossibly large when edge.mid != prev.min
    not.same.maxs <- (edge.mids != prev.maxs)*(length(maxima.idxs)+1)  # impossibly large when edge.mid != prev.max
    next.mins <- minima.idxs[pmin(min.idxs + 1, length(minima.idxs), prev.mins + not.same.mins)]  # next, last, or same
    next.maxs <- maxima.idxs[pmin(max.idxs + 1, length(maxima.idxs), prev.maxs + not.same.maxs)]  # next, last, or same
    # Compute widths for all edge intervals on this row:
    widths <- pmin(next.mins - prev.maxs,  next.maxs - prev.mins)  # choose the minimum width: maxmin or minmax
    total.width <- total.width + sum(widths)
    num.edges <- num.edges + length(edge.mids)
  }
  score.raw <- total.width/num.edges  # Inf when num.edges==0
  ## My own normalization
  score.gaussian <- limit(score.raw^1.04 - 4.5, 1, 10)  # From no blur (1) to distractive blur (10)
  score.jpeg2k <- limit(17.5*max(0,score.raw-3)^0.3 - 19.5, 1, 10)
  norm <- function(x) -((score.gaussian - 1) / (10-1) - 1)  # 1--10 (smaller is better) --> 0--1 (larger is better)
  return(list(score=score.raw, gaussian=norm(score.gaussian), jpeg2k=norm(score.jpeg2k)))
}


# A New No-reference Method for Color Image Quality Assessment
# Sonia Ouni, Ezzeddine Zagrouba, Majed Chambah, 2012
# International Journal of Computer Applications (0975 ? 8887)
# Volume 40? No.17, February 2012

# Note that there are errors in the RGB->HSV conversion formula. See Wikipedia instead.
# The article contains numerous errors, the examples don't match with definitions,
# and the exact methods are sometimes vague. I've done my best to develop a method
# that is both practical and close to the spirit of the authors' ideas.

# COMMENTS:

## Kappa is the same as the following (measure of "compactness"):
#C <- sum(cos(hues[mask]))
#S <- sum(sin(hues[mask]))
#R <- sqrt(C^2 + S^2)
#n <- length(hues[mask])
#kappa <- A1inv(R/n)

## The exact computation of DS leads to a huge computational task.
#k <- length(row.idx)
#ds <- 0
#for (i in 1:k) {
#  local.dist <- sum(sqrt((row.idx - row.idx[i])^2 + (col.idx - col.idx[i])^2)) / k
#  ds <- ds + local.dist
#}
#spatial.dispersion <- ds / k
## OR
#spatial.dispersion <- sum(sapply(1:length(row.idx), function(i) sqrt((row.idx - row.idx[i])^2 + (col.idx - col.idx[i])^2))) / npdc^2

## This is an average distance between pixels that have dominant color.
## I still think the value should be normalized somehow, for example dividing
## by an average distance between two pixels in the whole image, not just
## in dominant color. Otherwise image size has an effect on the absolute
## magnitude of the number. In the original article the values
## like 0.0020936 don't quite match with the definitions.

dispersionDominantColor <- function(img.hsv) {
  hues <- degreeToRadian(extractHSVChannel(img.hsv, HUE))  # convert from degrees to radians
  nr <- nrow(hues)
  nc <- ncol(hues)
  mask <- !is.na(hues)  # existing hues, excluding white and black

  #is.almost.gray <- sum(mask) < 0.01*nr*nc  # only a minority of pixels has a hue (is not grey)
  #if (is.almost.gray) return(list(mu=NA, kappa=NA, pi=NA, ds=NA, custom.ds=NA))

  mu <- atan2(sum(sin(hues[mask])), sum(cos(hues[mask])))%%(2*pi)  # angles between hues
  # Estimate kappa parameter on a von Mises distribution in a circular domain.
  kappa <- A1inv(mean(cos(hues[mask] - mu)))
  
  # NPDC and pi measures
  ######################
  # distances between hues (in radians)
  dist.raw <- abs(hues[mask] - mu)
  dist <- pmin(dist.raw, 2*pi - dist.raw)  # distance between hues, mod 2*pi
  npdc.mask <- dist <= 1/kappa  # article has a different interpretation of 1/kappa ??
  npdc <- sum(npdc.mask)  # pixels in a dominant color (within range kappa)

  # If we are too close to a grayscale image, then it's time to give up:
  if (npdc < 4) return(list(mu=NA, kappa=NA, pi=NA, ds=NA, custom.ds=NA))

  pi.measure <- npdc / length(hues[mask])
  px.row <- matrix(1:length(hues) %% nr, nrow=nr)
  px.col <- matrix(1:length(hues) %/% nr, nrow=nr)
  
  # Spatial dispersion, by approximation
  ######################################
  row.idx <- c(px.row[mask][npdc.mask])
  col.idx <- c(px.col[mask][npdc.mask])
  npdc.nr <- length(row.idx)
  npdc.nc <- length(col.idx)

  # Approximation for all Euclidean distances, for custom normalization
  num.samples <- min(10000, nr*nc/2)  # just enough to get a hint about the avg distance
  smpl <- function(values) sample(values, num.samples, replace=TRUE)
  approx.dist <- mean(sqrt((smpl(1:nr) - smpl(1:nr))^2 + (smpl(1:nc) - smpl(1:nc))^2))
  
  # NPDC Euclidean distance approximation by random sample
  px.1 <- smpl(1:npdc.nr)
  px.2 <- smpl(1:npdc.nr)
  npdc.approx.dist <- mean(sqrt((row.idx[px.1] - row.idx[px.2])^2 + (col.idx[px.1] - col.idx[px.2])^2))

  spatial.dispersion <- npdc.approx.dist
  custom.ds <- min(1, npdc.approx.dist / approx.dist)  # cap to 1 (exceeds by random variation)
  
  return(list(mu=mu, kappa=kappa, pi=pi.measure, ds=spatial.dispersion, custom.ds=custom.ds))
}



# A set of basic measures, as described in slideshow names
# Automatic photo quality assessment --- Taming subjective problems with hand-coded metrics

basicExposureLevel <- function(img) {
  # Using mean pixel intensity (interpreted via luminance) to determine
  # whether we have a "correct" exposure. An underexposured image has mean
  # luminance close to 0.0; that of an overexposured image is closer to 1.0.
  # Uses a scaled Beta distribution to set an exposure score (highest value
  # 1.0 at 0.5), so it punishes for underexposure and overexposure.
  SHAPE = 2.0
  rawScore <- function(x) dbeta(x, SHAPE, SHAPE)
  m <- mean(luminance(img))  # 0 is black, 1 is white
  exposure.score <- rawScore(m) / rawScore(0.5)  # normalize by maximum density
  return(exposure.score)
}

basicRmsContrast <- function(img) {
  # RMS contrast
  # https://en.wikipedia.org/wiki/Contrast_(vision)
  # i <- luminance(img)
  # return(sqrt(sum((i - mean(i))^2) / (nrow(i)*ncol(i))))
  return(sd(luminance(img)))
}

basicIntervalContrast <- function(img, interval=0.95) {
  # Slides
  tail.prob <- (1 - interval) / 2
  i <- sort(luminance(img))
  q <- quantile(i, c(tail.prob, 1 - tail.prob))
  return(as.numeric(diff(q)))
}

