##
## Collection of no-reference quality measures for photos and images.
##
## Esa Junttila, 2016-03-23


source('R/common.R')  # limit function

#' Blur annoyance measure.
#' 
#' \emph{Blur annoyance} measures the quality of photos
#' w.r.t perceptual blurriness. The method uses
#' a grayscaled version of given photo and applies both horizontal
#' \code{h = [1/n 1/n ... 1/n 1/n]} and
#' vertical filters to make the photo blurry. It then compares the
#' original and blurred images to measure the differences in
#' blurriness. The larger the perceptual difference, the less blurry
#' (and better in quality) the original photo is.
#' The blur quality in the original photo is then measured by a scale
#' from 0 to 1 (larger is better), based on an empirically-derived
#' non-linear conversion.
#' 
#' The method was introduced originally in:
#' 
#' Frédérique Crété-Roffet, Thierry Dolmiere, Patricia Ladret, and Marina Nicola.
#' "The Blur Effect: Perception and Estimation with a New No-Reference
#' Perceptual Blur Metric."
#' \emph{SPIE Electronic Imaging Symposium Conf Human Vision and Electronic Imaging},
#' Jan 2007, San Jose, United States. XII, pp.EI 6492-16, 2007.
#' URL: \url{https://hal.archives-ouvertes.fr/hal-00232709/document}
#' 
#' For some pathological photos the Blur Annoyance may be undefined,
#' returning \code{NaN}. Note that blurriness is normalized by the differences
#' in neighboring pixels in original photo, so if the photo is "flat"
#' there are no differences, which leads to division by 0.
#' 
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3).
#' @param f.len [integer], strength of the blur filter (length of
#' horizontal/vertical filters), default is 9 pixels in width.
#' @return numeric blur annoyance between 0 (worst) and 1 (best quality),
#' or \code{NaN} for pathological photos.
#' @seealso
#' There exists a reference implementation in Matlab.
#' Note that the original article describes a non-linear model for converting
#' blur measures \code{[0;1]} into blur annoyance qualities \code{[1;5]}.
#' The referred Matlab implementation, however, leaves the conversion out.
#' This R implementation has the conversion included, but also rescales the
#' range of blur annoyance from \code{[1;5]} into \code{[0;1]}.
#' \url{http://www.mathworks.com/matlabcentral/fileexchange/24676-image-blur-metric/content//blurMetric.m}
#' @examples
#' a <- array(seq(0,1,length.out=24), dim=c(2,4,3))
#' abs(blurAnnoyanceQuality(a) - 0.7208079) < 1e-6
#' r <- array(runif(100*200), dim=c(100,200,3))
#' abs(blurAnnoyanceQuality(r) - 0.95) < 1e-3
#' b <- array(0.6, dim=c(2,4,3))
#' is.nan(blurAnnoyanceQuality(b))
#' 
#' @export
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
  blur <- max(blur.vert, blur.horiz)
  # Minimum blur value 0 has best quality in terms of blur perception.
  # Maximum blur value 1 has worst blur perception. Maximum value 1 is
  # reported in the article, but how can it be reached? For example a
  # checkerboard photo has: sum.gray=1.0*h*w, sum.cmp=h*w*0.5 yields 0.5 (approx).
  
  # Convert into a quality parameter (range 1 to 5, larger is better)
  # Constants come from the interpolation model parameters, as reported in the article.
  quality <- 3.79/(1 + exp(10.72*blur - 4.55)) + 1.13  # from 1 to 5 (approx)
  ## My own normalization back to [0;1]
  normalized <- (quality-1)/4  # convert to 0--1 (larger is better)
  return(normalized)
}


# Helper function for the 'mdwe' method -- not a universal method.
# TODO: how to derive the magic constant? How many bins do we consider?
detectVerticalEdgeAreas <- function(img.gray, bins=20) {
  f.vert <- sobelVertical()
  filtered <- imfilter(img.gray, f.vert, pad='replicate')[ , , 1]
  edge.strengths <- abs(filtered)
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


#' MDWE Blur Measure.
#' 
#' MDWE is a measure of photo blurriness, originally labeled as
#' \emph{perceptual blur metric}. It is based on the analysis
#' of the spread of the vertical edges in an grayscale version of an image.
#' The MDWE Blur Measure maps a photo onto a quality value \code{[0;1]},
#' where larger value represents better blur quality in the photo.
#' The method measures quality on both Gaussian type blur and
#' JPEG200 artefacts.
#' 
#' The method was originally introduced in:
#' Pina Marziliano, Frederic Dufaux, Stefan Winkler and Touradj Ebrahimi.
#' A No-Reference Perceptual Blur Metric.
#' \emph{IEEE International Conference on Image Processing}, 2002, pages 57--60.
#' 
#' URL: \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.7.9921&rep=rep1&type=pdf}
#' 
#' The original article provides both subjective (expert) and
#' objective (computed) blur ratings for a small sample of photos.
#' For this R implementation I developed a simple estimator that predicts
#' a subjective blur rating from a MDWE raw score. The original range of
#' \code{[1;10]} for MDWE raw scores was also flipped and rescaled.
#' 
#' Future work includes developing an horizontal-edge version of the method,
#' and combining vertical and horizontal measures.
#' 
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3).
#' @return a list of edge-spanning blur quality values:
#'   \itemize{
#'     \item \code{"score"}:
#'       original blur score as described in the article
#'     \item \code{"gaussian"}:
#'       quality on Gaussian-type blur within \code{[0;1]} (larger is better).
#'     \item \code{"jpeg2k"}:
#'       quality on JPEG2000 artefacts within \code{[0;1]} (larger is better).
#'   }
#'   
#' @examples
#' set.seed(42)
#' r <- array(runif(100*200), dim=c(100,200,3))
#' all(abs(unlist(mdweVertical(r)) - c(1.917252, 1, 1)) < 1e-5)
#' b <- imfilter(r, gaussianFilter2D(10,4))
#' all(abs(unlist(mdweVertical(b)) - c(7.95, 0.6514027, 0.1359882)) < 1e-5)
#' 
#' @export
mdweVertical <- function(img) {
  img.gray <- luminance(img)  # use only the luminosity image
  # several consecutive pixels may be marked as belonging to an edge
  edges.areas <- detectVerticalEdgeAreas(img.gray)
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
  norm <- function(x) -((x - 1) / (10-1) - 1)  # 1--10 (smaller is better) --> 0--1 (larger is better)
  return(list(score=score.raw, gaussian=norm(score.gaussian), jpeg2k=norm(score.jpeg2k)))
}


#' Dominant color quality.
#' 
#' A measure of color quality in photos based on dominant colors.
#' The method uses a circular von Mises distribution to estimate the
#' dominant color and deviation on color hue dimension. It measures
#' the dominant colors and how concentrated the colors are within
#' the image.
#' 
#' The method was introduced in:
#' Sonia Ouni, Ezzeddine Zagrouba, Majed Chambah.
#' "A New No-reference Method for Color Image Quality Assessment."
#' \emph{International Journal of Computer Applications},
#' Vol. 40, No.17, February 2012.
#' URL: \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.259.1131&rep=rep1&type=pdf}
#' 
#' The article contains numerous errors: some of the examples don't match
#' with definitions and the exact methods are sometimes vague.
#' I've done my best to develop a method that is both practical and close
#' to the spirit of the authors' ideas.
#' Note that there are errors in the RGB-->HSV conversion formula the
#' article provides. See Wikipedia instead. 
#' 
#' This measure computes an average distance between pixels that have
#' dominant color. I still think the value should be normalized, as
#' in \code{custom.ds}, for example dividing by an average distance between
#' two pixels in the whole image, not just in dominant color. Otherwise
#' image size has an effect on the absolute magnitude of the number.
#' In the original article the values like 0.0020936 don't quite match
#' with the definitions.
#' @param img.hsv photo as an HSV image array (\emph{m} x \emph{n} x 3).
#' @return a list of estimated color model parameters:
#'   \itemize{
#'     \item \code{"mu"}: mean of von Mises distribution of color hues
#'       (in radians: \emph{0} is red, \emph{pi/2} is green,
#'       \emph{pi} is cyan, \emph{3*pi/2} is purple)
#'     \item \code{"kappa"}: compactness of von Mises distribution of
#'       color hues (a kind of "standard deviation"),
#'       with full hue circle as \emph{2*pi}. If hues are quite uniform,
#'       \code{kappa} may overshoot to >10.
#'     \item \code{"pi"}: proportion of pixels whose hue is within
#'       the dominant hue interval.
#'     \item \code{"ds"}: spatial dispersion of the dominant color.
#'     \item \code{"custom.ds"}: how much \code{ds} exceeds randomly
#'       distributed null hypothesis.
#'   }
#'   For an image that is (nearly) grayscale, the values are \code{NA} --
#'   no hue is defined for any pixel.
#' @examples
#' set.seed(42)
#' k <- 100*200
#' a <- array(runif(k), dim=c(100,200,3))  # all-gray image
#' all(is.na(dispersionDominantColor(toHSV(a))))
#' b <- array(c(0.3+0.6*runif(k), 0.55+0.4*runif(k), runif(k)), dim=c(100,200,3))  # lime-ish image
#' col.b <- dispersionDominantColor(toHSV(b))
#' all(abs(unlist(col.b) - c(1.844495, 1.048321, 0.57475, 81.17045, 1)) < 1e-4)
#' viewDDC(NULL, toHSV(b), col.b$mu, col.b$kappa)  # plot, source from "viewer.R"
#' 
#' @export

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
  
  # NPDC Euclidean distance approximation by random sample.
  px.1 <- smpl(1:npdc.nr)  # I think length(npdc.nr) == length(npdc.nc) so ok.
  px.2 <- smpl(1:npdc.nr)
  npdc.approx.dist <- mean(sqrt((row.idx[px.1] - row.idx[px.2])^2 + (col.idx[px.1] - col.idx[px.2])^2))

  spatial.dispersion <- npdc.approx.dist
  custom.ds <- min(1, npdc.approx.dist / approx.dist)  # cap to 1 (exceeds by random variation)
  
  return(list(mu=mu, kappa=kappa, pi=pi.measure, ds=spatial.dispersion, custom.ds=custom.ds))
}



# A set of basic measures, as described in slideshow names
# Automatic photo quality assessment --- Taming subjective problems with hand-coded metrics


#' Exposure quality measure (simple).
#' 
#' A simplistic measure of exposure quality in photos. The measurement is
#' based on the average luminosity on all pixels, mapped onto quality scores
#' by using a Beta distribution. Photos with middle-level exposures will
#' receive high scores (close to 1); photos that are on average overly
#' underexposed or overly overexposed will receive low scores (close to 0).
#' 
#' Inspiration from \emph{"Automatic photo quality assessment --
#' Taming subjective problems with hand-coded metrics"}:
#' 
#' \url{https://studylib.net/doc/17770925/automatic-photo-quality-assessment-taming-subjective-prob...}
#' 
#' Note that photos that contain only black and white may have medium-level
#' exposure on average.
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3).
#' @return numeric quality measure for photo exposure, ranging from 0 to 1 (best).
#' @examples
#' img <- array(sample(4:10, 40*80, replace=TRUE)/10, dim=c(40,80,3))  # high-exposure
#' quality <- basicExposureLevel(img)
#' abs(quality - 0.84) < 1e-2
#' scoref <- function(x) dbeta(x, 2, 2) / dbeta(0.5, 2, 2)
#' expo <- seq(0, 1, length.out=100)
#' plot(expo, scoref(expo), type='l', col='blue', xlab='exposure', ylab='quality')
#' abline(v=mean(luminance(img)), col='red', lty=3)
#' abline(h=quality, col='red', lty=3)
#' 
#' @export
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


#' RMS Contrast Measure.
#' 
#' Computes the \emph{RMS contrast} (Root Mean Square) for a photo. It is the
#' standard deviation of pixel intensities, where intensities are measured
#' by their luminance values.
#' 
#' Standard method, see for example
#' \url{https://en.wikipedia.org/wiki/Contrast_(vision)}
#' 
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3)
#'   with values from \code{[0;1]}.
#' @return numeric RMS contrast, ranging from 0 to 0.5 (most constrasty).
#' @examples
#' # high-exposure grayscale example
#' img <- array(sample(4:10, 40*80, replace=TRUE)/10, dim=c(40,80,3))
#' abs(basicRmsContrast(img) - 0.2) < 1e-3
#' 
#' # black and white maximum contrast example
#' img.bw <- array(sample(0:1, 40*80, replace=TRUE), dim=c(40,80,3))
#' abs(basicRmsContrast(img.bw) - 0.5) < 1e-3
#' 
#' @export
basicRmsContrast <- function(img) {
  # i <- luminance(img)
  # return(sqrt(sum((i - mean(i))^2) / (nrow(i)*ncol(i))))
  return(sd(luminance(img)))
}


#' Interval Contrast Measure.
#' 
#' Computes the \emph{Interval Contrast} measure for a photo. It is the
#' difference between lower and upper quantiles of the image's
#' luminance values --  it is related to dynamic range of the image.
#' 
#' The method is as follows.
#' Given an image and an interval size (proportion of pixels) \emph{s}:
#' \itemize{
#'   \item Compute luminance values \emph{L} for all pixels and
#'         sort them in increasing order.
#'   \item Find the lower quantile at \emph{(1 - s) / 2} on \emph{L}
#'   \item Find the upper quantile at \emph{1 - (1 - s) / 2} on \emph{L}
#'   \item Return the difference between upper and lower quantiles.
#' }
#' 
#' Inspiration from \emph{"Automatic photo quality assessment --
#' Taming subjective problems with hand-coded metrics"}:
#' 
#' \url{https://studylib.net/doc/17770925/automatic-photo-quality-assessment-taming-subjective-prob...}
#' 
#' @param img photo as an RGB image array (\emph{m} x \emph{n} x 3)
#'   with values from \code{[0;1]}.
#' @param interval interval size, as proportion of pixels of the image that
#'   belong to the interval. Default value is 0.95, which
#'   leaves 2.5\% of pixels to tails on both ends of luminance distribution.
#' @return numeric \emph{interval contrast}, ranging from 0 to 1 (most constrasty).
#' @examples
#' # Grayscale example of medium-level dynamic range
#' img <- array(runif(40*80, 0.4, 1.0), dim=c(40,80,3))
#' expected <- 0.95 * (1 - 0.4)  # 0.57
#' abs(basicIntervalContrast(img) - expected) < 1e-2
#' 
#' expected.half <- 0.5 * (1 - 0.4)  # 0.30
#' abs(basicIntervalContrast(img, 0.5) - expected.half) < 1e-2
#' 
#' # Range between maximum and minimum luminances
#' abs(basicIntervalContrast(img, interval=1.0) - 0.6) < 1e-2
#' 
#' # Maximum contrast example
#' img.max <- array(sample(0:1, 40*80, replace=TRUE), dim=c(40,80,3))
#' abs(basicIntervalContrast(img.max) - 1.0) < 1e-6
#' 
#' @export
basicIntervalContrast <- function(img, interval=0.95) {
  # Slides
  tail.prob <- (1 - interval) / 2  # two-sided distribution
  i <- sort(luminance(img))
  q <- quantile(i, c(tail.prob, 1 - tail.prob))
  return(as.numeric(diff(q)))
}

