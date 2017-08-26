# Download a DPChallenge photo dataset from the web,
# as instructed in http://ritendra.weebly.com/aesthetics-datasets.html
#
# Note that some photos have been removed from the web server, for example 65205.
#
# Esa Junttila, 2017-08-26

INPUT_FILE <- 'dpchallenge_dataset.txt'  # with added headers
OUTPUT_DIR <- 'originals/'
URL_BASE <- 'http://www.dpchallenge.com/image.php?IMAGE_ID='
JPG_REGEXP <- 'http:[^:]*\\.jpg'

scrape <- function(photo) {
  page.url <- paste0(URL_BASE, toString(photo))
  rawHTML <- paste(readLines(page.url), collapse="\n")
  matches <- gregexpr(JPG_REGEXP, rawHTML, perl=TRUE, ignore.case=TRUE)
  all.img.urls <- regmatches(rawHTML, matches)[[1]]
  photo.suffix <- paste0(toString(photo), '\\.jpg')
  filter.mask <- grep(photo.suffix, all.img.urls, perl=TRUE, ignore.case=TRUE)
  img.url <- all.img.urls[filter.mask]
  target.file <- paste0(OUTPUT_DIR, photo, '.jpg')
  res <- tryCatch(
    {download.file(img.url, target.file, mode='wb', quiet=TRUE)
     print(paste('Downloaded', target.file))},
    error = function(e) {print(paste('Failed', target.file))}
  )
}

## Main
df <- read.table(INPUT_FILE, header=TRUE, stringsAsFactors=FALSE)
photos <- df$Photo  # second column in the file
sapply(photos, scrape)
