# photo-quality
Measure the aesthetic quality of photographs, based on cutting-edge methodology.

## Motivation

An avid photographer may easily end up having dozens of
seemingly identical photos. Or perhaps they wish to pick
a handful of best photos from their collection. They might hope
that they could give a set of photos to someone, who then
ranks the photos and picks the best ones.

This is possible with modern image analysis and machine
learning techniques. I'm talking about assessing the aesthetic
quality of photographs -- automatically. Unfortunately, none
of the applications from Adobe, Apple, or others seem to
address this issue. This project aims to fill the gap, that is
aiding photo management by cutting-edge methodology.

## Methodology

The question of photo quality is subjective to certain extent.
However, most photos are taken by ordinary people who follow
general rules as it comes to photo quality. In this domain,
we want to see excellent sharpness, little motion blur,
*golden rule* compositions, no compression artefacts (think of
JPEG), people's faces that have their eyes open, etc.

In fact many photo features, both technical and aesthetic in nature,
can be modeled and measured. Indeed, the field of measuring image
quality is a well-known research topic. Most of the research, however,
has concentrated on measuring quality when the true image is known.
But in photography, all we have is a single photo with no ground-truth.
This fairly new field of study is called
*no-reference image quality assessment*, or *NR* for short, and it
is in this field where we wish to make a contribution.

## Setup

We have explored the scientific literature on NR and measures of
image quality, and we have implemented approximately sixty such
measures. We wish to use these to predict how humans rate photos.

As test data, we use the *DPChallange* dataset, which consists
of 16504 photos, of both color and black & white. Based on
measurements on photos, we train a statistical model to predict
average ratings that a photo has received. Models like linear
regression and RandomForest were considered.

## Results

Concentrating on color images, we achieved a correlation
of 0.423 between predicted ratings and average ratings.
In total 60% of data was used for training a RandomForest
model, and the rest of the data was used to assess the
prediction accuracy of the model.

Each photo was rated by several people, so a *true* rating
is in fact a distribution of ratings, rather than average.
When we sampled ratings randomly from these distributions,
as we assume a human would do, we achieved a correlation
of 0.468 between sampled ratings and averages. Based on
these correlations, we claim,
**our method works at 90% level when compared to a human performance**.

![Rating predictions](https://i.imgur.com/xhA8Irj.png)
<!--- ![Rating predictions](http://www.elisanet.fi/esajakatja/valokuvat/pq/rating_predictions.png) -->

We used the stdev of ratings to measure the prediction error.
On an average photo in the test set, the differences between
individual ratings and the average have stdev 1.628.
As the absolute prediction error is only 0.615 (rating
difference from average) on the test set, we end up with
relative error of only 0.380 stdev's.
With 95% confidence, the predictions fall within one
stdev [-0.86;1.07] around average rating, as shown in the plot:

![Prediction error](https://i.imgur.com/5sIcaO0.png)
<!--- ![Prediction error](http://www.elisanet.fi/esajakatja/valokuvat/pq/prediction-rel-err-sd.png) -->


## Example: puffin photo (by andreas612)

**Original 755 x 501 JPG photo:**

* Original puffin photo, which I find quite pleasing.
![Original](https://i.imgur.com/gyWFfwc.jpg)
<!--- ![Original](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin.jpg) -->

* Here are the blur-related  measurements we have for the photo.
The blur measurements show that the photo has good quality. In practice the photos
always contain some amount of blur -- here the radius of blur is estimated to be 3 pixels.

| *Blur measures*                     | *Value*   | *Explanation* |
| :----------------------------------- | :------- | ------------- |
|   Blur annoyance quality (0--1):   | 0.84759 | the blur does not distract too much |
|   MDWE horizontal blur width:      | 3.29405 | the amount of blur is still modest |
|   MDWE Gaussian quality (0--1):    | 1       | blurry appearance is still appealing |
|   MDWE JPEG2000 quality (0--1):    | 1       | there are not JPEG2000 artefacts to speak of |

By measuring the blurriness properties we can arrive at some conclusions on
a photo's attractiveness. From a set of dozen similar photos, automatically picking the sharpest
one saves us the trouble of pixel-peeping the minor differences between candidates.

**Hues in the photo:**
* The image below shows the hues in the photo. The block-like structure probably originates from the JPEG compression.
In the HSV color space the variation within those blocks is fully conveyed by the Saturation and
Value dimensions. Note that colors that are almost white or black are quite insensitive to hue.
![Hues](https://i.imgur.com/YMngr1D.png)
<!--- ![Hues](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-hue.png) -->

* Let's have a look at the hue amounts in a polar coordinates.
The image shows that the most dominant hues are lime green and orange.
A best-fit von Mises distribution (black line) is used to identify dominating hues.

![Hue Histogram](https://i.imgur.com/xJW1ze6.png)
<!--- ![Hue Histogram](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-hue-histogram.png) -->


| *Hue Measurements* | *Value* | *Explanation* |
| :--------------- | :---- | :--- |
| Color dispersion(mu) | 0.6655 (rad), or 38.13&deg; | dominant hue is orange |
| Color dispersion(kappa) | 1.3797 | the dominance is not particularly concentrated |
| Color dispersion(pi) | 0.72621 | dominant hues cover 72.6 % of the photo's pixels |
| Color dispersion(ds) | 323.5 px | average distance between two pixels that have dominant hue |
| Color dispersion(custom.ds) | 0.970699 | pixels with dominant hues are almost as scattered as any two pixels |

Detecting dominant colors may help to decide whether a photo has incorrect color temperature
or other type of color cast. By measuring the scatter, we are able to distinguish
natural colors, such as blue sky and green forest, from an incorrect color temperature.

**Sharpness and blurriness:**
* In some sense sharpness and blurriness are the opposites of each other.
A photograph may very well contain both, and still be of a high quality.
* In the sharpness image the sharpest pixels (white) are those at the strong edges:
![Sharpness](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-sharpness.png)
<!--- ![Sharpness](https://i.imgur.com/NWIbF4a.png) -->
* In the blurriness image the background of the photo appears as highly blurry (white):
![Blurriness](https://i.imgur.com/UtrSK2A.png)
<!--- ![Blurriness](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-blurriness.png) -->

**Image segmentation:**
* Image segmentation offers a way to split the image so that objects
and coherent areas can be detected. In many photos, we have some areas filled with more
or less the same color, while the areas may differ from each other. By detecting the segments, we
simultaneously split the image into distinguishable parts.
* In the following image we have decided that the puffin photo is naturally split into 26 clusters,
each of which corresponds to a single color. Using these clustered colors, the image has 20212 connected
color segments, especially below in the grass area. Here they are, each shown with a random
color to make them stand out:
![Clustered segments](https://i.imgur.com/TtNq9lw.png)
<!--- ![Clustered segments](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-connected-segments-20212-from-26-clusters.png) -->
* In the connected segments image, we show ten largest connected components against a black background:
![Largest segments](https://i.imgur.com/GtDj4sC.png)
<!--- ![Largest segments](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-largest-connected-components-10.png) -->
* In the reconstructed image we have used only 26 clusters, that is colors,
to reconstruct the original photo:
![Reconstructed from clusters](https://i.imgur.com/6YDwEeI.png)
<!--- ![Reconstructed from clusters](http://www.elisanet.fi/esajakatja/valokuvat/pq/puffin-reconstructed-26-clusters.png) -->

**Exposure and contrast:**

Three simple measures for detecting  whether image has pleasing exposure and contrast:

| *Measurement* | *Value* | *Explanation* |
| :------------ | :------ | :------------ |
| Basic exposure | 0.933 | Good exposure, near middle-grey |
| Basic RMS contrast | 0.136 | Natural luminance spread (dark vs. light, near 0.18) |
| Basic interval contrast | 0.570 | moderate dynamics: middle 95% of photo has 57% of full dynamic range |

**More measures:**

In Datta et. al (2006), the authors propose dozens of measures that are
potentially related to photo quality. We have implemented most of them,
and for the puffin image we obtained the following values.

Datta R., Joshi D., Li J. and Wang J.:
*Studying Aesthetics in Photographic Images Using a Computational Approach*, 2006.

| *Measurement* | *Datta variant* | *Value* | *Explanation* |
| :------------ | :-------- | :------ | :------------ |
| Average Intensity | 1 | 0.394 | Slightly darker than neutral grey |
| Colorfulness | 2 | 61.03 | Ordinary distance to equally distributed colors |
| Colorfulness-Grey | 2 (grey, n=6) | 25.14 |  Quite colorful (distance from a greyscale image) |
| Average saturation | 3 | 0.244 |  |
| Average hue | 4 | 103.02 | Average hue is green |
| Average central hue | 5 | 91.55 | In the image center we have green with a hint of yellow |
| Average central saturation | 6 | 0.196 | Center is less saturated than the rest of photo |
| Average central intensity | 7 | 0.475 | Center is brigher than edge areas |
| Texture, hue | 10,11,12;19 | 0.00803, 0.0259, 0.0921; sum 0.126 |  |
| Texture, saturation | 13,14,15;20 | 0.0125, 0.0316, 0.0747; sum 0.119 |  |
| Texture, value | 16,17,18;21 | 0.00955, 0.0248, 0.0584; sum 0.0927 |  |
| Size feature | 22 | 1256 | Modest-sized photo: sum of rows and columns |
| Aspect ratio | 23 | 0.664, 1.51 | Classical 3:2 aspect ratio |
| Number of large patches | 24 | 5 | Full number of large segments (max 5) |
| Number of clusters | 25 | 26 | Photo is complex: best described by 26 clusters (random variation) |
| Average patch HSV values | 26--40 | 272.6, 63.2, 30.1, 25.3, 224.4; 0.0586, 0.0867, 0.194, 0.137, 0.0438; 0.238, 0.276, 0.382, 0.329, 0.376 |  |
| Relative patch sizes | 41--45 | 0.112, 0.0616, 0.0604, 0.0454, 0.0389 | None of the segments fills the photo fully |
| Segment position codes | 48--52 | 12, 12, 21, 12, 21 |  Five largest segment centers reside on top, left-middle part of photo |
| Segment distances from center | Esa proxy 48--52 | 0.492, 0.446, 0.443, 0.324, 0.353 | Five largest segments not in the center nor corners |
| Only center within DOF for H, S, V | 53--55 rescaled | 0.113, -0.109, 0.109 | |


## Future roadmap

The project leans towards the scientific research on no-reference image quality. It combines the scattered results from numerous articles and provides an implementation for the methods that lack one. Together, a number of photographic features can be measured from a photo. To assess the quality of a photo, we use machine learning techniques on photo features to arrive at a fair quality assessment. As training data, we use publicly available photo collections that contain individual photo ratings.

1. Find relevant articles and measures for no-reference quality assessment.
2. Implement the most relevant methods.
3. Download a collection of photos with ratings from a public source.
4. Train a machine learning model to map measurements into a photo rating.
5. Learn and improve on what we have.


This software has been written in R version 3.4.3. It serves as a
prototype for photo quality assessment, and as a library of working implementations.
Keeping a future production version in mind, the code contains a lot of details
that support implementation work in lower-level languages.

In current code, we use the following R packages along with their dependencies:
emdist, waveslim, NbClust, jpeg, png, ggplot2, circular and mmand.
