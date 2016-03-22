# photo-quality
Measure the quality of photographs, based on cutting-edge methodology.

## Example: penguin photo (by andreas612)

**Original 755 x 501 JPG photo:**

* Original penguin photo, which I find quite pleasing.
![Original](http://www.elisanet.fi/esajakatja/valokuvat/pq/penguin.jpg)

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
![Hues](http://www.elisanet.fi/esajakatja/valokuvat/pq/penguin-hue.png)

* Let's look at the hue amounts in a polar coordinates.
The image shows that the most dominant hues are lime green and orange.
A best-fit von Mises distribution (black line) is used to identify dominating hues.

![Hue Histogram](http://www.elisanet.fi/esajakatja/valokuvat/pq/penguin-hue-histogram.png)


| *Hue Measurements* | *Value* | *Explanation* |
| :--------------- | :---- | :--- |
| Color dispersion(mu) | 0.6655 (rad), or 38.13&deg; | dominant hue is orange |
| Color dispersion(kappa) | 1.3797 | the dominance is not particularly concentrated |
| Color dispersion(pi) | 0.72621 | dominant hues cover 72.6 % of the photo's pixels |
| Color dispersion(ds) | 323.5 px | average distance between two pixels that have dominant hue |
| Color dispersion(custom.ds) | 0.970699 | pixels with dominant hues are almost as scattered as any two pixels |

Detecting dominant colors may help to decide whether a photo has incorrect color temperature
or other type of color cast. By measuring the scatter, we are able to distinguish natural colors, such as blue sky and green forest,
from an incorrect color temperature.

**Sharpness and blurriness:**
* In some sense sharpness and blurriness are the opposites of each other. A photograph may very well contain both, ands still be of a high quality.
* In the sharpness image the sharpest pixels (white) are those at the strong edges:
![Sharpness](http://www.elisanet.fi/esajakatja/valokuvat/pq/penguin-sharpness.png)
* In the blurriness image the background of the photo appears as highly blurry (white):
![Blurriness](http://www.elisanet.fi/esajakatja/valokuvat/pq/penguin-blurriness.png)

##Why?
It's tedious to manage a large collection of photographs. An avid photographer may easily end up having dozens of seemingly identical photos. Unfortunately, none of the applications from Adobe, Apple, or others seem to address this issue. It's time to ease the pain of photo management and spend the time taking photos instead.

Imagine that you could give a set of photos to someone, who then orders them and picks the best photos for you. This is possible with modern image analysis and machine learning techniques. I'm talking about assessing the quality of photographs -- automatically.

The question of photo quality is subjective to certain extent. However, most photos are taken by ordinary people, and more general rules apply to photo quality. In this domain, most people would favor photos that have excellent sharpness somewhere, avoid motion blur, lack compression artefacts (think of JPEG), have faces with their eyes open, follow the *golden rule*, etc.

In fact many photo features, both technical and aesthetic in nature, can be modeled and measured. The field of measuring image quality is a well-known research topic. Most of the research, however, has concentrated on measuring quality when the true image is known. But in photography, all we have is a single photo. This fairly new field of study is called *no-reference image quality assessment*, or *NR* for short.

## Future roadmap

The project leans towards the scientific research on no-reference image quality. It combines the scattered results from numerous articles and provides an implementation for the methods that lack one. Together, a number of photographic features can be measures from a photo. To assess the quality of a photo, we use machine learning techniques on photo features to arrive at a fair quality assessment. As training data, we use publicly available photo collections that contain individual photo ratings.

1. Find relevant articles and measures for no-reference quality assessment.
2. Implement the most relevant methods.
3. Download a collection of photos with ratings from a public source.
4. Train a machine learning model to map measurements into a photo rating.
5. Learn and improve on what we have.


This software has been written in R version 3.2.1. It serves as a
prototype for photo quality assessment, and as a library of working implementations.
Keeping a future production version in mind, the code contains a lot of details
that support implementation work in lower-level languages.

