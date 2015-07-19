package filter;

import photo.ChannelFloat;
import photo.Photo;
import photo.PhotoARGB;
import photo.PhotoARGB.ChannelType;

// FIXME: add a scientific definition of luminosity, including inverse gamma etc. 
public class Luminosity implements PhotoFilter {
	public ChannelFloat filter(Photo p) {
		if ( !(p instanceof PhotoARGB)) throw new UnsupportedOperationException(
			"Luminosity for other than ARGB photos");
		PhotoARGB photoRGB = (PhotoARGB) p;
		return ChannelMath.luminosity(
			photoRGB.getChannel(ChannelType.RED),
			photoRGB.getChannel(ChannelType.GREEN),
			photoRGB.getChannel(ChannelType.BLUE)
		);
	}
}

