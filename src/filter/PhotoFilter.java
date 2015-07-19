package filter;

import photo.Channel;
import photo.Photo;

public interface PhotoFilter {
	Channel filter(Photo p);
}
