package photo;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;

import javax.imageio.ImageIO;

import filter.Luminosity;
import util.Tools;

public class PhotoGreyscale implements Photo {
	public enum ChannelType { GREY }

	private final int width;
	private final int height;
	private ChannelByte greys;

	public PhotoGreyscale(String filename) throws IOException {
		this(new File(filename));
	}

	public PhotoGreyscale(File f) throws IOException {
		this(ImageIO.read(f));
	}

	public PhotoGreyscale(URL url) throws IOException {
		this(ImageIO.read(url));
	}

	public PhotoGreyscale(BufferedImage img) throws IOException {
		this(new Luminosity().filter(new PhotoARGB(img)));
	}
	public PhotoGreyscale(ChannelFloat greys) throws IOException {
		this(greys.toBytes());
	}
	public PhotoGreyscale(ChannelByte greys) throws IOException {
		Dim gDim = greys.getDimensions();
		this.height = gDim.getY();
		this.width = gDim.getX();
		this.greys = greys;
	}

	public ChannelByte getChannel(ChannelType type) {
		switch (type) {
			case GREY: return greys;
			default:
				throw new IllegalArgumentException("Invalid greyscale channel: " + type);
		}
	}

	public BufferedImage toImage() {
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		byte[][] vals = greys.getBytes();
		byte[][] alphas = new byte[height][width];
		for (int row = 0; row < height; row++) {
			Arrays.fill(alphas[row], (byte) Tools.NO_TRANSPARENCY);
		}
		int[][] argb = Tools.deriveARGB(alphas, vals, vals, vals);
		img.setRGB(0, 0, width, height, Tools.toArray1D(argb), 0, width);
		return img;
	}
	public String toString() {
		return
			"PhotoGreyscale(w="+width+", h="+height+"):\n" +
			"  greys=" + greys + "\n" +
			")";
	}
}
