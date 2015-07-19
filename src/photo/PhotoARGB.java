package photo;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.File;
import java.io.IOException;
import java.net.URL;

import javax.imageio.ImageIO;

import util.Tools;

public class PhotoARGB implements Photo {
	public enum ChannelType { ALPHA, RED, GREEN, BLUE }

	private static final int NO_TRANSPARENCY = 255;

	private final int width;
	private final int height;
	private ChannelByte reds;
	private ChannelByte greens;
	private ChannelByte blues;
	private ChannelByte alphas;

	public PhotoARGB(String filename) throws IOException {
		this(new File(filename));
	}

	public PhotoARGB(File f) throws IOException {
		this(ImageIO.read(f));
	}

	public PhotoARGB(URL url) throws IOException {
		this(ImageIO.read(url));
	}

	public PhotoARGB(BufferedImage img) throws IOException {
		int type = img.getType(); 
		if (type != BufferedImage.TYPE_3BYTE_BGR) 
			throw new RuntimeException("Unsupported image color type: " + type);
		this.width = img.getWidth();
		this.height = img.getHeight();
		boolean alphaExists = img.getAlphaRaster() != null;

		// Data in 1-dim array form, with separate channels.
		// For example size could be =width*height*3.
		final byte[] pixels = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();

		int[][] argbArray = new int[height][width];  // AlphaRBG pixels in 2-dim
		byte[][] alphaArray = new byte[height][width];
		byte[][] redArray = new byte[height][width];
		byte[][] greenArray = new byte[height][width];
		byte[][] blueArray = new byte[height][width];

		int p = 0;
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				// Signed bytes to unsigned ints: (for example -1 translates to 255).
				// Assumes (A)BGR ordering.
				byte a = (byte) (alphaExists ? pixels[p++] & Tools.MASK_8_BIT : NO_TRANSPARENCY);  // alpha
				byte b = (byte) (pixels[p++] & Tools.MASK_8_BIT);  // blue 
				byte g = (byte) (pixels[p++] & Tools.MASK_8_BIT);  // green
				byte r = (byte) (pixels[p++] & Tools.MASK_8_BIT);  // red
				int argbColor = Pixel.toARGB(a, r, g, b);
				argbArray[row][col] = argbColor;
				alphaArray[row][col] = a;
				blueArray[row][col] = b;
				greenArray[row][col] = g;
				redArray[row][col] = r;
			}
		}
		this.alphas = new ChannelByte(alphaArray);
		this.reds = new ChannelByte(redArray);
		this.greens = new ChannelByte(greenArray);
		this.blues = new ChannelByte(blueArray);
	}

	public PhotoARGB(ChannelByte alpha, ChannelByte red, ChannelByte green, ChannelByte blue) throws IOException {
		Dim aDim = alpha.getDimensions();
		Dim rDim = red.getDimensions();
		Dim gDim = green.getDimensions();
		Dim bDim = blue.getDimensions();
		if (aDim != rDim || aDim != gDim || aDim != bDim) {
			throw new IllegalArgumentException(
				"Mismatch in AlphaRGB dimensions: " +
				"A=" + aDim + ", R=" + rDim + ", G=" + gDim + ", B=" + bDim);
		}
		this.height = aDim.getY();
		this.width = aDim.getX();
		this.alphas = alpha;
		this.reds = red;
		this.greens = green;
		this.blues = blue;
	}

	public ChannelByte getChannel(ChannelType type) {
		switch (type) {
			case ALPHA: return alphas;
			case RED: return reds;
			case GREEN: return greens;
			case BLUE: return blues;
			default:
				throw new IllegalArgumentException("Invalid ARGB channel: " + type);
		}
	}

	public BufferedImage toImage() {
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int[][] argb = Tools.deriveARGB(alphas.getBytes(), reds.getBytes(), greens.getBytes(), blues.getBytes());
		img.setRGB(0, 0, width, height, Tools.toArray1D(argb), 0, width);
		return img;
	}

	public String toString() {
		return
			"PhotoARGB(w="+width+", h="+height+"):\n" +
			"  alpha=" + alphas + "\n" +
			"    red=" + reds + "\n" +
			"  green=" + greens + "\n" +
			"   blue=" + blues + "\n" +
			")";
	}
}
