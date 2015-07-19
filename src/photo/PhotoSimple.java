package photo;

import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.File;
import java.io.IOException;
import java.net.URL;

import javax.imageio.ImageIO;

import util.Tools;
import viewer.PhotoViewer;

public class PhotoSimple {
	private static final int NO_TRANSPARENCY = 255;
	private static final int MASK_8_BIT = 0xff;

	private final int width;
	private final int height;
	private ChannelByte reds;
	private ChannelByte greens;
	private ChannelByte blues;
	private ChannelByte alphas;

	public PhotoSimple(String filename) throws IOException {
		this(new File(filename));
	}

	public PhotoSimple(File f) throws IOException {
		this(ImageIO.read(f));
	}

	public PhotoSimple(URL url) throws IOException {
		this(ImageIO.read(url));
	}

	public PhotoSimple(BufferedImage img) throws IOException {
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
				byte a = (byte) (alphaExists ? pixels[p++] & MASK_8_BIT : NO_TRANSPARENCY);  // alpha
				byte b = (byte) (pixels[p++] & MASK_8_BIT);  // blue 
				byte g = (byte) (pixels[p++] & MASK_8_BIT);  // green
				byte r = (byte) (pixels[p++] & MASK_8_BIT);  // red
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

	public PhotoSimple(ChannelByte alpha, ChannelByte red, ChannelByte green, ChannelByte blue) throws IOException {
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

	public BufferedImage toImage() {
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		int[][] argb = Tools.deriveARGB(alphas.getValues(), reds.getValues(), greens.getValues(), blues.getValues());
		img.setRGB(0, 0, width, height, Tools.toArray1D(argb), 0, width);
		return img;
	}

	public void savePNG(String filename) throws IOException {
	    ImageIO.write(toImage(), "png", new File(filename));		
	}

	public void saveJPG(String filename) throws IOException {
		// FIXME: there is a bug with false colors in JPG file
		// that shows correctly on the screen
	    ImageIO.write(toImage(), "jpg", new File(filename));		
	}

	// Test method for development.
	public static void main(String[] args) throws IOException {
		if (args.length == 1) {
			String filename = args[0];
			PhotoBasic p = new PhotoBasic(filename);
			PhotoViewer.view(p);
			p.savePNG("saved.png");
			//p.saveJPG("saved.jpg");
		}
		else {
			System.err.println("Usage: java PhotoBasic <filename>");
		}
	}
}
