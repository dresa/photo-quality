/**
 * 
 */
package photo;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;

import javax.swing.*;

import util.Tools;
import viewer.PhotoViewer;

import java.awt.BorderLayout;
import javax.imageio.ImageIO;
import java.awt.image.DataBufferByte;

/**
 * Pixels in a photograph: 2-dimensional raster.
 * @author Esa Junttila
 */
public class PhotoBasic implements Photo {
	private static final int NO_TRANSPARENCY = 255;

	private final int width;
	private final int height;
	private int[][] reds;
	private int[][] greens;
	private int[][] blues;
	private int[][] alphas;
	private int[][] argb;    // AlphaRBG pixels in 2-dim
	private int[] flatArgb;  // AlphaRBG pixels in 1-dim

	public PhotoBasic(String filename) throws IOException {
		this(new File(filename));
	}

	public PhotoBasic(File f) throws IOException {
		this(ImageIO.read(f));
	}

	public PhotoBasic(URL url) throws IOException {
		this(ImageIO.read(url));
	}

	public PhotoBasic(BufferedImage img) throws IOException {
		int type = img.getType(); 
		if (type != BufferedImage.TYPE_3BYTE_BGR) 
			throw new RuntimeException("Unsupported image color type: " + type);
		this.width = img.getWidth();
		this.height = img.getHeight();
		boolean alphaExists = img.getAlphaRaster() != null;

		// Data in 1-dim array form, with separate channels.
		// For example size could be =width*height*3.
		final byte[] pixels = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();

		argb = new int[height][width];  // AlphaRBG pixels in 2-dim
		flatArgb = new int[height*width];  // AlphaRBG pixels in 1-dim
		alphas = new int[height][width];
		reds = new int[height][width];
		greens = new int[height][width];
		blues = new int[height][width];

		int p = 0;
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				// Signed bytes to unsigned ints: (for example -1 translates to 255).
				// Assumes (A)BGR ordering.
				int a = (alphaExists ? pixels[p++] & 0xff : NO_TRANSPARENCY);  // alpha
				int b = pixels[p++] & 0xff;  // blue 
				int g = pixels[p++] & 0xff;  // green
				int r = pixels[p++] & 0xff;  // red
				int argbColor = Pixel.toARGB(a, r, g, b);
				argb[row][col] = argbColor ;
				flatArgb[row*width+col] = argbColor; 
				alphas[row][col] = a;
				blues[row][col] = b;
				greens[row][col] = g;
				reds[row][col] = r;
			}
		}
	}

	public PhotoBasic(int[][] a, int[][] r, int[][] g, int[][] b) throws IOException {
		if (r.length != g.length || r.length != b.length || g.length != b.length) {
			throw new IllegalArgumentException("Mismatch in RGB arg dimensions");
		}
		this.height = r.length;
		this.width = r[0].length; 
		this.alphas = Tools.copy2DArray(a);
		this.reds = Tools.copy2DArray(r);
		this.greens = Tools.copy2DArray(g);
		this.blues = Tools.copy2DArray(b);
		initDerivedFields(alphas, reds, greens, blues);
	}

	private void initDerivedFields(int[][] a, int[][] r, int[][] g, int[][] b) {
		final int numRows = r.length;
		final int numCols = r[0].length;
		this.argb = new int[numRows][numCols];
		this.flatArgb = new int[numRows*numCols];
		for (int row = 0; row < numRows; row++) {
			for (int col = 0; col< numCols; col++) {
				int px = Pixel.toARGB(a[row][col], r[row][col], g[row][col], b[row][col]);
				argb[row][col] = px;
				flatArgb[row*numCols+col] = px;
			}
		}
	}

	public BufferedImage toImage() {
		BufferedImage img = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
		img.setRGB(0, 0, width, height, flatArgb, 0, width);
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
