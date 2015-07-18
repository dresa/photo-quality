/**
 * 
 */
package photo;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import javax.swing.*;
import java.awt.BorderLayout;
import javax.imageio.ImageIO;
import java.awt.image.DataBufferByte;

/**
 * Pixels in a photograph: 2-dimensional raster.
 * @author Esa Junttila
 */
public class Photo {
	private static final int NO_TRANSPARENCY = 255;

	private final int width;
	private final int height;
	private int[][] reds;
	private int[][] greens;
	private int[][] blues;
	private int[][] alphas;
	private int[][] argb;    // AlphaRBG pixels in 2-dim
	private int[] flatArgb;  // AlphaRBG pixels in 1-dim

	public Photo(String filename) throws IOException {
		this(new File(filename));
	}

	public Photo(File f) throws IOException {
		this(ImageIO.read(f));
	}

	public Photo(URL url) throws IOException {
		this(ImageIO.read(url));
	}

	public Photo(BufferedImage img) throws IOException {
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
				int argbColor = (a << 24) + (r << 16)  + (g << 8) + b;
				argb[row][col] = argbColor ;
				flatArgb[row*width+col] = argbColor; 
				alphas[row][col] = a;
				blues[row][col] = b;
				greens[row][col] = g;
				reds[row][col] = r;
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

	public void view() {
		viewImage(toImage());
	}

	/**
	 * Show an image within a frame.
	 * @param img image to be viewed
	 */
	private static void viewImage(BufferedImage img) {
		SwingUtilities.invokeLater(new Runnable() {  // FIXME: Convert into a Java 8 closure
			public void run() {
				JFrame photoFrame = new JFrame("Photo");
				photoFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
				ImageIcon imageIcon = new ImageIcon(img);
				JLabel jLabel = new JLabel();
				jLabel.setIcon(imageIcon);
				photoFrame.getContentPane().add(jLabel, BorderLayout.CENTER);
				photoFrame.pack();
				photoFrame.setLocationRelativeTo(null);
				photoFrame.setVisible(true);
		   }
		});
	}

	// Test method for development.
	public static void main(String[] args) throws IOException {
		Photo p = new Photo("examples/small_grid.png");
		p.view();
		p.savePNG("saved.png");
	}

}
