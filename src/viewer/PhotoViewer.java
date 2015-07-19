package viewer;

import java.awt.BorderLayout;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import filter.Luminosity;
import photo.Photo;
import photo.PhotoARGB;
import photo.PhotoGreyscale;

public class PhotoViewer {
	public static void view(Photo photo) {
		view(photo.toImage());
	}

	/**
	 * Show an image within a frame.
	 * @param img image to be viewed
	 */
	public static void view(BufferedImage img) {
		SwingUtilities.invokeLater(new Runnable() {  // FIXME: Convert into a Java 8 closure
			public void run() {
				JFrame photoFrame = new JFrame("Photo");
				photoFrame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
				ImageIcon imageIcon = new ImageIcon(img);
				JLabel jLabel = new JLabel();
				jLabel.setIcon(imageIcon);
				JScrollPane jsp = new JScrollPane(jLabel);
				photoFrame.getContentPane().add(jsp, BorderLayout.CENTER);
				photoFrame.pack();
				photoFrame.setLocationRelativeTo(null);
				photoFrame.setVisible(true);
		   }
		});
	}

	public static void savePNG(Photo p, String filename) throws IOException {
	    ImageIO.write(p.toImage(), "png", new File(filename));		
	}

	public static void saveJPG(Photo p, String filename) throws IOException {
		// FIXME: there is a bug with false colors in JPG file
		// that shows correctly on the screen
	    ImageIO.write(p.toImage(), "jpg", new File(filename));		
	}

	// Test method for development.
	public static void main(String[] args) throws IOException {
		if (args.length == 1) {
			String filename = args[0];
			PhotoARGB p = new PhotoARGB(filename);
			System.out.print(p);
			Photo luminosityPhoto = new PhotoGreyscale(new Luminosity().filter(p)); 
			System.out.print(luminosityPhoto);
			PhotoViewer.view(p);
			//PhotoViewer.view(luminosityPhoto);
			savePNG(p, "saved.png");
			savePNG(luminosityPhoto, "saved_lum.png");
			//saveJPG(p, "saved.jpg");
		}
		else {
			System.err.println("Usage: java PhotoBasic <filename>");
		}
	}
}
