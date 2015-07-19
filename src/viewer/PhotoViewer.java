package viewer;

import java.awt.BorderLayout;
import java.awt.image.BufferedImage;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.WindowConstants;

import photo.Photo;

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
}
