package filter;

import photo.Channel;
import photo.ChannelByte;
import photo.ChannelFloat;
import photo.Dim;
import util.Tools;

public class ChannelMath {
	public static ChannelFloat linear(Channel ch, float add, float mult) {
		float[][] f = ch.getFloats();
		for (int r = 0; r < f.length; r++) {
			for (int c = 0; c < f[0].length; c++) {
				f[r][c] = add + mult*f[r][c]; 
			}
		}
		return new ChannelFloat(f);
	}
	public static ChannelFloat add(Channel ch, float add) {
		float[][] f = ch.getFloats();
		for (int r = 0; r < f.length; r++) {
			for (int c = 0; c < f[0].length; c++) {
				f[r][c] = add + f[r][c]; 
			}
		}
		return new ChannelFloat(f);
	}
	public static ChannelFloat mult(Channel ch, float mult) {
		float[][] f = ch.getFloats();
		for (int r = 0; r < f.length; r++) {
			for (int c = 0; c < f[0].length; c++) {
				f[r][c] = mult*f[r][c]; 
			}
		}
		return new ChannelFloat(f);
	}
	public static ChannelFloat weighted(Channel[] chs, float[] weights) {
		if (chs.length != weights.length) throw new IllegalArgumentException(
			"Mismatch: num channels and num weights");
		Dim dim = chs[0].getDimensions();
		float[][] target = new float[dim.getY()][dim.getX()];
		for (int i = 0; i < chs.length; i++) {
			float[][] f = chs[i].getFloats();
			float w = weights[i];
			for (int r = 0; r < f.length; r++) {
				for (int c = 0; c < f[0].length; c++) {
					target[r][c] += w*f[r][c]; 
				}
			}
		}
		return new ChannelFloat(target);
	}
	public static ChannelFloat luminosity(ChannelByte red, ChannelByte green, ChannelByte blue) {
		Dim dim = red.getDimensions();
		final int width = dim.getX();
		final int height = dim.getY();
		byte[][] reds = red.getBytes();
		byte[][] greens = green.getBytes();
		byte[][] blues = blue.getBytes();
		float[][] target = new float[height][width];
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				final int r = reds[row][col] & Tools.MASK_8_BIT;
				final int g = greens[row][col] & Tools.MASK_8_BIT;
				final int b = blues[row][col] & Tools.MASK_8_BIT;
				target[row][col] = (float) Math.sqrt(0.299*r*r + 0.587*g*g + 0.114*b*b);				
			}
		}
		return new ChannelFloat(target);
	}
}

