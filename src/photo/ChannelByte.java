package photo;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import util.Tools;

public class ChannelByte implements Channel {

	final int width;
	final int height;
	private byte[][] vals;

	public ChannelByte(byte[][] values) {
		this.height = values.length;
		this.width = values[0].length;
		this.vals = Tools.copy2DArray(values);		
	}

	public ChannelByte(byte[] values, int width) {
		final int n = values.length;
		this.width = width;
		this.height = n / width;
		if (n % width != 0) throw new IllegalArgumentException(
			"Number of channel values " + n + " not divisible by columns " + width);
		this.vals = new byte[height][width];
		for (int i = 0; i < n; i++) {
			this.vals[i/width][i%width] = values[i];  
		}
	}
	public byte[][] getBytes() {
		return Tools.copy2DArray(this.vals);
	}
	@Override
	public float[][] getFloats() {
		return Tools.toFloatArray2D(this.vals);
	}
	@Override
	public Dim getDimensions() { return new Dim(width, height); }

	@Override
	public String toString() {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		Tools.write2DArray(vals, new PrintStream(baos));
		return baos.toString();
		/*
		StringBuilder sb = new StringBuilder();
		Formatter form = new Formatter(sb, Locale.US);
		sb.append("ChannelBytes(w="+width+", h="+height+"):\n");
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				byte v = vals[row][col];
				form.format("%4d", v).toString();
				if ( !(col == width-1 && row == height-1) ) sb.append(", ");
				if (col == width-1) sb.append('\n');
			}
		}
		form.close();
		return sb.toString();
		*/
	}
}
