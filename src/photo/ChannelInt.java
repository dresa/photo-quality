package photo;

import util.Tools;

public class ChannelInt implements Channel {

	final int width;
	final int height;
	private int[][] vals;

	public ChannelInt(int[][] values) {
		this.height = values.length;
		this.width = values[0].length;
		this.vals = Tools.copy2DArray(values);		
	}

	public ChannelInt(int[] values, int width) {
		final int n = values.length;
		this.width = width;
		this.height = n / width;
		if (n % width != 0) throw new IllegalArgumentException(
			"Number of channel values " + n + " not divisible by columns " + width);
		this.vals = new int[height][width];
		for (int i = 0; i < n; i++) {
			this.vals[i/width][i%width] = values[i];  
		}
	}
	public int[][] getValues() {
		return Tools.copy2DArray(this.vals);
	}
	@Override
	public Dim getDimensions() { return new Dim(width, height); }
}
