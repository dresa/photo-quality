package photo;

import util.Tools;

public class ChannelFloat implements Channel {

	final int width;
	final int height;
	private float[][] vals;

	public ChannelFloat(float[][] values) {
		this.height = values.length;
		this.width = values[0].length;
		this.vals = Tools.copy2DArray(values);		
	}

	public ChannelFloat(float[] values, int width) {
		final int n = values.length;
		this.width = width;
		this.height = n / width;
		if (n % width != 0) throw new IllegalArgumentException(
			"Number of channel values " + n + " not divisible by columns " + width);
		this.vals = new float[height][width];
		for (int i = 0; i < n; i++) {
			this.vals[i/width][i%width] = values[i];
		}
	}
	public ChannelByte toBytes() {
		return new ChannelByte(Tools.roundToByte(vals));
	}
	public float[][] getFloats() {
		return Tools.copy2DArray(this.vals);
	}
	@Override
	public Dim getDimensions() { return new Dim(width, height); }
}
