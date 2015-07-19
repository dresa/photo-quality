package util;

import java.io.PrintWriter;

import photo.Pixel;

public final class Tools {
	public static byte[] toArray1D(byte[][] src) {
		final int height = src.length;
		final int width = src[0].length;
		byte[] target = new byte[height*width];
		for (int row = 0; row < height; row++) {
			System.arraycopy(src[row], 0, target, row*width, width);
		}
		return target;
	}
	public static int[] toArray1D(int[][] src) {
		final int height = src.length;
		final int width = src[0].length;
		int[] target = new int[height*width];
		for (int row = 0; row < height; row++) {
			System.arraycopy(src[row], 0, target, row*width, width);
		}
		return target;
	}
	public static byte[][] toArray2D(byte[] src, int width) {
		final int n = src.length;
		final int height = n/width;
		if (n%width != 0) throw new IllegalArgumentException(
			"Number of values " + n + " not divisible by width " + width);
		byte[][] target = new byte[height][width];
		for (int row = 0; row < height; row++) {
			System.arraycopy(src, row*width, target[row], 0, width);
		}
		return target;
	}
	public static int[][] toArray2D(int[] src, int width) {
		final int n = src.length;
		final int height = n/width;
		if (n%width != 0) throw new IllegalArgumentException(
			"Number of values " + n + " not divisible by width " + width);
		int[][] target = new int[height][width];
		for (int row = 0; row < height; row++) {
			System.arraycopy(src, row*width, target[row], 0, width);
		}
		return target;
	}
	public static byte[] copy1DArray(byte[] src) {
	    int n = src.length;
	    byte[] target = new byte[n];
        System.arraycopy(src, 0, target, 0, n);
	    return target;
	}
	public static byte[][] copy2DArray(byte[][] src) {
	    int length = src.length;
	    byte[][] target = new byte[length][src[0].length];
	    for (int i = 0; i < length; i++) {
	        System.arraycopy(src[i], 0, target[i], 0, src[i].length);
	    }
	    return target;
	}
	public static int[][] copy2DArray(int[][] src) {
	    int length = src.length;
	    int[][] target = new int[length][src[0].length];
	    for (int i = 0; i < length; i++) {
	        System.arraycopy(src[i], 0, target[i], 0, src[i].length);
	    }
	    return target;
	}
	public static void write2DArray(int[][] array, PrintWriter pw) {
		int height = array.length;
		for (int r=0; r<height; r++) {
			int width = array[r].length;
			for (int c=0; c<width; c++) {
				pw.print(array[r][c] + " ");
			}
			pw.println();
		}
	}

	public static int[][] deriveARGB(byte[][] a, byte[][] r, byte[][] g, byte[][] b) {
		final int numRows = r.length;
		final int numCols = r[0].length;
		int[][] argb = new int[numRows][numCols];
		for (int row = 0; row < numRows; row++) {
			for (int col = 0; col< numCols; col++) {
				int px = Pixel.toARGB(a[row][col], r[row][col], g[row][col], b[row][col]);
				argb[row][col] = px;
			}
		}
		return argb;
	}
}
