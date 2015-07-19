package util;

import java.io.PrintStream;
import java.util.Formatter;
import java.util.Locale;

import photo.Pixel;

public final class Tools {
	public static final int MASK_8_BIT = 0xff;  // FIXME: move somewhere
	public static final byte NO_TRANSPARENCY = (byte) 255;  // = -1

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
	public static float[][] toFloatArray2D(byte[][] src) {
	    int numRows = src.length;
	    int numCols = src[0].length;
	    float[][] target = new float[numRows][numCols];
	    for (int row = 0; row < numRows; row++) {
		    for (int col = 0; col < numCols; col++) {
		    	target[row][col] = src[row][col] & Tools.MASK_8_BIT;
		    }
	    }
	    return target;
	}
	public static float[][] toFloatArray2D(int[][] src) {
	    int numRows = src.length;
	    int numCols = src[0].length;
	    float[][] target = new float[numRows][numCols];
	    for (int row = 0; row < numRows; row++) {
		    for (int col = 0; col < numCols; col++) {
		    	target[row][col] = src[row][col];
		    }
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
	public static float[][] copy2DArray(float[][] src) {
	    int length = src.length;
	    float[][] target = new float[length][src[0].length];
	    for (int i = 0; i < length; i++) {
	        System.arraycopy(src[i], 0, target[i], 0, src[i].length);
	    }
	    return target;
	}
	public static void write2DArray(int[][] array, PrintStream ps) {
		int height = array.length;
		for (int r=0; r<height; r++) {
			int width = array[r].length;
			for (int c=0; c<width; c++) {
				ps.print(array[r][c] + " ");
			}
			ps.println();
		}
	}
	public static void write2DArray(byte[][] array, PrintStream ps) {
		int height = array.length;
		int width = array[0].length;
		StringBuilder sb = new StringBuilder();
		Formatter form = new Formatter(sb, Locale.US);
		sb.append("ChannelBytes(w="+width+", h="+height+"):\n");
		for (int row = 0; row < height; row++) {
			for (int col = 0; col < width; col++) {
				byte v = array[row][col];
				form.format("%4d", v).toString();
				if ( !(col == width-1 && row == height-1) ) sb.append(", ");
				if (col == width-1) sb.append('\n');
			}
		}
		form.close();
		ps.print(sb.toString());
	}
	public static void write2DArray(float[][] array, PrintStream ps) {
		int height = array.length;
		for (int r=0; r<height; r++) {
			int width = array[r].length;
			for (int c=0; c<width; c++) {
				ps.print((int)Math.round(array[r][c]) + " ");  // FIXME: toInt
			}
			ps.println();
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

	public static byte[][] roundToByte(float[][] src) {
		int numRows = src.length;
		int numCols = src[0].length;
		byte[][] target = new byte[numRows][numCols];
		for (int row = 0; row < numRows; row++) {
			for (int col= 0; col< numCols; col++) {
				target[row][col] = (byte) Math.round(src[row][col]);
			}
		}
		return target;
	}

}
