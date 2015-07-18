package util;

import java.io.PrintWriter;

public class Tools {
	private static void showArray(int[][] array, PrintWriter pw) {
		int height = array.length;
		for (int r=0; r<height; r++) {
			int width = array[r].length;
			for (int c=0; c<width; c++) {
				pw.print(array[r][c] + " ");
			}
			pw.println();
		}
	}
}
