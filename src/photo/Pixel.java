package photo;

import util.Tools;

public class Pixel {
	public static int toARGB(int a, int r, int g, int b) {
		int m = Tools.MASK_8_BIT;
		return ((a&m) << 24) + ((r&m) << 16)  + ((g&m) << 8) + (b&m);
	}
}
