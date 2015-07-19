package photo;

public class Pixel {
	public static int toARGB(int a, int r, int g, int b) {
		return (a << 24) + (r << 16)  + (g << 8) + b;
	}
}
