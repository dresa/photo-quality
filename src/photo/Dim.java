package photo;

public class Dim {
	private int x;
	private int y;
	
	public Dim(int x, int y) {
		this.x = x;
		this.y = y;
	}
	public Dim(Dim d) {
		this.x = d.x;
		this.y = d.y;
	}

	public int getX() { return x; }
	public int getY() { return y; }
	public boolean equals(Dim other) {
		return other.x == this.x && other.y == this.y;
	}
	public String toString() { return "Dim(x=" + x + ", y=" + y + ")"; }
}
