package util;

public class Vertex2D {

    int[] coordinates;

    public Vertex2D() {
        this.coordinates = new int[2];
    }

    public void set(int x, int y) {
        setX(x);
        setY(y);
    }

    public int getX() {
        return coordinates[0];
    }

    public void setX(int x) {
        coordinates[0] = x;
    }

    public int getY() {
        return coordinates[1];
    }

    public void setY(int y) {
        coordinates[1] = y;
    }
}