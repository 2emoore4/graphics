package util;

public class Vertex {

    double[] coordinates;

    public Vertex() {
        this.coordinates = new double[6];
    }

    public Vertex(double x, double y, double z) {
        this.coordinates = new double[6];
        set(0, x);
        set(1, y);
        set(2, z);
    }

    public Vertex(double x, double y, double z, double nx, double ny, double nz) {
        this.coordinates = new double[6];
        set(0, x);
        set(1, y);
        set(2, z);
        set(3, nx);
        set(4, ny);
        set(5, nz);
    }

    public void set(int pos, double value) {
        this.coordinates[pos] = value;
    }

    public double[] toArray() {
        return this.coordinates;
    }
}