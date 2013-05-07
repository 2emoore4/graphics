package util;

public class RayShape {

    double[] pos;
    Material material;

    public RayShape() {
        this.pos = new double[] {0.0, 0.0, 0.0};
        this.material = new Material();
    }

    public void translate(double x, double y, double z) {
        this.pos[0] += x;
        this.pos[1] += y;
        this.pos[2] += z;
    }

    public void setPosition(double x, double y, double z) {
        this.pos[0] = x;
        this.pos[1] = y;
        this.pos[2] = z;
    }

    public double[] getPosition() {
        return this.pos;
    }

    public Material getMaterial() {
        return this.material;
    }
}