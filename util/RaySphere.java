package util;

public class RaySphere extends RayShape {

    double radius;

    public RaySphere() {
        this.radius = 0.1;
    }

    public void setRadius(double r) {
        this.radius = r;
    }

    public double getRadius() {
        return this.radius;
    }
}