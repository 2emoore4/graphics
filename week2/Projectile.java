package week2;

import java.awt.*;
import java.lang.Math.*;

public class Projectile extends Drawable {

    int visibility;
    double scale;

    public Projectile(int xPos, int yPos, int size) {
        super(xPos, yPos, size);
        scale = 0.01;
        visibility = 0;
    }

    public void enlarge() {
        scale += 0.02;
        visibility += 5;
    }

    public boolean scaling() {
        return (scale < 1.0);
    }

    public double getScale() {
        return scale;
    }

    public int getVisibility() {
        return visibility;
    }
}