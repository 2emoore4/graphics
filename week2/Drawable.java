package week2;

import java.awt.*;
import java.lang.Math.*;

public class Drawable extends Polygon {

    int w, h, size;
    Point loc;

    public Drawable() {
        super();
    }

    // Constructor for rectangular drawable
    public Drawable(int xPos, int yPos, int size) {
        super();

        w = size;
        h = size;
        this.size = size;

        loc = new Point(xPos, yPos);

        int halfSize = size / 2;
        addPoint(xPos - halfSize, yPos - halfSize);
        addPoint(xPos + halfSize, yPos - halfSize);
        addPoint(xPos + halfSize, yPos + halfSize);
        addPoint(xPos - halfSize, yPos + halfSize);
    }

    public void moveTo(int x, int y) {
        double xDiff = x - loc.getX();
        double yDiff = y - loc.getY();
        translate((int) xDiff, (int) yDiff);
    }

    public void move(int xDiff, int yDiff){
        loc.translate(xDiff, yDiff);
        translate((int) xDiff, (int) yDiff);
    }

    public Point getLoc() {
        return loc;
    }
}