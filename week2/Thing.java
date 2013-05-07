package week2;

import java.awt.*;

public class Thing extends Drawable {

    boolean alive;
    Point dest;

    public Thing(int xPos, int yPos, int size) {
        super(xPos, yPos, size);
        this.alive = true;
        this.dest = new Point(xPos, yPos);
    }

    public void moveToDest() {
        double xDiff = dest.getX() - loc.getX();
        double yDiff = dest.getY() - loc.getY();

        move((int) xDiff / 10, (int) yDiff / 10);
    }

    //GETTERS AND SETTERS
    public void setDest(int x, int y) {
        dest.setLocation(x, y);
    }

    public void transDest(int xDiff, int yDiff) {
        dest.translate(xDiff, yDiff);
    }

    public Point getDest() {
        return dest;
    }
}