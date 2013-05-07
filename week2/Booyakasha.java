package week2;

import java.awt.*;
import java.util.ArrayList;
import java.lang.Math.*;
import util.marchten.BufferedApplet;

public class Booyakasha extends BufferedApplet {

    int w = 0, h = 0;
    Color bgColor = new Color(167, 219, 216, 150);
    Color thingColor = new Color(85, 85, 85, 255);
    Color bombColor = new Color(105, 210, 231, 255);

    double startTime;
    Point mouse = new Point();
    ArrayList<Drawable> drawableList;

    Thing thingy = new Thing(100, 100, 50);

    public boolean mouseDown(Event e, int x, int y) {
        drawableList.add(new Projectile(x, y, 60));
        return true;
    }

    public boolean mouseMove(Event e, int x, int y) {
        mouse.setLocation(x, y);
        return true;
    }

    public boolean keyDown(Event e, int key) {
        switch (key) {
            case 'c':
                drawableList.add(new Thing((int) mouse.getX(), (int) mouse.getY(), 50));
                break;
            case 'r':
                drawableList.clear();
                break;
        }

        return true;
    }

    public void render(Graphics g) {

        if (w == 0) {

            // INIT STUFF
            startTime = System.currentTimeMillis() / 1000.0;
            drawableList = new ArrayList<Drawable>();

            drawableList.add(thingy);
        }

        w = getWidth();
        h = getHeight();

        // FILL BACKGROUND
        g.setColor(bgColor);
        g.fillRect(0, 0, w, h);

        // TIME COUNTER
        // g.setColor(bombColor);
        // double time = System.currentTimeMillis() / 1000.0 - startTime;
        // g.drawString("" + time, w - 45, h - 2);

        // FIND NEW POSITIONS
        moveDrawables();

        // REDRAW EVERYTHING
        repaintDrawables(g);
    }

    private void repaintDrawables(Graphics g) {
        for (Drawable d : drawableList) {

            if (d instanceof Projectile) {
                // DO PROJECTILE STUFF
                Projectile p = (Projectile) d;
                g.setColor(new Color(224, 228, 204, p.getVisibility()));
                int xSize = (int) (p.w * p.scale);
                int ySize = (int) (p.h * p.scale);
                g.fillOval((int) p.getLoc().getX() - xSize / 2, (int) p.getLoc().getY() - ySize / 2,
                    xSize, ySize);
            } else if (d instanceof Thing) {
                g.setColor(thingColor);
                g.drawPolygon(d);
            }
            
        }
    }

    private void moveDrawables() {
        int bombIndex = -1;
        for (Drawable d : drawableList) {
            if (d instanceof Thing) {
                // LOOK FOR ENEMIES, SET NEW DESTINATION                 
                Thing t = (Thing) d;
                avoidEnemies(t);
                t.moveToDest();
            } else if (d instanceof Projectile) {
                Projectile p = (Projectile) d;
                if (p.scaling()) {
                    p.enlarge();
                } else {
                    // EXPLODE AT THE LOCATION OF THE BOMB
                    explode(p.getLoc());
                    bombIndex = drawableList.indexOf(d);
                }
            }
        }

        if (bombIndex != -1) {
            drawableList.remove(bombIndex);
        }
    }

    private void avoidEnemies(Thing t) {

        // AVOID WALLS
        avoidWalls(t);

        // AVOID CURSOR
        avoidCursor(t);

        // AVOID PROJECTILES
        avoidBombs(t);
    }

    private void avoidWalls(Thing t) {

        Point pos = t.getLoc();
        int x = (int) pos.getX();
        int y = (int) pos.getY();

        // IF LEFT
        if ((x - (t.w / 2)) < 10) {
            t.transDest(20, 0);
        }

        // IF RIGHT
        if ((x + (t.w / 2)) > w - 10) {
            t.transDest(-20, 0);
        }

        // IF TOP
        if ((y - (t.h / 2)) < 10) {
            t.transDest(0, 20);
        }

        // IF BOTTOM
        if ((y + (t.h / 2)) > h - 10) {
            t.transDest(0, -20);
        }
    }

    private void avoidCursor(Thing t) {
        repel(t, mouse, 0.75);
    }

    private void avoidBombs(Thing t) {
        for (Drawable d : drawableList) {
            if (d instanceof Projectile) {
                repel(t, d.getLoc(), 1.25);
            }
        }
    }

    private void repel(Thing t, Point p, double severity) {
        double dist = distanceBetween(p, t.getLoc());
        int force = fleeForce(dist, severity);
        double xTrans = force * xBetweenNorm(p, t.getLoc());
        double yTrans = force * yBetweenNorm(p, t.getLoc());
        t.transDest((int) xTrans, (int) yTrans);
    }

    private void explode(Point p) {
        // EXPLODE AT POINT P

        for (Drawable d : drawableList) {
            if (d instanceof Thing) {
                Thing t = (Thing) d;
                // FIND DISTANCE FROM BOMB, THEN TRANSLATE THE DESTINATION BASED ON THAT DISTANCE
                repel(t, p, 100.0);
            }
        }
    }

    private int fleeForce(double dist, double range) {
        // CALCULATE EXPLOSION MULTIPLIER BASED ON DISTANCE FROM BOMB

        int force;

        if (range > 2.0) {
            // IT'S A BOMB
            if (dist < 100.0) {
                force =  0;
            } else if (dist < 200.0) {
                force =  100;
            } else if (dist < 250.0) {
                force =  75;
            } else if (dist < 300.0) {
                force =  50;
            } else if (dist < 400.0) {
                force =  25;
            } else {
                force = 0;
            }
        } else {
            if (dist < 50.0 * range) {
                force =  5;
            } else if (dist < 65.0 * range) {
                force =  4;
            } else if (dist < 85.0 * range) {
                force =  3;
            } else if (dist < 115.0 * range) {
                force =  2;
            } else {
                force =  0;
            }
        }

        return force;
    }

    private double distanceBetween(Point a, Point b) {
        double xDiff = b.getX() - a.getX();
        double yDiff = b.getY() - a.getY();

        return Math.sqrt((xDiff * xDiff) + (yDiff * yDiff));
    }

    private int xBetweenNorm(Point a, Point b) {
        if (b.getX() - a.getX() > 0) {
            return 1;
        } else {
            return -1;
        }
    }

    private int yBetweenNorm(Point a, Point b) {
        if (b.getY() - a.getY() > 0) {
            return 1;
        } else {
            return -1;
        }
    }
}