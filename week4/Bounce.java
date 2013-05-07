package week4;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import util.marchten.Geometry;
import util.marchten.BufferedApplet;
import util.marchten.Matrix;

public class Bounce extends BufferedApplet implements KeyListener {

    public void keyPressed(KeyEvent e) {}

    // RESPOND TO KEY PRESSES
    public void keyTyped(KeyEvent e) {
        if (e.getKeyChar() == 'x') {
            // ENABLE/DISABLE X MOVEMENT
            xMove = !xMove;
        } else if (e.getKeyChar() == '1') {
            bouncer.sphere();
        } else if (e.getKeyChar() == '2') {
            bouncer.cylinder();
        } else if (e.getKeyChar() == '3') {
            bouncer.torus(1.0, 0.3);
        } else if (e.getKeyChar() == '4') {
            bouncer.cube();
        }
    }

   public void keyReleased(KeyEvent e) {}

    Color bgColor = new Color(167, 219, 216, 150);
    Color thingColor = new Color(85, 85, 85, 255);

    int w, h;
    double FL = 10.0;

    double[] point0 = new double[3];
    double[] point1 = new double[3];

    int[] a = new int[2];
    int[] b = new int[2];

    boolean xMove = true;

    Geometry bouncer = new Geometry();
    Geometry ground = new Geometry();

    ArrayList<Geometry> shapeList;

    Matrix t = new Matrix();
    Matrix r = new Matrix();

    double startTime = System.currentTimeMillis() / 1000.0;

    public void render(Graphics g) {

        double time = System.currentTimeMillis() / 1000.0 - startTime;

        if (w == 0) {
            addKeyListener(this);
            requestFocus();

            shapeList = new ArrayList<Geometry>();
            bouncer.sphere();
            shapeList.add(bouncer);
            ground.cube();
            shapeList.add(ground);

            // INITIALIZE BOUNCING SHAPE
            t.identity();
            if (xMove) {
                t.translate(-8, 3.75, 0);
            } else {
                t.translate(0, 3.75, 0);
            }
            bouncer.getMatrix().rightMultiply(t);

            //INITIALIZE GROUND BOX SHAPE
            t.identity();
            t.scale(3, 0.2, 3);
            t.translate(0, -3, -0.5);
            ground.getMatrix().rightMultiply(t);
        }

        w = getWidth();
        h = getHeight();

        g.setColor(bgColor);
        g.fillRect(0, 0, w, h);

        g.setColor(thingColor);

        if (time < 1.25) {

            // MOVE THE BOUNCING SHAPE BASED ON THE CURRENT TIME
            t.identity();

            if (xMove) {
                t.translate(0.425, bounceVelocity(time), 0);
            } else {
                t.translate(0, bounceVelocity(time), 0);
            }

            // SCALE THE BOUNCING SHAPE BASED ON THE CURRENT TIME
            t.scale(squashXZ(time), squashY(time), squashXZ(time));

            // bouncer.getMatrix().rightMultiply(t);

            // r.identity();
            // r.rotateZ(0.01);
            bouncer.getMatrix().rightMultiply(t);
            // bouncer.getMatrix().rightMultiply(r);
        } else if (time > 2) {

            // RESET THE BOUNCING SHAPE AND RESET THE ANIMATION TIMER
            bouncer.getMatrix().identity();
            t.identity();

            if (xMove) {
                t.translate(-8, 3.75, 0);
            } else {
                t.translate(0, 3.75, 0);
            }

            bouncer.getMatrix().rightMultiply(t);

            startTime = System.currentTimeMillis() / 1000.0;
        }

        for (Geometry s : shapeList) {
            for (int f = 0 ; f < s.getFaces().length; f++) {
                for (int e = 0; e < s.getFace(f).length - 1; e++) {
                    int i = s.getFace(f)[e];
                    int j = s.getFace(f)[e + 1];

                    s.getMatrix().transform(s.getVertex(i), point0);
                    s.getMatrix().transform(s.getVertex(j), point1);

                    projectPoint(point0, a);
                    projectPoint(point1, b);

                    g.drawLine(a[0], a[1], b[0], b[1]);
                }

                int i = s.getFace(f)[s.getFace(f).length - 1];
                int j = s.getFace(f)[0];

                s.getMatrix().transform(s.getVertex(i), point0);
                s.getMatrix().transform(s.getVertex(j), point1);

                projectPoint(point0, a);
                projectPoint(point1, b);

                g.drawLine(a[0], a[1], b[0], b[1]);
            }
        }
    }

    private double bounceVelocity(double time) {
        // CALCULATE THE Y VELOCITY OF THE BOUNCING SHAPE BASED ON TIME
        if (time < 0.625) {
            return -1.05 * time;
        } else {
            return 1.087 - (0.87 * time);
        }
    }

    private double squashY(double time) {
        // CALCULATE THE Y DIMENSION SCALE FACTOR BASED ON TIME
        if (time > 0.5 && time < 0.625) {
            return 0.85;
        } else if (time > 0.625 && time < 0.75) {
            return 1.18;
        }
        return 1.0;
    }

    private double squashXZ(double time) {
        // CALCULATE THE X AND Z DIMENSION SCALE FACTORS BASED ON TIME
        if (time > 0.5 && time < 0.625) {
            return 1.05;
        } else if (time > 0.625 && time < 0.75) {
            return 0.95;
        }
        return 1.0;
    }

    private void projectPoint(double[] xyz, int[] pxy) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        pxy[0] = w / 2 + (int)(h * x / (FL - z));
        pxy[1] = h / 2 - (int)(h * y / (FL - z));
    }

}