package week3;

import java.awt.*;
import java.util.ArrayList;
import util.marchten.BufferedApplet;
import util.marchten.Matrix;

public class Safe extends BufferedApplet {

    Color bgColor = new Color(167, 219, 216, 150);
    Color thingColor = new Color(85, 85, 85, 255);

    int w, h;
    double FL = 10.0;

    double[] point0 = new double[3];
    double[] point1 = new double[3];

    int[] a = new int[2];
    int[] b = new int[2];

    // VERTICES FOR WEIRD SPIKY SHAPE
    double[][] vertices1 = {
        { -1,-1,-1}, { 1,-1,-1}, {-1, 1,-1}, {1, 1,-1},
        { -1,-1, 1}, { 1,-1, 1}, {-1, 1, 1}, {1, 1, 1},
        { 0, 1.5, 0}, {1.5, 0, 0}, {0, -1.5, 0}, {-1.5, 0, 0}
    };

    // VERTICES FOR KINDA INVERTED SHAPE
    double[][] vertices2 = {
        { -1,-1,-1}, { 1,-1,-1}, {-1, 1,-1}, {1, 1,-1},
        { -1,-1, 1}, { 1,-1, 1}, {-1, 1, 1}, {1, 1, 1},
        { 0, 0, 0}
    };

    // EDGES FOR WEIRD SPIKY SHAPE
    int[][] edges1 = {
        {0,1},{2,3},{4,5},{6,7},
        {0,2},{1,3},{4,6},{5,7},
        {0,4},{1,5},{2,6},{3,7},
        {7,8},{6,8},{3,8},{2,8},
        {1,9},{3,9},{5,9},{7,9},
        {0,10},{1,10},{4,10},{5,10},
        {0,11},{2,11},{4,11},{6,11},
    };

    // EDGES FOR KINDA INVERTED SHAPE
    int[][] edges2 = {
        {0,1},{2,3},{4,5},{6,7},
        {0,2},{1,3},{4,6},{5,7},
        {0,4},{1,5},{2,6},{3,7},
        {0,8},{1,8},{2,8},{3,8},{4,8},{6,8},{6,8},{7,8}
    };

    Shape cube1 = new Shape(vertices2, edges2);
    Shape cube2 = new Shape(vertices1, edges1);
    Shape cube3 = new Shape(vertices2, edges2);
    Shape cube4 = new Shape(vertices1, edges1);

    ArrayList<Shape> shapeList;

    Matrix t = new Matrix();

    double startTime = System.currentTimeMillis() / 1000.0;

    public void render(Graphics g) {

        double time = System.currentTimeMillis() / 1000.0 - startTime;

        if (w == 0) {
            shapeList = new ArrayList<Shape>();
            shapeList.add(cube1);
            t.identity();
        }

        w = getWidth();
        h = getHeight();

        g.setColor(bgColor);
        g.fillRect(0, 0, w, h);

        g.setColor(thingColor);

        if (time < 3) {
            // TRANSLATE CUBE1 TO UPPER LEFT CORNER
            t.identity();
            t.translate(-0.03, 0.03, 0.0);
            cube1.getMatrix().rightMultiply(t);
        } else if (time < 6) {
            // TRANSLATE CUBE2 TO UPPER RIGHT CORNER
            if (!shapeList.contains(cube2)) {
                shapeList.add(cube2);
            }
            t.identity();
            t.translate(0.03, 0.03, -0.03);
            cube2.getMatrix().rightMultiply(t);
        } else if (time < 9) {
            // TRANSLATE CUBE3 TO BOTTOM RIGHT CORNER
            if (!shapeList.contains(cube3)) {
                shapeList.add(cube3);
            }
            t.identity();
            t.translate(0.03, -0.03, 0.01);
            cube3.getMatrix().rightMultiply(t);
        } else if (time < 12) {
            // TRANSLATE CUBE4 TO BOTTOM LEFT CORNER
            if (!shapeList.contains(cube4)) {
                shapeList.add(cube4);
            }
            t.identity();
            t.translate(-0.03, -0.03, -0.06);
            cube4.getMatrix().rightMultiply(t);
        } else if (time < 15) {
            // SCALE CUBE 1
            t.identity();
            t.scale(0.99, 0.99, 0.99);
            cube1.getMatrix().rightMultiply(t);
            // SCALE CUBE 3
            t.identity();
            t.scale(1.01, 1.01, 1.01);
            cube3.getMatrix().rightMultiply(t);
            // ROTATE CUBE 2
            t.identity();
            t.rotateZ(0.05);
            cube2.getMatrix().rightMultiply(t);
            // ROTATE CUBE 4
            t.identity();
            t.rotateY(0.05);
            cube4.getMatrix().rightMultiply(t);
        }

        // RENDER SHAPE
        for (Shape s : shapeList) {
            for (int e = 0 ; e < s.getEdges().length ; e++) {
                int i = s.getEdge(e, 0);
                int j = s.getEdge(e, 1);

                s.getMatrix().transform(s.getVertices(i), point0);
                s.getMatrix().transform(s.getVertices(j), point1);

                projectPoint(point0, a);
                projectPoint(point1, b);

                g.drawLine(a[0], a[1], b[0], b[1]);
            }
        }
    }

    public void projectPoint(double[] xyz, int[] pxy) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        pxy[0] = w / 2 + (int)(h * x / (FL - z));
        pxy[1] = h / 2 - (int)(h * y / (FL - z));
    }
}