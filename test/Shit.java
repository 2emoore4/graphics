package test;

import java.awt.*;
import java.awt.event.*;
import java.util.ArrayList;
import util.Geometry;
import util.Moveable;
import util.BufferedApplet;
import util.Matrix;
import util.Face;
import util.ProjectedFace;
import util.Vertex2D;

public class Shit extends BufferedApplet {

    Color bgColor = new Color(167, 219, 216, 150);
    Color thingColor = new Color(85, 85, 85, 255);

    int w, h;
    double FL = 10.0;

    double[] point0 = new double[3];
    double[] point1 = new double[3];

    int[] a = new int[2];
    int[] b = new int[2];

    ZRenderer renderer;

    Moveable bouncer = new Moveable();

    Matrix t = new Matrix();
    Matrix r = new Matrix();

    double startTime = System.currentTimeMillis() / 1000.0;

    public void render(Graphics g) {

        double time = System.currentTimeMillis() / 1000.0 - startTime;

        if (w == 0) {
            renderer = new ZRenderer(getWidth(), getHeight());
            renderer.setGraphics(g);
            bouncer.cube();
            renderer.add(bouncer);
        }

        w = getWidth();
        h = getHeight();

        g.setColor(bgColor);
        g.fillRect(0, 0, w, h);

        g.setColor(thingColor);

        bouncer.rotateY(0.01);

        renderer.renderWorld();
    }
}

class ZRenderer {

    // TEMPORARY VARIABLES
    double[] tmpVertex;
    int[] tmpPoint, tmpColor, tmpBgColor;
    Vertex2D tmpPointD;
    ProjectedFace tmpProjectedFace;
    Vertex2D[][] tmpTriangles;
    Vertex2D[][] tmpTrapezoids;

    int[][] pix;

    Geometry world;

    int w, h;
    double FL = 10.0;

    public ZRenderer(int w, int h) {

        tmpVertex = new double[3];
        tmpPoint = new int[2];
        tmpPointD = new Vertex2D();
        tmpProjectedFace = new ProjectedFace();
        tmpTriangles = new Vertex2D[2][3];
        tmpTrapezoids = new Vertex2D[2][4];

        for (int i = 0; i < tmpTriangles.length; i++) {
            for (int j = 0; j < tmpTriangles[0].length; j++) {
                tmpTriangles[i][j] = new Vertex2D();
            }
        }

        this.pix = new int[w * h][3];

        world = new Geometry();
        this.w = w;
        this.h = h;

        tmpColor = new int[] {0, 0, 0};
        tmpBgColor = new int[] {255, 255, 255};
    }

    public int[] getPix(int x, int y) {
        return pix[x + w * y];
    }

    public void resetPix() {
        for (int i = 0; i < pix.length; i++) {
            pix[i][0] = tmpBgColor[0];
            pix[i][1] = tmpBgColor[1];
            pix[i][2] = tmpBgColor[2];
        }
    }

    public void renderWorld() {
        for (int i = 0; i < world.getNumChildren(); i++) {
            // System.err.println("a child");
            Geometry child = world.getChild(i);
            child.transformByParent(world);
            renderGeometry(child);
        }
    }

    private void renderGeometry(Geometry geo) {

        if (geo instanceof Moveable) {
            Moveable m = (Moveable) geo;
            m.renderPrep();
        }

        if (geo.hasVertex()) {
            processFaces(geo);
        }

        for (int i = 0; i < geo.getNumChildren(); i++) {
            Geometry child = geo.getChild(i);
            child.transformByParent(geo);
            renderGeometry(child);
        }
    }

    private void processFaces(Geometry geo) {

        int testFace = 1;

        for (int f = 0 ; f < geo.getNumFaces(); f++) {
        // for (int f = testFace; f < testFace + 1; f++) {


            // System.err.println("face number " + f);

            for (int v = 0; v < geo.getFace(f).length(); v++) {

                int i = geo.getFace(f).getVertex(v);
                geo.getGlobMatrix().transform(geo.getVertex(i).toArray(), tmpVertex);
                projectPoint(tmpVertex, tmpPoint);

                tmpProjectedFace.setVertex(v, tmpPoint[0], tmpPoint[1]);
            }

            // FACE TO TRIANGLES
            trianglize();
            fuckTriangles();
            buildTrapezoids();
            scanTrapezoids();
        }
    }

    private void trianglize() {

            // tmpProjectedFace.toTriangles(tmpTriangles);

        tmpTriangles[0][0].setX(tmpProjectedFace.getVertex(0).getX());
        tmpTriangles[0][0].setY(tmpProjectedFace.getVertex(0).getY());

        tmpTriangles[0][1].setX(tmpProjectedFace.getVertex(1).getX());
        tmpTriangles[0][1].setY(tmpProjectedFace.getVertex(1).getY());

        tmpTriangles[0][2].setX(tmpProjectedFace.getVertex(2).getX());
        tmpTriangles[0][2].setY(tmpProjectedFace.getVertex(2).getY());


        tmpTriangles[1][0].setX(tmpProjectedFace.getVertex(0).getX());
        tmpTriangles[1][0].setY(tmpProjectedFace.getVertex(0).getY());

        tmpTriangles[1][1].setX(tmpProjectedFace.getVertex(2).getX());
        tmpTriangles[1][1].setY(tmpProjectedFace.getVertex(2).getY());

        tmpTriangles[1][2].setX(tmpProjectedFace.getVertex(3).getX());
        tmpTriangles[1][2].setY(tmpProjectedFace.getVertex(3).getY());


            // System.err.println("triangle 1");
            sortTriangle(tmpTriangles[0]);
            // printTriangle(tmpTriangles[0]);

            // System.err.println("\ntriangle 2");
            sortTriangle(tmpTriangles[1]);
            // printTriangle(tmpTriangles[1]);

            // System.err.println();
    }

    private double area(Vertex2D[] triangle) {
        return areaPiece(0, 1, triangle) + areaPiece(1, 2, triangle) - areaPiece(2, 0, triangle);
    }

    private double areaPiece(int a, int b, Vertex2D[] triangle) {
        int xA = triangle[a].getX();
        int xB = triangle[b].getX();
        int yA = triangle[a].getY();
        int yB = triangle[b].getY();

        return (xA - xB) * (yA + yB) / 2.0;
    }

    private void sortTriangle(Vertex2D[] triangle) {

        if (triangle[0].getY() > triangle[1].getY()) {
            swapVertex(0, 1, triangle);
        }

        if (triangle[1].getY() > triangle[2].getY()) {
            swapVertex(1, 2, triangle);
        }

        if (triangle[0].getY() > triangle[1].getY()) {
            swapVertex(0, 1, triangle);
        }
    }

    private void swapVertex(int i, int j, Vertex2D[] triangle) {
        tmpPoint[0] = triangle[i].getX();
        tmpPoint[1] = triangle[i].getY();
        triangle[i].set(triangle[j].getX(), triangle[j].getY());
        triangle[j].set(tmpPoint[0], tmpPoint[1]);
    }

    private void buildTrapezoids() {

        for (int i = 0; i < tmpTriangles.length; i++) {

            Vertex2D a = tmpTriangles[i][0];
            Vertex2D b = tmpTriangles[i][1];
            Vertex2D c = tmpTriangles[i][2];

            int yA = a.getY();
            int yB = b.getY();
            int yC = c.getY();

            /* I realize that this is terrible, but I have much bigger problems right now,
               so I have to get rid of these divide by zero errors. If I fix everything else,
               I'll come back to this. */
            int t;
            if ((yC - yA) != 0) {
                t = (yB - yA) / (yC - yA);
            } else {
                t = (int) ((yB - yA) / 0.000001);
            }

            int xA = a.getX();
            int xB = b.getX();
            int xC = c.getX();

            int xD = xA + t * (xC - xA);
            int yD = b.getY();
            tmpPointD.setX(xD);
            tmpPointD.setY(yD);

            if (xB < xD) {
                // FIRST SCENARIO
                tmpTrapezoids[0][0] = a;
                tmpTrapezoids[0][1] = a;
                tmpTrapezoids[0][2] = b;
                tmpTrapezoids[0][3] = tmpPointD;

                tmpTrapezoids[1][0] = b;
                tmpTrapezoids[1][1] = tmpPointD;
                tmpTrapezoids[1][2] = c;
                tmpTrapezoids[1][3] = c;
            } else {
                // SECOND SCENARIO
                tmpTrapezoids[0][0] = a;
                tmpTrapezoids[0][1] = a;
                tmpTrapezoids[0][2] = tmpPointD;
                tmpTrapezoids[0][3] = b;

                tmpTrapezoids[1][0] = tmpPointD;
                tmpTrapezoids[1][1] = b;
                tmpTrapezoids[1][2] = c;
                tmpTrapezoids[1][3] = c;
            }
        }
    }

    private void scanTrapezoids() {

        for (int i = 0; i < tmpTrapezoids.length; i++) {

            Vertex2D XLT = tmpTrapezoids[i][0];
            Vertex2D XRT = tmpTrapezoids[i][1];
            Vertex2D XLB = tmpTrapezoids[i][2];
            Vertex2D XRB = tmpTrapezoids[i][3];

            int YT = XLT.getY();
            int YB = XLB.getY();

            // FOR EACH SCANLINE
            for (int y = YT; y < YB; y++) {

                int t = (y - YT) / (YB - YT);
                int XL = (XLT.getX() + t * (XLB.getX() - XLT.getX()));
                int XR = (XRT.getX() + t * (XRB.getX() - XRT.getX()));

                for (int x = XL; x <= XR; x++) {
                    pix[x + w * y][0] = tmpColor[0];
                    pix[x + w * y][1] = tmpColor[1];
                    pix[x + w * y][2] = tmpColor[2];
                }
            }
        }
    }

    private void projectPoint(double[] xyz, int[] pxy) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        pxy[0] = w / 2 + (int)(h * x / (FL - z));
        pxy[1] = h / 2 - (int)(h * y / (FL - z));
    }

    //DEBUGGING SHIT
    public void printTriangle(Vertex2D[] triangle) {
        for (int i = 0; i < triangle.length; i++) {
            System.err.println("v" + i + ": " + triangle[i].getX() + "," + triangle[i].getY());
        }
    }

    // DESPERATE DEBUGGING SHIT
    Graphics g;
    public void setGraphics(Graphics g) {
        this.g = g;
    }

    public void fuckTriangles() {
        for (int i = 0; i < tmpTriangles.length; i++) {
            Vertex2D a = tmpTriangles[i][0];
            Vertex2D b = tmpTriangles[i][1];
            Vertex2D c = tmpTriangles[i][2];

            g.drawLine(a.getX(), a.getY(), b.getX(), b.getY());
            g.drawLine(c.getX(), c.getY(), b.getX(), b.getY());
            g.drawLine(a.getX(), a.getY(), c.getX(), c.getY());
        }
    }

    public void setFocalLength(double l) {
        this.FL = l;
    }

    public void add(Geometry geo) {
        world.add(geo);
    }

    public void remove(Geometry geo) {
        world.remove(geo);
    }
}