package util;

import java.awt.*;

public class Renderer {

    Geometry world;
    Graphics g;

    double[] point0;
    double[] point1;

    double FL = 10.0;

    int[] a;
    int[] b;

    int w, h;

    public Renderer(Graphics g, int w, int h) {
        world = new Geometry();
        this.g = g;
        this.w = w;
        this.h = h;
        point0 = new double[3];
        point1 = new double[3];
        a = new int[2];
        b = new int[2];
    }

    public void renderWorld() {
        for (int i = 0; i < world.getNumChildren(); i++) {
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
            for (int f = 0 ; f < geo.getNumFaces(); f++) {
                for (int e = 0; e < geo.getFace(f).length() - 1; e++) {
                    int i = geo.getFace(f).getVertex(e);
                    int j = geo.getFace(f).getVertex(e + 1);

                    geo.getGlobMatrix().transform(geo.getVertex(i).toArray(), point0);
                    geo.getGlobMatrix().transform(geo.getVertex(j).toArray(), point1);

                    projectPoint(point0, a);
                    projectPoint(point1, b);

                    g.drawLine(a[0], a[1], b[0], b[1]);
                }

                int i = geo.getFace(f).getVertex(geo.getFace(f).length() - 1);
                int j = geo.getFace(f).getVertex(0);

                geo.getGlobMatrix().transform(geo.getVertex(i).toArray(), point0);
                geo.getGlobMatrix().transform(geo.getVertex(j).toArray(), point1);

                projectPoint(point0, a);
                projectPoint(point1, b);

                g.drawLine(a[0], a[1], b[0], b[1]);
            }
        }

        for (int i = 0; i < geo.getNumChildren(); i++) {
            Geometry child = geo.getChild(i);
            child.transformByParent(geo);
            renderGeometry(child);
        }
    }

    private void projectPoint(double[] xyz, int[] pxy) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        pxy[0] = w / 2 + (int)(h * x / (FL - z));
        pxy[1] = h / 2 - (int)(h * y / (FL - z));
    }

    public void updateDims(int w, int h) {
        this.w = w;
        this.h = h;
    }

    public void add(Geometry geo) {
        world.add(geo);
    }

    public void remove(Geometry geo) {
        world.remove(geo);
    }
}