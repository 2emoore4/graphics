package week7;

import java.awt.*;
import java.lang.Math.*;
import util.Geometry;
import util.MISApplet;
import util.Moveable;
import util.ZRenderer;

public class ZBuf extends MISApplet {

    ZRenderer zrenderer;
    Moveable center, first, second, centerOrbit, firstOrbit, secondOrbit;

    @Override
    public void initialize() {
        zrenderer = new ZRenderer(W, H);
        zrenderer.showNormals(true);
        centerOrbit = new Moveable();
        firstOrbit = new Moveable();
        secondOrbit = new Moveable();

        center = new Moveable();
        center.cylinder();
        centerOrbit.add(center);

        first = new Moveable();
        first.cube();
        firstOrbit.add(first);
        first.scale(1, 1, 1);
        first.translate(3, 0, 0);

        second = new Moveable();
        second.sphere();
        secondOrbit.add(second);
        second.scale(1, 1, 1);
        second.translate(0, -3, 0);

        zrenderer.add(centerOrbit);
        zrenderer.add(firstOrbit);
        zrenderer.add(secondOrbit);

        zrenderer.setBGColor(new int[] {167, 219, 216});

        zrenderer.resetPix();
        zrenderer.renderWorld();
    }
   
    @Override
    public void initFrame(double time) {
        zrenderer.resetPix();
        zrenderer.renderWorld();

        centerOrbit.rotateX(0.1);
        firstOrbit.rotateY(0.1);
        secondOrbit.rotateY(0.15);
        secondOrbit.rotateZ(0.15);
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int[] color = zrenderer.getPix(x, y);

        rgb[0] = color[0];
        rgb[1] = color[1];
        rgb[2] = color[2];
    }
}