package week8;

import java.awt.*;
import java.lang.Math.*;
import util.Geometry;
import util.MISApplet;
import util.Moveable;
import util.ZRenderer;

public class Phong extends MISApplet {

    ZRenderer zrenderer;
    Moveable center, first, second, centerOrbit, firstOrbit, secondOrbit;

    @Override
    public void initialize() {
        zrenderer = new ZRenderer(W, H);
        centerOrbit = new Moveable();
        firstOrbit = new Moveable();
        secondOrbit = new Moveable();

        center = new Moveable();
        center.cylinder();
        centerOrbit.add(center);
        center.scale(1, 1, 1.5);
        center.rotateY(0.2);
        center.getMaterial().setAmbientColor(0.19125,0.0735,0.0225);
        center.getMaterial().setDiffuseColor(0.7038,0.27048,0.0828);
        center.getMaterial().setSpecularColor(0.256777,0.137622,0.086014);
        center.getMaterial().setSpecularPower(12.8);

        first = new Moveable();
        first.sphere();
        firstOrbit.add(first);
        first.translate(3, 0, 0);
        first.getMaterial().setAmbientColor(0.19225,0.19225,0.19225);
        first.getMaterial().setDiffuseColor(0.50754,0.50754,0.50754);
        first.getMaterial().setSpecularColor(0.508273,0.508273,0.508273);
        first.getMaterial().setSpecularPower(76.8);

        second = new Moveable();
        second.cube();
        secondOrbit.add(second);
        second.translate(0, -3, 0);
        second.getMaterial().setAmbientColor(0,0,0);
        second.getMaterial().setDiffuseColor(0.01,0.01,0.01);
        second.getMaterial().setSpecularColor(0.5,0.5,0.5);
        second.getMaterial().setSpecularPower(32.0);

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
        secondOrbit.rotateY(0.05);
        secondOrbit.rotateZ(0.05);
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int[] color = zrenderer.getPix(x, y);

        rgb[0] = color[0];
        rgb[1] = color[1];
        rgb[2] = color[2];
    }
}