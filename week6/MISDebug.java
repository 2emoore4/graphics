package week6;

import java.awt.*;
import util.MISApplet;
import util.Moveable;
import util.ZRenderer;

public class MISDebug extends MISApplet {

    ZRenderer zrenderer;
    Moveable shape;

    int[] bgColor = new int[] {255, 255, 255};

    @Override
    public void initialize() {

        zrenderer = new ZRenderer(W, H);
        shape = new Moveable();
        shape.torus(1.0, 0.2);
        zrenderer.add(shape);

        zrenderer.setBGColor(new int[] {167, 219, 216});
        zrenderer.setDrawColor(new int[] {85, 85, 85});

        zrenderer.resetPix();
        zrenderer.renderWorld();
    }
   
    @Override
    public void initFrame(double time) {

        zrenderer.resetPix();
        zrenderer.renderWorld();

        shape.rotateY(0.05);
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int[] color = zrenderer.getPix(x, y);

        rgb[0] = color[0];
        rgb[1] = color[1];
        rgb[2] = color[2];
    }
}