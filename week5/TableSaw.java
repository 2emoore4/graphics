package week5;

import java.awt.*;
import java.util.ArrayList;
import util.marchten.Geometry;
import util.marchten.BufferedApplet;
import util.marchten.Matrix;
import util.marchten.Renderer;
import util.marchten.Moveable;

public class TableSaw extends BufferedApplet {

    Color bgColor = new Color(167, 219, 216, 255);
    Color thingColor = new Color(85, 85, 85, 255);

    int w, h;

    Renderer renderer;

    Moveable tableJoint = new Moveable();
    Moveable tableTop = new Moveable();
    Saw saw = new Saw();
    Wood block = new Wood();

    double startTime = System.currentTimeMillis() / 1000.0;

    public void render(Graphics g) {

        double time = System.currentTimeMillis() / 1000.0 - startTime;

        if (w == 0) {
            renderer = new Renderer(g, getWidth(), getHeight());
            renderer.add(tableJoint);

            tableJoint.add(tableTop.cube());
            tableJoint.add(saw);
            tableJoint.add(block);
            tableJoint.translate(0, -3, -5);
            tableJoint.rotateX(0.2);
            
            tableTop.scale(7, 0.2, 4);

            saw.translate(2, -0.5, -2.5);

            block.translate(-5, 0.15, 0);
        }

        w = getWidth();
        h = getHeight();
        renderer.updateDims(getWidth(), getHeight());

        g.setColor(bgColor);
        g.fillRect(0, 0, w, h);

        saw.rotateZ(0.08);

        if (block.getPosX() < 2) {
            block.translate(0.1, 0, 0);
        } else if (block.getPosX() < 6.75) {
            if (!block.cut) {
                block.cut();
            }

            block.split();
            block.translate(0.1, 0, 0);
        } else if (block.getPosX() < 8) {
            block.translate(0.01, 0, 0);
            block.fall();
        } else if (block.getPosY() > -10) {
            block.translate(0, -0.25, 0);
        } else {
            tableJoint.remove(block);
            block = new Wood();
            tableJoint.add(block);
            block.translate(-5, 0.15, 0);
        }

        g.setColor(thingColor);

        renderer.renderWorld();
    }
}