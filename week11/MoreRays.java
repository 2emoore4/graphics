package week11;

import java.awt.*;
import java.awt.event.*;
import util.MISApplet;
import util.RaySphere;
import util.RayTracer;

public class MoreRays extends MISApplet implements KeyListener {

    public void keyPressed(KeyEvent e) {}

    // RESPOND TO KEY PRESSES
    public void keyTyped(KeyEvent e) {
        if (e.getKeyChar() == 'a') {
            rayTracer.toggleAA();
        } else if (e.getKeyChar() == 's') {
            rayTracer.toggleSS();
        } else if (e.getKeyChar() == 't') {
            rayTracer.toggleTextures();
        } else if (e.getKeyChar() == 'y') {
            rayTracer.setTexture(0);
        } else if (e.getKeyChar() == 'u') {
            rayTracer.setTexture(1);
        } else if (e.getKeyChar() == 'i') {
            rayTracer.setTexture(2);
        } else if (e.getKeyChar() == '2') {
            rayTracer.setSSRate(2);
        } else if (e.getKeyChar() == '4') {
            rayTracer.setSSRate(4);
        } else if (e.getKeyChar() == '6') {
            rayTracer.setSSRate(6);
        } else if (e.getKeyChar() == '8') {
            rayTracer.setSSRate(8);
        } else if (e.getKeyChar() == 'r') {
            rayTracer.resetStuff();
        } else if (e.getKeyChar() == 'l') {
            lights = !lights;;
        } else if (e.getKeyChar() == 'h') {
            rayTracer.toggleShadows();
        }
    }

    public void keyReleased(KeyEvent e) {}

    RayTracer rayTracer;
    RaySphere sphere, sphere2, sphere3;
    boolean lights = false;

    @Override
    public void initialize() {

        addKeyListener(this);
        requestFocus();

        rayTracer = new RayTracer(W, H);
        rayTracer.setBGColor(new int[] {152, 205, 202});

        // RED PLASTIC
        sphere = new RaySphere();
        sphere.setRadius(0.075);
        sphere.getMaterial().setAmbientColor(0.0, 0.0, 0.0);
        sphere.getMaterial().setDiffuseColor(0.5,0.0,0.0);
        sphere.getMaterial().setSpecularColor(0.7,0.6,0.6);
        sphere.getMaterial().setSpecularPower(32.0);
        rayTracer.add(sphere);

        // CYAN PLASTIC
        sphere2 = new RaySphere();
        sphere2.setRadius(0.075);
        sphere2.getMaterial().setAmbientColor(0.0, 0.1, 0.06);
        sphere2.getMaterial().setDiffuseColor(0.0,0.50980392,0.50980392);
        sphere2.getMaterial().setSpecularColor(0.50196078,0.50196078,0.50196078);
        sphere2.getMaterial().setSpecularPower(32.0);
        rayTracer.add(sphere2);

        // COPPER
        sphere3 = new RaySphere();
        sphere3.getMaterial().setAmbientColor(0.19125, 0.0735, 0.0225);
        sphere3.getMaterial().setDiffuseColor(0.7038, 0.27048, 0.0828);
        sphere3.getMaterial().setSpecularColor(0.256777, 0.137622, 0.086014);
        sphere3.getMaterial().setSpecularPower(12.8);
        rayTracer.add(sphere3);

        rayTracer.resetPix();
        rayTracer.renderWorld();
    }

    @Override
    public void initFrame(double time) {

        if (lights) {
            rayTracer.getLights()[2][0][0] = Math.sin(time * 0.5);
            rayTracer.getLights()[2][0][1] = Math.cos(time * 0.5);

            rayTracer.getLights()[1][1][0] += Math.sin(time * 2) * 0.25;
            rayTracer.getLights()[1][1][2] += Math.cos(time * 2) * 0.25;
        }

        sphere.setPosition(0.0,
                           -0.2 * Math.cos(time),
                           -2.2 * Math.sin(time));
        sphere2.setPosition(0.2 * Math.cos(time),
                            0.0,
                            2.2 * Math.sin(time));
        rayTracer.resetPix();
        rayTracer.renderWorld();
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int[] color = rayTracer.getPix(x, y);

        rgb[0] = color[0];
        rgb[1] = color[1];
        rgb[2] = color[2];
    }
}