package week6;

import java.awt.*;
import java.awt.event.*;
import util.MISApplet;
import util.Moveable;
import util.ZRenderer;

public class MISAttempt extends MISApplet implements KeyListener {

    public void keyPressed(KeyEvent e) {}

    // RESPOND TO KEY PRESSES
    public void keyTyped(KeyEvent e) {
        if (e.getKeyChar() == 'x') {
            // ENABLE/DISABLE X MOVEMENT
            xMove = !xMove;
        } else if (e.getKeyChar() == '2') {
            bouncer.sphere();
        } else if (e.getKeyChar() == '3') {
            bouncer.cylinder();
        } else if (e.getKeyChar() == '4') {
            bouncer.torus(1.0, 0.3);
        } else if (e.getKeyChar() == '1') {
            bouncer.cube();
        }
    }

    public void keyReleased(KeyEvent e) {}

    boolean xMove = true;

    ZRenderer zrenderer;
    Moveable bouncer;
    Moveable ground;

    @Override
    public void initialize() {
        addKeyListener(this);
        requestFocus();

        zrenderer = new ZRenderer(W, H);
        bouncer = new Moveable();
        bouncer.sphere();
        zrenderer.add(bouncer);
        ground = new Moveable();
        ground.cube();
        zrenderer.add(ground);

        if (xMove) {
            bouncer.translate(-8, 3.75, 0);
        } else {
            bouncer.translate(0, 3.75, 0);
        }

        ground.scale(3, 0.2, 3);
        ground.translate(0, -3, -0.5);

        zrenderer.setBGColor(new int[] {167, 219, 216});
        zrenderer.setDrawColor(new int[] {85, 85, 85});
    }
   
    @Override
    public void initFrame(double time) {

        zrenderer.resetPix();
        zrenderer.renderWorld();

        if (time < 1.25) {

            // MOVE THE BOUNCING SHAPE BASED ON THE CURRENT TIME
            if (xMove) {
                bouncer.translate(0.425, bounceVelocity(time), 0);
            } else {
                bouncer.translate(0, bounceVelocity(time), 0);
            }

            // SCALE THE BOUNCING SHAPE BASED ON THE CURRENT TIME
            bouncer.scale(squashXZ(time), squashY(time), squashXZ(time));
        } else if (time > 2) {

            bouncer.reset();

            // RESET THE BOUNCING SHAPE AND RESET THE ANIMATION TIMER
            if (xMove) {
                bouncer.translate(-8, 3.75, 0);
            } else {
                bouncer.translate(0, 3.75, 0);
            }

            restartClock();
        }
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int[] color = zrenderer.getPix(x, y);

        rgb[0] = color[0];
        rgb[1] = color[1];
        rgb[2] = color[2];
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
            return 1.15;
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
}