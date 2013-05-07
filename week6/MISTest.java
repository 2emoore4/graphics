package week6;

import util.ImprovedNoise;
import util.MISApplet;

public class MISTest extends MISApplet {

    double t = 0;
   
    @Override
    public void initFrame(double time) {
        t = 3 * time;
    }

    @Override
    public void setPixel(int x, int y, int rgb[]) {

        int halfX = W / 2;
        int halfY = H / 2;

        int dx = W - x;
        int dy = H - y;

        if (x < halfX) {
            // LEFT SIDE
            if (y < halfY) {
                // TOP LEFT
                rgb[0] = (int)(Math.sin(dx / 5 + t) * 255);
                rgb[1] = (int)(Math.cos(dx / 5 + t) * 255);
                rgb[2] = (int)(Math.tan(dx / 5 + t) * 255);
            } else {
                // BOTTOM LEFT
                rgb[0] = (int)(Math.sin(dy / 5 + t) * 255);
                rgb[1] = (int)(Math.cos(dy / 5 + t) * 255);
                rgb[2] = (int)(Math.tan(dy / 5 + t) * 255);
            }
        }
        else {
            // RIGHT SIDE
            if (y < halfY) {
                // TOP RIGHT
                rgb[0] = (int)(ImprovedNoise.noise(dx, dy, Math.sin(t)) * 255);
                rgb[1] = (int)(ImprovedNoise.noise(dx, dy, Math.cos(t)) * 255);
                rgb[2] = (int)(ImprovedNoise.noise(dx, dy, Math.tan(t)) * 255);
            } else {
                // BOTTOM RIGHT
                rgb[0] = (int)(ImprovedNoise.noise(dx, dy, Math.sin(t)) * 255);
                rgb[1] = (int)(ImprovedNoise.noise(dx, dy, Math.cos(t)) * 255);
                rgb[2] = (int)(ImprovedNoise.noise(dx, dy, Math.tan(t)) * 255);
            }
        }
    }
}