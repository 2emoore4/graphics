package week5;

import util.marchten.Moveable;

public class Wood extends Moveable {

    Moveable full, halfOne, halfTwo;
    public boolean cut;

    public Wood() {
        super();
        initShape();
    }

    public void initShape() {
        full = new Moveable();
        full.cube();
        full.scale(2, 0.3, 1);
        add(full);
        cut = false;
    }

    public void cut() {
        remove(full);

        halfOne = new Moveable();
        halfOne.cube();
        halfOne.scale(2, 0.3, 0.5);
        halfOne.translate(0, 0, 0.25);

        halfTwo = new Moveable();
        halfTwo.cube();
        halfTwo.scale(2, 0.3, 0.5);
        halfTwo.translate(0, 0, -0.25);

        add(halfOne);
        add(halfTwo);

        cut = true;
    }

    public void split() {
        halfOne.translate(0, 0, 0.025);
        halfTwo.translate(0, 0, -0.025);
        halfOne.rotateY(-0.001);
        halfTwo.rotateY(0.001);
    }

    public void fall() {
        halfOne.rotateZ(-0.008);
        halfTwo.rotateZ(-0.008);
    }

    public void restart() {
        remove(halfOne);
        remove(halfTwo);
        initShape();
    }
}