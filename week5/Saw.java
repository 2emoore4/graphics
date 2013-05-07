package week5;

import java.lang.Math.*;
import util.marchten.Moveable;

public class Saw extends Moveable {

    public Saw() {
        super();
        initShape();
    }

    private void initShape() {
        cylinder();
        scale(1.2, 1.2, 0.2);

        for (int i = 0; i < 15; i++) {
            double angle = (i / 15.0) * 2 * Math.PI;

            Moveable tooth = new Moveable();
            tooth.cube();
            tooth.scale(0.1, 0.1, 0.2);
            tooth.rotateZ(angle);
            tooth.translate(Math.cos(angle) * 1.1, Math.sin(angle) * 1.1, 0);
            
            add(tooth);
        }
    }
}