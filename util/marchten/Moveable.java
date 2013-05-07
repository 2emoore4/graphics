package util.marchten;

public class Moveable extends Geometry {

    double[] rotate, translate, scale;
    Matrix t;

    public Moveable() {
        super();
        this.rotate = new double[3];
        this.translate = new double[3];
        this.scale = new double[3];
        this.t = new Matrix();

        reset();
    }

    public void renderPrep() {

        matrix.identity();

        if (translate[0] != 0 || translate[1] != 0 || translate[2] != 0) {
            t.identity();
            t.translate(translate[0], translate[1], translate[2]);
            matrix.rightMultiply(t);
        }

        if (rotate[0] != 0) {
            t.identity();
            t.rotateX(rotate[0]);
            matrix.rightMultiply(t);
        }

        if (rotate[1] != 0) {
            t.identity();
            t.rotateY(rotate[1]);
            matrix.rightMultiply(t);
        }

        if (rotate[2] != 0) {
            t.identity();
            t.rotateZ(rotate[2]);
            matrix.rightMultiply(t);
        }

        if (scale[0] != 1 || scale[1] != 1 || scale[2] != 1) {
            t.identity();
            t.scale(scale[0], scale[1], scale[2]);
            matrix.rightMultiply(t);
        }
    }

    public void reset() {
        for (int i = 0; i < rotate.length; i++) {
            rotate[i] = 0;
        }

        for (int i = 0; i < translate.length; i++) {
            translate[i] = 0;
        }

        for (int i = 0; i < scale.length; i++) {
            scale[i] = 1;
        }
    }

    public void setPosition(double x, double y, double z) {
        translate[0] = x;
        translate[1] = y;
        translate[2] = z;
    }

    public void translate(double x, double y, double z) {
        translate[0] += x;
        translate[1] += y;
        translate[2] += z;
    }

    public void rotateX(double radians) {
        rotate[0] += radians;
    }

    public void rotateY(double radians) {
        rotate[1] += radians;

        // matrix.identity();
        // t.identity();
        // t.rotateY(radians);
        // matrix.rightMultiply(t);
    }

    public void rotateZ(double radians) {
        rotate[2] += radians;
    }

    public void scale(double x, double y, double z) {
        scale[0] = scale[0] * x;
        scale[1] = scale[1] * y;
        scale[2] = scale[2] * z;
    }

    public double getPosX() {
        return translate[0];
    }

    public double getPosY() {
        return translate[1];
    }

    public double getPosZ() {
        return translate[2];
    }
}