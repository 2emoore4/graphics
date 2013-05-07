package util.marchten;

import java.lang.Math.*;
import java.util.ArrayList;

public class Geometry implements IGeometry {

    int meshM = 15;
    int meshN = 15;
    double[][] vertices;
    int[][] faces;
    Matrix matrix, globMatrix;
    ArrayList<Geometry> children;

    public Geometry() {
        this.matrix = new Matrix();
        matrix.identity();
        this.globMatrix = new Matrix();
        globMatrix.identity();
        children = new ArrayList<Geometry>();
    }

    public void transformByParent(Geometry parent) {
        this.globMatrix.identity();
        parent.getGlobMatrix().copyTo(this.globMatrix);
        this.globMatrix.rightMultiply(this.matrix);
    }

    public Geometry cube() {
        this.vertices = new double[][] {
            { -1,-1,-1}, { 1,-1,-1}, {-1, 1,-1}, {1, 1,-1},
            { -1,-1, 1}, { 1,-1, 1}, {-1, 1, 1}, {1, 1, 1}
        };

        this.faces = new int[][] {
            {0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4},
            {2, 3, 7, 6}, {0, 4, 6, 2}, {1, 5, 7, 3}
        };
        return this;
    }

    private void meshToFaces() {
        for (int m = 0; m < meshM; m++) {
            for (int n = 0; n < meshN; n++) {
                faces[m + (meshM * n)][0] = pointToVertex(m, n);
                faces[m + (meshM * n)][1] = pointToVertex(m + 1, n);
                faces[m + (meshM * n)][2] = pointToVertex(m + 1, n + 1);
                faces[m + (meshM * n)][3] = pointToVertex(m, n + 1);
            }
        }
    }

    public int pointToVertex(int m, int n) {
        return m + ((meshM + 1) * n);
    }

    public Geometry sphere() {
        this.vertices = new double[(meshM + 1) * (meshN + 1)][3];
        this.faces = new int[meshM * meshN][4];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)][0] = Math.cos(2 * Math.PI * u) * Math.cos(Math.PI * (v - 0.5));
                vertices[pointToVertex(m, n)][1] = Math.sin(2 * Math.PI * u) * Math.cos(Math.PI * (v - 0.5));
                vertices[pointToVertex(m, n)][2] = Math.sin(Math.PI * (v - 0.5));
            }
        }

        return this;
    }

    public Geometry torus(double bigR, double littleR) {
        this.vertices = new double[(meshM + 1) * (meshN + 1)][3];
        this.faces = new int[meshM * meshN][4];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)][0] = Math.cos(2 * Math.PI * u) * (bigR + littleR * Math.cos(2 * Math.PI * v));
                vertices[pointToVertex(m, n)][1] = Math.sin(2 * Math.PI * u) * (bigR + littleR * Math.cos(2 * Math.PI * v));
                vertices[pointToVertex(m, n)][2] = littleR * Math.sin(2 * Math.PI * v);
            }
        }

        return this;
    }

    public Geometry cylinder() {
        this.vertices = new double[(meshM + 1) * (meshN + 1)][3];
        this.faces = new int[meshM * meshN][4];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)][0] = Math.cos(2 * Math.PI * u) * r(v);
                vertices[pointToVertex(m, n)][1] = Math.sin(2 * Math.PI * u) * r(v);
                if (v < 0.5) {
                    vertices[pointToVertex(m, n)][2] = -1;
                } else {
                    vertices[pointToVertex(m, n)][2] = 1;
                }
            }
        }

        return this;
    }

    private int r(double v) {
        return v == 0 || v == 1 ? 0 : 1;
    }

    public boolean hasVertex() {
        return this.vertices == null ? false : true;
    }

    public void add(Geometry child) {
        children.add(child);
    }

    public Geometry getChild(int i) {
        return children.get(i);
    }

    public int getNumChildren() {
        return children.size();
    }

    public void remove(Geometry child) {
        children.remove(child);
    }

    public int[] getFace(int i) {
        return faces[i];
    }

    public int[][] getFaces() {
        return this.faces;
    }

    public double[] getVertex(int i) {
        return vertices[i];
    }

    public Matrix getMatrix() {
        return matrix;
    }

    public Matrix getGlobMatrix() {
        return globMatrix;
    }
}