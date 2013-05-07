package util;

import java.lang.Math.*;
import java.util.ArrayList;

public class Geometry implements IGeometry {

    int meshM = 20;
    int meshN = 20;
    Vertex[] vertices;
    Face[] faces;
    Matrix matrix, globMatrix;
    ArrayList<Geometry> children;
    Material material;

    public Geometry() {
        this.matrix = new Matrix();
        matrix.identity();
        this.globMatrix = new Matrix();
        globMatrix.identity();
        children = new ArrayList<Geometry>();
        material = new Material();
    }

    public void transformByParent(Geometry parent) {
        this.globMatrix.identity();
        parent.getGlobMatrix().copyTo(this.globMatrix);
        this.globMatrix.rightMultiply(this.matrix);
    }

    public Geometry cube() {

        this.vertices = new Vertex[] {
            new Vertex(-1, -1, 1, 0, 0, 1), new Vertex(1, -1, 1, 0, 0, 1), new Vertex(1, 1, 1, 0, 0, 1), new Vertex(-1, 1, 1, 0, 0, 1),
            new Vertex(-1, 1, -1, 0, 0, -1), new Vertex(1, 1, -1, 0, 0, -1), new Vertex(1, -1, -1, 0, 0, -1), new Vertex(-1, -1, -1, 0, 0, -1),
            new Vertex(-1, 1, 1, 0, 1, 0), new Vertex(1, 1, 1, 0, 1, 0), new Vertex(1, 1, -1, 0, 1, 0), new Vertex(-1, 1, -1, 0, 1, 0),
            new Vertex(-1, -1, 1, 0, -1, 0), new Vertex(-1, -1, -1, 0, -1, 0), new Vertex(1, -1, -1, 0, -1, 0), new Vertex(1, -1, 1, 0, -1, 0),
            new Vertex(1, -1, 1, 1, 0, 0), new Vertex(1, -1, -1, 1, 0, 0), new Vertex(1, 1, -1, 1, 0, 0), new Vertex(1, 1, 1, 1, 0, 0),
            new Vertex(-1, -1, -1, -1, 0, 0), new Vertex(-1, -1, 1, -1, 0, 0), new Vertex(-1, 1, 1, -1, 0, 0), new Vertex(-1, 1, -1, -1, 0, 0),
        };

        this.faces = new Face[] {
            new Face(0, 1, 2, 3),
            new Face(4, 5, 6, 7),
            new Face(8, 9, 10, 11),
            new Face(12, 13, 14, 15),
            new Face(16, 17, 18, 19),
            new Face(20, 21, 22, 23),
        };

        return this;
    }

    private void meshToFaces() {
        for (int m = 0; m < meshM; m++) {
            for (int n = 0; n < meshN; n++) {
                int currentFace = m + (meshM * n);
                faces[currentFace] = new Face();
                faces[currentFace].set(0, pointToVertex(m, n));
                faces[currentFace].set(1, pointToVertex(m + 1, n));
                faces[currentFace].set(2, pointToVertex(m + 1, n + 1));
                faces[currentFace].set(3, pointToVertex(m, n + 1));
            }
        }
    }

    public int pointToVertex(int m, int n) {
        return m + ((meshM + 1) * n);
    }

    public Geometry sphere() {
        this.vertices = new Vertex[(meshM + 1) * (meshN + 1)];
        this.faces = new Face[meshM * meshN];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)] = new Vertex();
                double x = Math.cos(2 * Math.PI * u) * Math.cos(Math.PI * (v - 0.5));
                double y = Math.sin(2 * Math.PI * u) * Math.cos(Math.PI * (v - 0.5));
                double z = Math.sin(Math.PI * (v - 0.5));
                vertices[pointToVertex(m, n)].set(0, x);
                vertices[pointToVertex(m, n)].set(1, y);
                vertices[pointToVertex(m, n)].set(2, z);
                vertices[pointToVertex(m, n)].set(3, x);
                vertices[pointToVertex(m, n)].set(4, y);
                vertices[pointToVertex(m, n)].set(5, z);
            }
        }

        return this;
    }

    public Geometry torus(double bigR, double littleR) {
        this.vertices = new Vertex[(meshM + 1) * (meshN + 1)];
        this.faces = new Face[meshM * meshN];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)] = new Vertex();

                double x = Math.cos(2 * Math.PI * u) * (bigR + littleR * Math.cos(2 * Math.PI * v));
                double y = Math.sin(2 * Math.PI * u) * (bigR + littleR * Math.cos(2 * Math.PI * v));
                double z = littleR * Math.sin(2 * Math.PI * v);

                vertices[pointToVertex(m, n)].set(0, x);
                vertices[pointToVertex(m, n)].set(1, y);
                vertices[pointToVertex(m, n)].set(2, z);
                vertices[pointToVertex(m, n)].set(3, x);
                vertices[pointToVertex(m, n)].set(4, y);
                vertices[pointToVertex(m, n)].set(5, z);
            }
        }

        return this;
    }

    public Geometry cylinder() {
        this.vertices = new Vertex[(meshM + 1) * (meshN + 1)];
        this.faces = new Face[meshM * meshN];
        meshToFaces();

        for (int m = 0; m < meshM + 1; m++) {
            for (int n = 0; n < meshN + 1; n++) {
                double u = m * 1.0 / meshM;
                double v = n * 1.0 / meshN;
                vertices[pointToVertex(m, n)] = new Vertex();

                double x = Math.cos(2 * Math.PI * u) * r(v);
                double y = Math.sin(2 * Math.PI * u) * r(v);
                double z;

                if (v < 0.5) {
                    z = -1.0;
                } else {
                    z = 1.0;
                }

                vertices[pointToVertex(m, n)].set(0, x);
                vertices[pointToVertex(m, n)].set(1, y);
                vertices[pointToVertex(m, n)].set(2, z);

                if (v == 0.0) {
                    vertices[pointToVertex(m, n)].set(3, 0);
                    vertices[pointToVertex(m, n)].set(4, 0);
                    vertices[pointToVertex(m, n)].set(5, -1);
                } else if (v == 1.0) {
                    vertices[pointToVertex(m, n)].set(3, 0);
                    vertices[pointToVertex(m, n)].set(4, 0);
                    vertices[pointToVertex(m, n)].set(5, 1);
                } else {
                    vertices[pointToVertex(m, n)].set(3, Math.cos(2 * Math.PI * u));
                    vertices[pointToVertex(m, n)].set(4, Math.sin(2 * Math.PI * u));
                    vertices[pointToVertex(m, n)].set(5, 0);
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

    public Face getFace(int i) {
        return faces[i];
    }

    public int getNumFaces() {
        return this.faces.length;
    }

    public Vertex getVertex(int i) {
        return vertices[i];
    }

    public Matrix getMatrix() {
        return matrix;
    }

    public Matrix getGlobMatrix() {
        return globMatrix;
    }

    public Material getMaterial() {
        return material;
    }
}