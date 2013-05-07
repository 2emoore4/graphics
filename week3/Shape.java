package week3;

import util.marchten.Matrix;

public class Shape {

    private double[][] vertices;
    private int[][] edges;
    private Matrix matrix;

    public Shape(double[][] vertices, int[][] edges) {
        this.vertices = vertices;
        this.edges = edges;
        this.matrix = new Matrix();
        matrix.identity();
    }

    // THESE GETTERS ARE USED FOR RENDERING
    public int getEdge(int i, int j) {
        return edges[i][j];
    }

    public int[][] getEdges() {
        return this.edges;
    }

    public double[] getVertices(int i) {
        return vertices[i];
    }

    // THIS GETTER IS USED FOR TRANSFORMATIONS
    public Matrix getMatrix() {
        return matrix;
    }
}