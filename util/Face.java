package util;

public class Face {

    int[] vertices;

    public Face() {
        this.vertices = new int[4];
    }

    public Face(int zero, int one, int two, int three) {
        this.vertices = new int[] {zero, one, two, three};
    }

    public void set(int pos, int vertex) {
        vertices[pos] = vertex;
    }

    public int getVertex(int pos) {
        return vertices[pos];
    }

    public int length() {
        return vertices.length;
    }
}