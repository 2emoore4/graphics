package util;

public interface IGeometry {
    public void add(Geometry child);
    public Geometry getChild(int i);
    public Matrix getMatrix();
    public int getNumChildren();
    public void remove(Geometry child);
}