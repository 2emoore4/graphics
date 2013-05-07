package util.marchten;

/*
Note:

   set(row,col,value) sets the value for a single row and column.
   get(row,col) retrieves the value at a single row and column.

   There are two possible ways to multiply this matrix by another matrix:

      leftMultiply  means:  this <= other times this
      rightMultiply means:  this <= this times other

   transform(src,dst) is how you apply the matrix transformation
   to individual points.
*/

public interface IMatrix
{
   public void identity();
   public void set(int col, int row, double value);
   public double get(int col, int row);
   public void translate(double x, double y, double z);
   public void rotateX(double radians);
   public void rotateY(double radians);
   public void rotateZ(double radians);
   public void scale(double x, double y, double z);
   public void leftMultiply(Matrix other);
   public void rightMultiply(Matrix other);
   public void transform(double[] src, double[] dst);
}

