package util;

import java.lang.Math.*;

public class ZRenderer {

    // TEMPORARY VARIABLES
    double[] tmp2DVertex, tmp3DVertex, tmpDoubleVertex, tmpVec, tmpCol;
    int[] tmpColor, tmpBgColor;
    double[][] tmpProjectedFace;
    double[][][] tmpTriangles, tmpTrapezoids;
    boolean[] positive;
    boolean showNormals;

    // TEMPORARY PIXEL ARRAY AND ZBUF ARRAY
    int[][] pix;
    double[][] zbuf;

    // ARRAY OF LIGHTS
    double[][][] lights = {
        { { 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0} },
        { {-1.0,-1.0,-1.0}, {1.0, 1.0, 1.0} },
    };
    double[] eyeDir = {0.0, 0.0, 1.0};
    Material tmpMaterial;

    Geometry world;

    int w, h;
    double FL = 10.0;

    public ZRenderer(int w, int h) {

        // INITIALIZE TEMP VARIABLES
        tmp3DVertex = new double[6];
        tmpDoubleVertex = new double[9];
        tmpVec = new double[3];
        tmpCol = new double[3];
        tmp2DVertex = new double[9];
        tmpProjectedFace = new double[4][9];
        tmpTriangles = new double[2][3][9];
        tmpTrapezoids = new double[4][4][9];

        positive = new boolean[2];

        showNormals = false;

        this.pix = new int[w * h][3];
        this.zbuf = new double[w][h];

        world = new Geometry();
        this.w = w;
        this.h = h;
    }

    // USED BY MISAPPLET TO GET THE RGB FOR A SPECIFIC PIXEL
    public int[] getPix(int x, int y) {
        return pix[x + w * y];
    }

    // RESETS THE PIXEL ARRAY TO THE BACKGROUND COLOR
    public void resetPix() {
        for (int i = 0; i < pix.length; i++) {
            pix[i][0] = tmpBgColor[0];
            pix[i][1] = tmpBgColor[1];
            pix[i][2] = tmpBgColor[2];
        }

        for (int i = 0; i < zbuf.length; i++) {
            for (int j = 0; j < zbuf[0].length; j++) {
                zbuf[i][j] = -FL;
            }
        }
    }

    // RENDERS EACH CHILD OF THE 'WORLD' GEOMETRY
    public void renderWorld() {
        for (int i = 0; i < world.getNumChildren(); i++) {
            Geometry child = world.getChild(i);
            child.transformByParent(world);
            renderGeometry(child);
        }
    }

    // RENDERS THE FACES OF A GEOMETRY, THEN RENDERS ITS CHILDREN
    private void renderGeometry(Geometry geo) {

        tmpMaterial = geo.getMaterial();

        if (geo instanceof Moveable) {
            Moveable m = (Moveable) geo;
            m.renderPrep();
        }

        if (geo.hasVertex()) {
            processFaces(geo);
        }

        for (int i = 0; i < geo.getNumChildren(); i++) {
            Geometry child = geo.getChild(i);
            child.transformByParent(geo);
            renderGeometry(child);
        }
    }

    // PROJECTS FACES, CONVERTS TO TRIANGLES, CONVERTS TO TRAPEZOIDS, SCANS TRAPEZOIDS
    private void processFaces(Geometry geo) {

        for (int f = 0 ; f < geo.getNumFaces(); f++) {

            // PROJECT FACES
            for (int v = 0; v < geo.getFace(f).length(); v++) {

                int i = geo.getFace(f).getVertex(v);
                geo.getGlobMatrix().transform(geo.getVertex(i).toArray(), tmp3DVertex);
                projectPoint(tmp3DVertex, tmp2DVertex);
                lightPoint(tmp2DVertex, tmpMaterial);

                copy2DVertex(tmp2DVertex, tmpProjectedFace[v]);
            }

            trianglize();
            buildTrapezoids();
            scanTrapezoids();

        }
    }

    // CONVERTS PROJECTED FACES TO TRIANGLES
    private void trianglize() {

        double[] vertexZero = tmpProjectedFace[0];
        double[] vertexOne = tmpProjectedFace[1];
        double[] vertexTwo = tmpProjectedFace[2];
        double[] vertexThree = tmpProjectedFace[3];

        double[][] triangleZero = tmpTriangles[0];
        double[][] triangleOne = tmpTriangles[1];

        copy2DVertex(vertexZero, triangleZero[0]);
        copy2DVertex(vertexOne, triangleZero[1]);
        copy2DVertex(vertexTwo, triangleZero[2]);

        copy2DVertex(vertexZero, triangleOne[0]);
        copy2DVertex(vertexTwo, triangleOne[1]);
        copy2DVertex(vertexThree, triangleOne[2]);

        if (area(triangleZero) < 0) {
            positive[0] = true;
        } else {
            positive[0] = false;
        }
        if (area(triangleOne) < 0) {
            positive[1] = true;
        } else {
            positive[1] = false;
        }

        sortTriangle(triangleZero);
        sortTriangle(triangleOne);
    }

    // FINDS THE AREA OF A TRIANGLE
    private double area(double[][] triangle) {

        double xA = triangle[0][0];
        double xB = triangle[1][0];
        double yA = triangle[0][1];
        double yB = triangle[1][1];

        double a1 = (xA - xB) * (yA + yB) / 2.0;

        xA = triangle[1][0];
        xB = triangle[2][0];
        yA = triangle[1][1];
        yB = triangle[2][1];

        double a2 = (xA - xB) * (yA + yB) / 2.0;

        xA = triangle[2][0];
        xB = triangle[0][0];
        yA = triangle[2][1];
        yB = triangle[0][1];

        double a3 = (xA - xB) * (yA + yB) / 2.0;

        return a1 + a2 + a3;
    }

    // SORTS THE VERTICES OF A TRIANGLE FROM LOW TO HIGH IN THE Y DIRECTION
    private void sortTriangle(double[][] triangle) {

        if (triangle[0][1] > triangle[1][1]) {
            swapVertex(0, 1, triangle);
        }

        if (triangle[1][1] > triangle[2][1]) {
            swapVertex(1, 2, triangle);
        }

        if (triangle[0][1] > triangle[1][1]) {
            swapVertex(0, 1, triangle);
        }
    }

    // HELPER SWAP METHOD
    private void swapVertex(int i, int j, double[][] triangle) {
        copy2DVertex(triangle[i], tmpDoubleVertex);
        copy2DVertex(triangle[j], triangle[i]);
        copy2DVertex(tmpDoubleVertex, triangle[j]);
    }

    // CONVERTS TRIANGLES TO TRAPEZOIDS
    private void buildTrapezoids() {

        for (int i = 0; i < tmpTriangles.length; i++) {

            // ONLY MAKE TRAPEZOIDS IF TRIANGLE HAS POSITIVE AREA
            if (positive[i]) {
                double[] a = tmpTriangles[i][0];
                double[] b = tmpTriangles[i][1];
                double[] c = tmpTriangles[i][2];

                double yA = a[1];
                double yB = b[1];
                double yC = c[1];

                double t = (b[1] - a[1]) / (c[1] - a[1]);

                double xD = interpolate(t, a[0], c[0]);
                double yD = b[1];
                double zD = interpolate(t, a[2], c[2]);
                double nxD = interpolate(t, a[3], c[3]);
                double nyD = interpolate(t, a[4], c[4]);
                double nzD = interpolate(t, a[5], c[5]);
                double rD = interpolate(t, a[6], c[6]);
                double gD = interpolate(t, a[7], c[7]);
                double bD = interpolate(t, a[8], c[8]);

                double[] d = new double[] {xD, yD, zD, nxD, nyD, nzD, rD, gD, bD};

                // SETS THE TWO TRAPEZOIDS FOR EACH TRIANGLE
                double[][] trapZero = tmpTrapezoids[2 * i];
                double[][] trapOne = tmpTrapezoids[(2 * i) + 1];

                copy2DVertex(a, trapZero[0]);
                copy2DVertex(a, trapZero[1]);

                copy2DVertex(c, trapOne[2]);
                copy2DVertex(c, trapOne[3]);

                if (b[0] < d[0]) {
                    // FIRST SCENARIO FROM NOTES
                    copy2DVertex(b, trapZero[2]);
                    copy2DVertex(d, trapZero[3]);

                    copy2DVertex(b, trapOne[0]);
                    copy2DVertex(d, trapOne[1]);

                } else {
                    // SECOND SCENARIO FROM NOTES
                    copy2DVertex(d, trapZero[2]);
                    copy2DVertex(b, trapZero[3]);

                    copy2DVertex(d, trapOne[0]);
                    copy2DVertex(b, trapOne[1]);
                }
            }
        }
    }

    // SCANS THE TRAPEZOIDS AND SETS THE PIXELS IN THE TEMP PIXEL ARRAY
    private void scanTrapezoids() {

        for (int i = 0; i < tmpTrapezoids.length; i++) {

            // ONLY SCAN TRAPEZOIDS IF TRIANGLE IS POSITIVE
            if (positive[i / 2]) {
                double[] LT = tmpTrapezoids[i][0];
                double[] RT = tmpTrapezoids[i][1];
                double[] LB = tmpTrapezoids[i][2];
                double[] RB = tmpTrapezoids[i][3];

                double XLT = LT[0];
                double XRT = RT[0];
                double XLB = LB[0];
                double XRB = RB[0];

                double YT = LT[1];
                double YB = LB[1];

                // FOR EACH SCANLINE
                for (int y = (int) YT; y < (int) YB; y++) {

                    double yPercent = (y - YT) / (YB - YT);
                    double XL = (XLT + yPercent * (XLB - XLT));
                    double XR = (XRT + yPercent * (XRB - XRT));

                    double pzXL = interpolate(yPercent, LT[2], LB[2]);
                    double rXL = interpolate(yPercent, LT[6], LB[6]);
                    double gXL = interpolate(yPercent, LT[7], LB[7]);
                    double bXL = interpolate(yPercent, LT[8], LB[8]);

                    double pzXR = interpolate(yPercent, RT[2], RB[2]);
                    double rXR = interpolate(yPercent, RT[6], RB[6]);
                    double gXR = interpolate(yPercent, RT[7], RB[7]);
                    double bXR = interpolate(yPercent, RT[8], RB[8]);

                    for (int x = (int) XL; x <= (int) XR; x++) {

                        double xPercent = (x - XL) / (XR - XL);

                        // ONLY DRAW PIXEL IF IT'S IN THE BOUNDS OF THE WINDOW
                        if ((x >= 0 && x < w) && (y >= 0 && y < h)) {

                            double pz = interpolate(xPercent, pzXL, pzXR);

                            if (pz > zbuf[x][y]) {
                                pix[x + w * y][0] = (int) interpolate(xPercent, rXL, rXR);
                                pix[x + w * y][1] = (int) interpolate(xPercent, gXL, gXR);
                                pix[x + w * y][2] = (int) interpolate(xPercent, bXL, bXR);

                                zbuf[x][y] = pz;
                            }
                        }
                    }
                }
            }
        }
    }

    // KEN PERLIN ORIGINALLY WROTE THIS, MODIFIED FOR Z BUFFER STUFF
    private void projectPoint(double[] xyz, double[] pxyz) {
        double x = xyz[0];
        double y = xyz[1];
        double z = xyz[2];

        // SET PROJECTED X, Y, AND Z
        pxyz[0] = w / 2 + (int)(h * x / (FL - z));
        pxyz[1] = h / 2 - (int)(h * y / (FL - z));
        pxyz[2] = (FL * z) / (FL - z);

        pxyz[3] = xyz[3];
        pxyz[4] = xyz[4];
        pxyz[5] = xyz[5];
    }

    // APPLY LIGHTING CALCULATIONS TO VERTEX USING ITS MATERIAL
    private void lightPoint(double[] vertex, Material material) {

        // NORMAL AT THE VERTEX
        double[] normal = {vertex[3], vertex[4], vertex[5]};

        // RESET TEMPORARY COLOR ARRAY
        for (int i = 0; i < tmpCol.length; i++) {
            tmpCol[i] = 0.0;
        }

        // FOR EACH LIGHT SOURCE
        for (int i = 0; i < lights.length; i++) {

            // LIGHT DIRECTION
            double[] lDir = lights[i][0];

            // REFLECTION DIRECTION
            for (int k = 0; k < lDir.length; k++) {
                tmpVec[k] = 2 * (dot(lDir, normal)) * normal[k] - lDir[k];
            }
            double[] rDir = tmpVec;
            normalize(rDir);

            // FOR EACH COLOR IN RGB
            for (int j = 0; j < 3; j++) {
                tmpCol[j] += tmpMaterial.getAmbientColor()[j] + lights[i][1][j] *
                            (tmpMaterial.getDiffuseColor()[j] * Math.max(0.0, dot(lDir, normal)) +
                             tmpMaterial.getSpecularColor()[j] * Math.pow(Math.max(0.0, dot(rDir, eyeDir)), tmpMaterial.getSpecularPower()));
            }
        }

        // GAMMA CORRECTION

        if (!showNormals) {
            vertex[6] = 255.0 * Math.pow(tmpCol[0], 0.45);
            vertex[7] = 255.0 * Math.pow(tmpCol[1], 0.45);
            vertex[8] = 255.0 * Math.pow(tmpCol[2], 0.45);
        } else {
            vertex[6] = mapRGB(normal[0]);
            vertex[7] = mapRGB(normal[1]);
            vertex[8] = mapRGB(normal[2]);
        }
    }

    private void normalize(double[] vector) {
        double length = Math.sqrt(vector[0] * vector[0] +
                                  vector[1] * vector[1] +
                                  vector[2] * vector[2]);
        vector[0] /= length;
        vector[1] /= length;
        vector[2] /= length;
    }

    private double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    // MAP COMPONENT OF SURFACE NORMAL TO COMPONENT OF RGB COLOR
    private int mapRGB(double normal) {
        return (int) ((normal + 1.0) * 255 / 2.0);
    }

    private double interpolate(double percent, double one, double two) {
        return percent * (two - one) + one;
    }

    private void copy2DVertex(double[] src, double[] dest) {
        dest[0] = src[0];
        dest[1] = src[1];
        dest[2] = src[2];
        dest[3] = src[3];
        dest[4] = src[4];
        dest[5] = src[5];
        dest[6] = src[6];
        dest[7] = src[7];
        dest[8] = src[8];
    }

    private void set2DVertex(double[] dest, double x, double y) {
        dest[0] = x;
        dest[1] = y;
    }

    public void setFocalLength(double l) {
        this.FL = l;
    }

    public void add(Geometry geo) {
        world.add(geo);
    }

    public void remove(Geometry geo) {
        world.remove(geo);
    }

    public void setBGColor(int[] rgb) {
        this.tmpBgColor = rgb;
    }

    public void setDrawColor(int[] rgb) {
        this.tmpColor = rgb;
    }

    public void showNormals(boolean show) {
        this.showNormals = show;
    }
}