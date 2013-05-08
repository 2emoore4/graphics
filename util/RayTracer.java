package util;

import java.lang.Math.*;
import java.util.ArrayList;

public class RayTracer {

    int width, height, xRes, yRes;
    int[] tmpBgColor, tmpBgColor2, tmpIntColor;
    int[][] pix, antialiasBuf;
    double FL;
    double[] v, w, s, v_s, nn, t, tempColor, tempVector;
    double[][] zbuf;
    ArrayList<RayShape> sceneChildren;
    Material tempMaterial;

    boolean ss = true;
    boolean aa = false;
    boolean textures = false;
    int type = 1;

    int subSampleRate = 2;

    // ARRAY OF LIGHTS
    double[][][] lights = {
        { { 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0} },
        { {-0.25,-0.25,-0.25}, {1.0, 1.0, 1.0} },
        { {0.0, 0.0, 1.0},  {1.0, 1.0, 1.0} },
    };
    double[] eyeDir = {0.0, 0.0, 1.0};

    public RayTracer(int w, int h) {

        this.width = w;
        this.height = h;
        this.xRes = 1;
        this.yRes = 1;

        this.tmpBgColor2 = new int[3];
        this.tmpIntColor = new int[3];

        this.pix = new int[width * height][3];
        this.antialiasBuf = new int[17][3];

        this.FL = 10;

        this.v = new double[3];
        this.w = new double[3];
        this.s = new double[4];
        this.v_s = new double[3];
        this.nn = new double[3];
        this.t = new double[2];
        this.tempColor = new double[3];
        this.tempVector = new double[3];

        this.zbuf = new double[width][height];

        this.sceneChildren = new ArrayList<RayShape>();
    }

    // USED BY MISAPPLET TO GET THE RGB FOR A SPECIFIC PIXEL
    public int[] getPix(int x, int y) {
        return pix[x + width * y];
    }

    // RESETS THE PIXEL ARRAY TO THE BACKGROUND COLOR
    public void resetPix() {

        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                pix[x + width * y][0] = (int) interpolate(y / (double) height, tmpBgColor[0], tmpBgColor2[0]);
                pix[x + width * y][1] = (int) interpolate(y / (double) height, tmpBgColor[1], tmpBgColor2[1]);
                pix[x + width * y][2] = (int) interpolate(y / (double) height, tmpBgColor[2], tmpBgColor2[2]);
            }
        }

        for (int i = 0; i < zbuf.length; i++) {
            for (int j = 0; j < zbuf[0].length; j++) {
                zbuf[i][j] = Double.MAX_VALUE;
            }
        }
    }

    // MAIN LOOP
    public void renderWorld() {

        castRay(0, 0);

        for (int y = subSampleRate; y < height; y += subSampleRate) {

            castRay(0, y);

            for (int x = subSampleRate; x < width; x += subSampleRate) {

                if (ss) {
                    if (edgeDetect(x, y, subSampleRate)) {
                        // EDGE DETECTED. CAST MORE RAYS.
                        for (int ySub = y - subSampleRate; ySub < y; ySub++) {
                            for (int xSub = x - subSampleRate; xSub < x; xSub++) {
                                castRay(xSub, ySub);
                            }
                        }
                    } else {
                        // NO EDGE DETECTED. JUST INTERPOLATE.
                        castRay(x - subSampleRate, y - subSampleRate);

                        int[] q12 = pix[(x - subSampleRate) + width * y];
                        int[] q21 = pix[x + width * (y - subSampleRate)];
                        int[] q11 = pix[(x - subSampleRate) + width * (y - subSampleRate)];
                        int[] q22 = pix[x + width * y];

                        for (int ySub = y - subSampleRate; ySub < y; ySub++) {
                            double upPercent = (y - ySub) / (double) subSampleRate;
                            for (int xSub = x - subSampleRate; xSub < x; xSub++) {
                                double leftPercent = (x - xSub) / (double) subSampleRate;

                                pix[xSub + width * ySub][0] = (int) interpolate(q11[0], q12[0], q21[0], q22[0], xSub, x - subSampleRate, x, ySub, y - subSampleRate, y);
                                pix[xSub + width * ySub][1] = (int) interpolate(q11[1], q12[1], q21[1], q22[1], xSub, x - subSampleRate, x, ySub, y - subSampleRate, y);
                                pix[xSub + width * ySub][2] = (int) interpolate(q11[2], q12[2], q21[2], q22[2], xSub, x - subSampleRate, x, ySub, y - subSampleRate, y);
                            }
                        }
                    }
                } else {
                    for (int ySub = y - subSampleRate; ySub < y; ySub++) {
                        for (int xSub = x - subSampleRate; xSub < x; xSub++) {
                            castRay(xSub, ySub);
                        }
                    }
                }
            }
        }
    }

    // CASTS A RAY AT (I,J)
    private void castRay(int i, int j) {

        double x = 0.8 * ((i + 0.5) / width - 0.5);
        double y = 0.8 * ((j + 0.5) / width - 0.5 * height / width);

        set(v, 0, 0, FL);
        set(w, x, y, -FL);

        normalize(w);

        for (int it = 0; it < sceneChildren.size(); it++) {
            RaySphere sphere = (RaySphere) sceneChildren.get(it);
            tempMaterial = sphere.getMaterial();

            s[0] = sphere.getPosition()[0];
            s[1] = sphere.getPosition()[1];
            s[2] = sphere.getPosition()[2];
            s[3] = sphere.getRadius();

            if (rayTrace(v, w, t)) {

                if (t[0] < zbuf[i][j]) {

                    for (int k = 0; k < 3; k++) {
                        nn[k] = v[k] + t[0] * w[k] - s[k];
                    }
                    normalize(nn);
                    computeShading(nn);

                    pix[i + width * j][0] = (int) tempColor[0];
                    pix[i + width * j][1] = (int) tempColor[1];
                    pix[i + width * j][2] = (int) tempColor[2];

                    zbuf[i][j] = t[0];
                }
            }
        }
    }

    // RETURNS TRUE IF AN EDGE IS DETECTED IN THE 4X4 PX SQUARE WITH X,Y AS THE BOTTOM RIGHT CORNER
    private boolean edgeDetect(int x, int y, int diff) {

        castRay(x, y);

        int[] left = pix[(x - diff) + width * y];
        int[] up = pix[x + width * (y - diff)];
        int[] upLeft = pix[(x - diff) + width * (y - diff)];
        int[] current = pix[x + width * y];

        int leftDiff = (current[0] - left[0]) * (current[0] - left[0]) + (current[1] - left[1]) * (current[1] - left[1]) + (current[2] - left[2]) * (current[2] - left[2]);
        int upDiff = (current[0] - up[0]) * (current[0] - up[0]) + (current[1] - up[1]) * (current[1] - up[1]) + (current[2] - up[2]) * (current[2] - up[2]);
        int upLeftDiff = (current[0] - upLeft[0]) * (current[0] - upLeft[0]) + (current[1] - upLeft[1]) * (current[1] - upLeft[1]) + (current[2] - upLeft[2]) * (current[2] - upLeft[2]);

        if (leftDiff > 5000 || upDiff > 5000 || upLeftDiff > 5000) {
            return true;
        } else {
            return false;
        }
    }

    // DETERMINE IF A GIVEN RAY INTERSECTS A GIVEN SPHERE
    private boolean rayTrace(double[] v, double[] w, double[] t) {
        diff(v, s, v_s);

        double a = 1.0;
        double b = 2 * dot(w, v_s);
        double c = dot(v_s, v_s) - s[3] * s[3];

        return solveQuadraticEquation(a, b, c, t);
    }

    // DOES WHAT IT SAYS
    private boolean solveQuadraticEquation(double a, double b, double c, double[] t) {
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) {
            return false;
        }

        double d = Math.sqrt(discriminant);
        t[0] = (-b - d) / (2 * a);
        t[1] = (-b + d) / (2 * a);
        return true;
    }

    // COMPUTE SHADING FOR INDIVIDUAL PIXEL, DEPENDING ON NORMAL DIRECTION
    private void computeShading(double[] nn) {

        // RESET TEMPORARY COLOR ARRAY
        for (int i = 0; i < tempColor.length; i++) {
            tempColor[i] = tempMaterial.getAmbientColor()[i];
        }

        if (textures) {
            double f0 = f(nn[0]      ,nn[1]      ,nn[2]      ),
                   fx = f(nn[0]+.0001,nn[1]      ,nn[2]      ),
                   fy = f(nn[0]      ,nn[1]+.0001,nn[2]      ),
                   fz = f(nn[0]      ,nn[1]      ,nn[2]+.0001);

            // SUBTRACT THE FUNCTION'S GRADIENT FROM THE SURFACE NORMAL

            nn[0] -= (fx - f0) / .0001;
            nn[1] -= (fy - f0) / .0001;
            nn[2] -= (fz - f0) / .0001;
        }

        // FOR EACH LIGHT SOURCE
        for (int i = 0; i < lights.length; i++) {

            // LIGHT DIRECTION
            double[] lDir = lights[i][0];

            // REFLECTION DIRECTION
            for (int k = 0; k < 3; k++) {
                tempVector[k] = 2 * (dot(lDir, nn)) * nn[k] - lDir[k];
            }
            double[] rDir = tempVector;
            normalize(rDir);

            // FOR EACH COLOR IN RGB
            for (int j = 0; j < 3; j++) {
                tempColor[j] += lights[i][1][j] * (tempMaterial.getDiffuseColor()[j] * Math.max(0.0, dot(lDir, nn)) +
                                tempMaterial.getSpecularColor()[j] * Math.pow(Math.max(0.0, dot(rDir, eyeDir)), tempMaterial.getSpecularPower()));
            }
        }

        // GAMMA CORRECTION
        tempColor[0] = 255.0 * Math.pow(tempColor[0], 0.45);
        tempColor[1] = 255.0 * Math.pow(tempColor[1], 0.45);
        tempColor[2] = 255.0 * Math.pow(tempColor[2], 0.45);
    }

    // SETS VALUES FOR A SINGLE VECTOR
    private void set(double[] vector, double x, double y, double z) {
        vector[0] = x;
        vector[1] = y;
        vector[2] = z;
    }

    // NORMALIZES A VECTOR
    private void normalize(double[] vector) {
        double length = Math.sqrt(vector[0] * vector[0] +
                                  vector[1] * vector[1] +
                                  vector[2] * vector[2]);
        vector[0] /= length;
        vector[1] /= length;
        vector[2] /= length;
    }

    // SUBTRACTS V2 FROM V1 AND SETS TO V3
    private void diff(double[] v1, double[] v2, double[] v3) {
        v3[0] = v1[0] - v2[0];
        v3[1] = v1[1] - v2[1];
        v3[2] = v1[2] - v2[2];
    }

    // RETURNS THE DOT PRODUCT OF A AND B
    private double dot(double[] a, double[] b) {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    // ADDS A SHAPE TO THE SCENE
    public void add(RayShape shape) {
        sceneChildren.add(shape);
    }

    // REMOVES A SHAPE FROM THE SCENE
    public void remove(RayShape shape) {
        sceneChildren.remove(shape);
    }

    // SETS THE BACKGROUND COLOR, AS WELL AS THE GRADIENT COLOR
    public void setBGColor(int[] rgb) {
        this.tmpBgColor = rgb;

        for (int i = 0; i < rgb.length; i++) {
            this.tmpBgColor2[i] = this.tmpBgColor[i] + 50;
        }
    }

    // RETURNS THE ARRAY OF LIGHTS IN THE SCENE
    public double[][][] getLights() {
        return lights;
    }

    // SETS THE X RESOLUTION OF THE SCENE
    public void setXRes(int x) {
        this.xRes = x;
    }

    // SETS THE Y RESOLUTION OF THE SCENE
    public void setYRes(int y) {
        this.yRes = y;
    }

    // INTERPOLATE BETWEEN TWO POINTS
    private double interpolate(double percent, double one, double two) {
        return percent * (two - one) + one;
    }

    // INTERPOLATE BETWEEN THREE POINTS
    private double interpolate(int q11, int q12, int q21, int q22, int x, int x1, int x2, int y, int y1, int y2) {
        return (1 / (double) ((x2 - x1) * (y2 - y1))) * ((q11 * (x2 - x) * (y2 - y)) +
                                                         (q21 * (x - x1) * (y2 - y)) +
                                                         (q12 * (x2 - x) * (y - y1)) +
                                                         (q22 * (x - x1) * (y - y1)));
    }

    double f(double x, double y, double z) {
        switch (type) {
        case 0:  return  .06 * noise(x/2.0,y/2.0,z/2.0, 8);
        case 1:  return  .01 * stripes(x + 2*turbulence(x,y,z,1), 1.6);
        default: return -.10 * turbulence(x,y,z, 1);
        }
    }

    // STRIPES TEXTURE (GOOD FOR MAKING MARBLE)

    double stripes(double x, double f) {
       double t = .5 + .5 * Math.sin(f * 2*Math.PI * x);
       return t * t - .5;
    }

    // TURBULENCE TEXTURE

    double turbulence(double x, double y, double z, double freq) {
       double t = -.5;
       for ( ; freq <= width/12 ; freq *= 2)
          t += Math.abs(noise(x,y,z,freq) / freq);
       return t;
    }

    // NOISE TEXTURE

    double noise(double x, double y, double z, double freq) {
       double x1, y1, z1;
       x1 = .707*x-.707*z;
       z1 = .707*x+.707*z;
       y1 = .707*x1+.707*y;
       x1 = .707*x1-.707*y;
       return ImprovedNoise.noise(freq*x1 + 100, freq*y1, freq*z1);
    }

    private int[] castIntermediateRay(double i, double j, int centerX, int centerY) {

        double x = 0.8 * ((i + 0.5) / width - 0.5);
        double y = 0.8 * ((j + 0.5) / width - 0.5 * height / width);

        set(v, 0, 0, FL);
        set(w, x, y, -FL);

        normalize(w);

        for (int it = 0; it < sceneChildren.size(); it++) {
            RaySphere sphere = (RaySphere) sceneChildren.get(it);
            tempMaterial = sphere.getMaterial();

            s[0] = sphere.getPosition()[0];
            s[1] = sphere.getPosition()[1];
            s[2] = sphere.getPosition()[2];
            s[3] = sphere.getRadius();

            if (rayTrace(v, w, t)) {

                if (t[0] < zbuf[centerX][centerY]) {

                    for (int k = 0; k < 3; k++) {
                        nn[k] = v[k] + t[0] * w[k] - s[k];
                    }
                    normalize(nn);
                    computeShading(nn);

                    for (int k = 0; k < 3; k++) {
                        tmpIntColor[k] = (int) tempColor[k];
                    }

                    // zbuf[centerX][centerY] = t[0];

                    return tmpIntColor;
                }
            }
        }

        return null;
    }

    private void antialias(int i, int j) {
        double x = 0.8 * ((i + 0.5) / width - 0.5);
        double y = 0.8 * ((j + 0.5) / width - 0.5 * height / width);

        set(v, 0, 0, FL);
        set(w, x, y, -FL);

        normalize(w);

        for (int it = 0; it < sceneChildren.size(); it++) {
            RaySphere sphere = (RaySphere) sceneChildren.get(it);
            tempMaterial = sphere.getMaterial();

            s[0] = sphere.getPosition()[0];
            s[1] = sphere.getPosition()[1];
            s[2] = sphere.getPosition()[2];
            s[3] = sphere.getRadius();

            if (rayTrace(v, w, t)) {

                if (t[0] < zbuf[i][j]) {

                    for (int k = 0; k < 3; k++) {
                        nn[k] = v[k] + t[0] * w[k] - s[k];
                    }
                    normalize(nn);
                    computeShading(nn);

                    pix[i + width * j][0] = (int) tempColor[0];
                    pix[i + width * j][1] = (int) tempColor[1];
                    pix[i + width * j][2] = (int) tempColor[2];

                    zbuf[i][j] = t[0];
                }
            }
        }

        // antialiasBuf[0] = castIntermediateRay(x, y, x, y);
        // antialiasBuf[1] = castIntermediateRay(x - 0.25, y - 0.25, x, y);
        // antialiasBuf[2] = castIntermediateRay(x + 0.25, y - 0.25, x, y);
        // antialiasBuf[3] = castIntermediateRay(x - 0.25, y + 0.25, x, y);
        // antialiasBuf[4] = castIntermediateRay(x + 0.25, y + 0.25, x, y);
        // antialiasBuf[5] = castIntermediateRay(x - 0.5, y - 0.5, x, y);
        // antialiasBuf[6] = castIntermediateRay(x + 0.5, y - 0.5, x, y);
        // antialiasBuf[7] = castIntermediateRay(x - 0.5, y + 0.5, x, y);
        // antialiasBuf[8] = castIntermediateRay(x + 0.5, y + 0.5, x, y);
        // antialiasBuf[9] = castIntermediateRay(x - 0.75, y - 0.75, x, y);
        // antialiasBuf[10] = castIntermediateRay(x + 0.75, y - 0.75, x, y);
        // antialiasBuf[11] = castIntermediateRay(x - 0.75, y + 0.75, x, y);
        // antialiasBuf[12] = castIntermediateRay(x + 0.75, y + 0.75, x, y);
        // antialiasBuf[13] = castIntermediateRay(x - 1.0, y - 1.0, x, y);
        // antialiasBuf[14] = castIntermediateRay(x + 1.0, y - 1.0, x, y);
        // antialiasBuf[15] = castIntermediateRay(x - 1.0, y + 1.0, x, y);
        // antialiasBuf[16] = castIntermediateRay(x + 1.0, y + 1.0, x, y);

        // int rTotal = 0;
        // int gTotal = 0;
        // int bTotal = 0;
        // int num = 0;

        // for (int i = 0; i < antialiasBuf.length; i++) {
        //     if (antialiasBuf[i] != null) {
        //         rTotal += antialiasBuf[i][0];
        //         gTotal += antialiasBuf[i][1];
        //         bTotal += antialiasBuf[i][2];
        //     } else {
        //         rTotal += (int) interpolate(y / (double) height, tmpBgColor[0], tmpBgColor2[0]);
        //         gTotal += (int) interpolate(y / (double) height, tmpBgColor[1], tmpBgColor2[1]);
        //         bTotal += (int) interpolate(y / (double) height, tmpBgColor[2], tmpBgColor2[2]);
        //     }
        // }

        // tempColor[0] = (double) rTotal / (double) antialiasBuf.length;
        // tempColor[1] = (double) gTotal / (double) antialiasBuf.length;
        // tempColor[2] = (double) bTotal / (double) antialiasBuf.length;

        // pix[x + width * y][0] = (int) tempColor[0];
        // pix[x + width * y][1] = (int) tempColor[1];
        // pix[x + width * y][2] = (int) tempColor[2];

        // // pix[x + width * y][0] = 0;
        // // pix[x + width * y][1] = 0;
        // // pix[x + width * y][2] = 0;
    }

    // private int preciseEdgeLocation(int x, int y, int diff) {
    //     if (diff == 1) {
    //         // System.err.println("diff is 1");
    //         if (edgeDetect(x, y, 1)) {
    //             // System.err.println("edge is right");
    //             return x;
    //         } else if (edgeDetect(x - 1, y, 1)) {
    //             // System.err.println("edge is left");
    //             return x - 1;
    //         }
    //     } else {
    //         if (edgeDetect(x - diff / 2, y, diff / 2)) {
    //             // System.err.println("going left");
    //             return preciseEdgeLocation(x - diff / 2, y, diff / 4);
    //         } else if (edgeDetect(x, y, diff / 2)) {
    //             // System.err.println("going right");
    //             return preciseEdgeLocation(x, y, diff / 4);
    //         }
    //     }

    //     return -1;
    // }
}