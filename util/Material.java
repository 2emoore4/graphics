package util;

public class Material {

    double[] ambientColor, diffuseColor, specularColor;
    double specularPower;

    public Material() {
        ambientColor = new double[3];
        diffuseColor = new double[3];
        specularColor = new double[3];
    }

    public void setAmbientColor(double r, double g, double b) {
        ambientColor[0] = r;
        ambientColor[1] = g;
        ambientColor[2] = b;
    }

    public double[] getAmbientColor() {
        return ambientColor;
    }

    public void setDiffuseColor(double r, double g, double b) {
        diffuseColor[0] = r;
        diffuseColor[1] = g;
        diffuseColor[2] = b;
    }

    public double[] getDiffuseColor() {
        return diffuseColor;
    }

    public void setSpecularColor(double r, double g, double b) {
        specularColor[0] = r;
        specularColor[1] = g;
        specularColor[2] = b;
    }

    public double[] getSpecularColor() {
        return specularColor;
    }

    public void setSpecularPower(double p) {
        specularPower = p;
    }

    public double getSpecularPower() {
        return specularPower;
    }
}