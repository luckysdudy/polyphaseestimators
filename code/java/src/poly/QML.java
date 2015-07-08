package poly;

import Jama.Matrix;
import pubsim.Complex;

/**
 * This is the Quasi maximum likelihood (QML) estimator from
 * 
 * Djurovic, I., Simeunovic, M., "Resolving aliasing effect in the QML estimation of PPSs"
 * 
 * Makes an estimate of the instantaneous frequency using the short term Fourier
 * transform and applies phase unwrapping to this.
 * 
 * @author Robby McKilliam
 */
public class QML extends AbstractPolynomialPhaseEstimator {
    
    public final int N;
    public final double h;

    protected final double[] wh; //this stores the maxima of the short-term Fourier transform
    
    public QML(int N, int m, double h){
        super(m);
        this.N = N;
        this.h = h;
        wh = new double[N];
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        
        //compute maximum of short term Fourier transform at each sample.
        
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    /**
     * Returns the maximiser of the absolute value of the short term Fourier transform at time t
     * with window size h.
     */
    public static double max_stft(double t, double h, Complex[] x) {
        double fstep = 1.0/h/4.0; // 4 times oversampling for obtain course maximiser
        double cbest = Double.NEGATIVE_INFINITY;
        double cf = -0.5;
        for(double f = -0.5; f < 0.5; f+=fstep){
            double stftabs = stft(t,f,h,x).abs2();
            if(stftabs > cbest) {
                cbest = stftabs;
                cf = f;
            }
        }
        return cf;  //could apply Brent's method or otherwise to optimise this, but's it's unlikely to be worth it
    }
    
    /** 
     * Short term Fourier transform of sequence x at time t, frequency f, and window size h.
     * The sample rate is assumed to be 1.
     */
    public static Complex stft(double t, double f, double h, Complex[] x){
        int N = x.length;
        int to = (int) Math.min(N-1, Math.floor(t + h/2));
        int from = (int) Math.max(0, Math.ceil(t - h/2));
        double w = f*2*Math.PI; //frequency in radians per second
        Complex sum = new Complex(0,0);
        for(int n = from; n <= to; n++) {
            Complex term = new Complex.UnitCircle(-w*n).multiply(x[n]);
            sum = sum.add(term);
        }
        return sum;
    }
    
    /**
     * Transformation matrix to standard polynomial basis from the origin
     * centered basis used in the papers of Djurovic
     * @param N
     * @param m
     * @return 
     */
    public static Matrix QMLcenteredtransformation(int N, int m) {
        double k = -N/2.0 - 1.0;
        Matrix T = new Matrix(m+1,m+1);
        for(int t = 0; t <=m; t++) {
            for(int i = t; i <= m; i++) {
                T.set(t,i, pubsim.Util.binom(i, t) * Math.pow(k, i-t) / 2.0 / Math.PI );          
            }
        }
        return T;
    }
    
    // Transform from origin centered basis to standard basis. 
    public static double[] toStandardBasis(double[] a, int N, int m){
        Matrix B = QMLcenteredtransformation(N, m);
        return pubsim.VectorFunctions.matrixMultVector(B,a);
    }
    
}
