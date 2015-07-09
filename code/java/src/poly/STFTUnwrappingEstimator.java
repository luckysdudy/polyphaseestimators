package poly;

import Jama.Matrix;
import org.mckilliam.lattices.Vnmstar.HilbertMatrix;
import org.mckilliam.optimisation.NewtonRaphson;
import pubsim.Complex;
import pubsim.VectorFunctions;

/**
 * This is part of the QML estimator from
 * 
 * Djurovic, I., Simeunovic, M., "Resolving aliasing effect in the QML estimation of PPSs"
 * 
 * Makes an estimate of the instantaneous frequency using the short term Fourier
 * transform and applies phase unwrapping to this.
 * 
 * The QML estimator runs this estimator for a range of different short term Fourier transform
 * widths.
 * 
 * @author Robby McKilliam
 */
public class STFTUnwrappingEstimator extends AbstractPolynomialPhaseEstimator {
    
    public final int N;
    public final double h;
    
    protected final Complex[] x;
    protected double[] afreq;
    protected double[] a; //array for storing returned polynomial phase parameters
    protected final double[] wh; //this stores the maxima of the short-term Fourier transform
    protected final double[] wh_unwrapped; //store unwrapped maxima
    protected final Matrix K; //projection matrix into the space of polynomials (for polynomial regression)
    protected final Matrix Kfreq; //projection matrix for polynomial regression on instantaneous frequency
    
    /**
     * @param m polynomial order
     * @param N number of samples
     * @param h width of the short term Fourier transform
     */
    public STFTUnwrappingEstimator(int m, int N, double h){
        super(m);
        this.N = N;
        this.h = h;
        wh = new double[N];
        x = new Complex[N];
        wh_unwrapped = new double[N];
        K = new HilbertMatrix(m+1,N).KDouble();
        Kfreq = new HilbertMatrix(m,N).KDouble();
        afreq = new double[m];
        a = new double[m+1];
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        for(int n = 0; n < N; n++) x[n] = new Complex(real[n], imag[n]);
        
        //compute maximum of short term Fourier transform at each sample.
        for(int n = 0; n < N; n++) wh[n] = max_stft(n,h,x); 
        
        //unwrap wh
        unwrap(wh, wh_unwrapped); 
        
        //polynomial regression on unwrapped instantaneous frequency
        VectorFunctions.matrixMultVector(Kfreq, wh_unwrapped, afreq); 
        for(int i = 0; i < m; i++) afreq[i] = afreq[i]/(i+1); //undo multipliers from derivative
        
        //dechirp signal to obtain phase estimate
        Complex sum = new Complex(0,0);
        for(int n = 0; n < N; n++){
            Complex dn = x[n];
            for(int i = 1; i <= m; i++) {
                Complex phase = new Complex.UnitCircle(-2*Math.PI*Math.pow(n+1, i)*afreq[i-1]);
                dn = phase.multiply(dn);
            }
            sum = sum.add(dn);
         }
        a[0] = sum.phase()/2/Math.PI;
        for(int i = 1; i <= m; i++) a[i] = afreq[i-1];
        
        //a = refine(a, real, imag);
        
        //mod back to identifiable region and return estimates.
        return ambiguityRemover.disambiguate(a);
    }

    /// Compute uwrapped x and return in y
    protected void unwrap(double[] x, double[] y) {
        y[0] = x[0];
        double c = 0;
        for(int n = 1; n < N; n++) {
            double delta = x[n] - x[n-1];
            if(delta > 0.5) c--;
            else if(delta < -0.5) c++;
            y[n] = x[n] + c;
        }
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
    
    /// Returns refined polynomial phase estimator.  Uses Newton-Raphson method
    protected double[] refine(double[] a, double[] real, double[] imag) throws ArithmeticException {
        MaximumLikelihood.PolynomialPhaseLikelihood func
                = new MaximumLikelihood.PolynomialPhaseLikelihood(real, imag);
        NewtonRaphson newtonRaphson
                = new NewtonRaphson(func);
        
        //refine the best parameter using Newton's method
        Matrix params = new Matrix(m+1,1);
        for(int i = 0; i <= m; i++) params.set(i,0,a[i]);
        try {
            Matrix p =  newtonRaphson.maximise(params);
            return VectorFunctions.unpackRowise(p);
        }catch(Exception e){
            //if Newton's method fails to converge just return the initial parameters.
            //This usually occurs when the initial guess is nowhere near the true parameters.
            return a; 
        }
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
