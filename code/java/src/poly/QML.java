package poly;

import Jama.Matrix;

/**
 * This is part of the QML estimator from
 * 
 * Djurovic, I., Simeunovic, M., "Resolving aliasing effect in the QML estimation of PPSs"
 * 
 * Makes an estimate of the instantaneous frequency using the short term Fourier
 * transform and applies phase unwrapping to this.
 * 
 * Iterates of Short term Fourier transform widths and selects width based on the QML
 * optimisation criterion.
 * 
 * @author Robby McKilliam
 */
public class QML extends AbstractPolynomialPhaseEstimator {
    
    protected final STFTUnwrappingEstimator[] stftestimators; 
    protected final int[] H;
    
    public QML(int m, int N, int[] H){
        super(m);
        this.H = H;
        stftestimators = new STFTUnwrappingEstimator[H.length];
        for( int k = 0; k < H.length; k++ ) stftestimators[k] = new STFTUnwrappingEstimator(m,N,H[k]);
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        MaximumLikelihood.PolynomialPhaseLikelihood func = new MaximumLikelihood.PolynomialPhaseLikelihood(real, imag);
        double QMLbest = Double.NEGATIVE_INFINITY;
        double[] pbest = new double[m+1];
        Matrix params = new Matrix(m+1,1);
        for(int k = 0; k < H.length; k++){
            double[] p = stftestimators[k].estimate(real, imag);
            for(int i = 0; i <= m; i++) params.set(i,0,p[i]);
            double likelihood = func.value(params);
            if( likelihood > QMLbest ) {
                QMLbest = likelihood;
                pbest = p;
            }
        }
        pbest = MaximumLikelihood.refine(pbest, real, imag);
        return ambiguityRemover.disambiguate(pbest);
    }
    
    
}
