/*
 */
package poly;

import Jama.Matrix;
import pubsim.VectorFunctions;
import org.mckilliam.lattices.ClosestVectorInterface;
import org.mckilliam.lattices.Vnmstar.HilbertMatrix;
import org.mckilliam.lattices.Vnmstar.VnmStar;
import org.mckilliam.lattices.reduction.LLL;
import org.mckilliam.lattices.reduction.LatticeReduction;

/**
 * Uses the Babai nearest plane algorithm
 * @author Robby McKilliam
 */
public class Babai extends AbstractPolynomialPhaseEstimator {

    final protected double[] ya,  p;
    final protected int n;
    final protected VnmStar lattice;
    protected ClosestVectorInterface npalgorithm;
    final protected Matrix M,  K;
    
    public Babai(int m, int n){
        this(m,n,new LLL());
    } 
    
    /** 
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Babai(int m, int n, LatticeReduction lr) {
        super(m);
        lattice = new VnmStar(m, n - m - 1);
        npalgorithm = new org.mckilliam.lattices.cvp.Babai(lattice, lr);
        ya = new double[n];
        p = new double[m+1];
        this.n = n;
        M = lattice.getMMatrix();
        //System.out.println(Mt.times(M).inverse().cond());
        K = new HilbertMatrix(m+1,n).KDouble();
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length) throw new RuntimeException("Data length does not equal " + n);
        
        for (int i = 0; i < real.length; i++) {
            ya[i] = Math.atan2(imag[i], real[i]) / (2 * Math.PI);
        }
        npalgorithm.closestPoint(ya);
        double[] u = npalgorithm.getIndex();

        double[] ymu = new double[ya.length];
        for (int i = 0; i < u.length; i++) {
            ymu[i] = ya[i] - u[i];
        }
        System.arraycopy(ya, u.length, ymu, u.length, ya.length - u.length);

        //compute the parameters
        VectorFunctions.matrixMultVector(K, ymu, p); 
        
        return ambiguityRemover.disambiguate(p);
    }

}
