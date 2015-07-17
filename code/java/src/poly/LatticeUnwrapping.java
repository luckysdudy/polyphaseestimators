package poly;

import Jama.Matrix;
import org.mckilliam.lattices.cvp.ClosestVectorInterface;
import org.mckilliam.lattices.Vnmstar.HilbertMatrix;
import pubsim.VectorFunctions;

/**
 * General class for polynomial phase estimators bases on closest lattice point
 * algorithms.  Both MAP and LSU estimators can be implemented by setting the cvp
 * argument to the desired closest lattice point algorithm during construction.
 * @author Robby McKilliam
 */
public class LatticeUnwrapping extends AbstractPolynomialPhaseEstimator {
                
    /** Transformation from lattice basis to unwrapping variables */
    final Matrix T;
    
    /** Closest vector algorithm for this lattice */
    final protected ClosestVectorInterface cvp;
    
    /** Projection onto the space of polynomials */
    final protected Matrix K;
    
    /** Length of the polynomial phase signal */
    final protected int n;
    
    ///working memory
    final protected double[] ya,  p, w;

    public LatticeUnwrapping(int m, int n, ClosestVectorInterface cvp, Matrix T) {
        super(m);
        this.n = n;
        this.T = T;
        this.cvp = cvp;
        //System.out.println(Mt.times(M).inverse().cond());
        K = new HilbertMatrix(m+1,n).KDouble();
        ya = new double[n];
        w = new double[n];
        p = new double[m+1];
    }

    @Override
    public double[] estimate(double[] real, double[] imag) {
        if(n != real.length || n != imag.length) throw new RuntimeException("Data length does not equal " + n);
        for (int i = 0; i < n; i++) ya[i] = Math.atan2(imag[i], real[i]) / (2 * Math.PI); //phase of observations divided by 2pi
        cvp.closestPoint(ya); //compute closest lattice point
        double[] u = cvp.getIndex();
        VectorFunctions.matrixMultVector(T, u, w); //compute the unwrapping variables
        for (int i = 0; i < n; i++) ya[i] = ya[i] - w[i]; //subtract unwrapping variables
        VectorFunctions.matrixMultVector(K, ya, p);  //compute  parameters
        return ambiguityRemover.disambiguate(p); //return parameters wrapped back into identifiable region
    }
    
    
    
    
    
}
