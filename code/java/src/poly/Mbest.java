package poly;

import org.mckilliam.lattices.reduction.LLL;
import org.mckilliam.lattices.reduction.LatticeReduction;

/**
 * @author Robby McKilliam
 */
public class Mbest extends Babai {

    /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Mbest(int m, int n, int M) {
        this(m,n,M,new LLL());
    }
    
     /**
     * You must set the polynomial order in the constructor
     * @param m = polynomial order
     */
    public Mbest(int m, int n, int M, LatticeReduction lr) {
        super(m,n);
        npalgorithm = new org.mckilliam.lattices.cvp.Mbest(lattice,M,lr);
    }
}