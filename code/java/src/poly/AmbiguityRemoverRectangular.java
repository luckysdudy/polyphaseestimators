/*
 */
package poly;

import org.mckilliam.lattices.Lattice;
import org.mckilliam.lattices.cvp.Babai;
import org.mckilliam.lattices.reduction.None;

/**
 * Disambiguates polynomial phase signals into the rectangular region in the paper
 * @author Robby McKilliam
 */
public class AmbiguityRemoverRectangular extends AmbiguityRemover {

    public AmbiguityRemoverRectangular(int m) {
        this.m = m;
        p = new double[m+1];
        M = constructBasisMatrix();
        Lattice lattice = new Lattice(M);
        np = new Babai(lattice, new None());
    }
    
}
