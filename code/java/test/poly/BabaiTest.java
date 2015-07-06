package poly;

import Jama.Matrix;
import org.mckilliam.distributions.Gaussian;
import org.mckilliam.lattices.Vnmstar.VnmStar;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.VectorFunctions;
import static org.junit.Assert.*;

/**
 *
 * @author Robby McKilliam
 */
public class BabaiTest {

    public BabaiTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

            /**
     * Test of estimate method, of class DPTEstimator.
     */
    @Test
    public void testEstimate() {
        System.out.println("testEstimate");

        int n = 24;
        double[] params = {0.11, 0.05002, 0.0205, 0.0001};
        int m = params.length-1;

        Matrix M = VnmStar.getGeneratorMatrix(m, n-m-1);
        System.out.println(VectorFunctions.print(M));

        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n,new Gaussian(0, 0.00001),params);

        siggen.generateReceivedSignal();
        PolynomialPhaseEstimatorInterface inst = new Babai(m,n);
        double[] p = inst.estimate(siggen.real(), siggen.imag());

        System.out.println(VectorFunctions.print(p));

        assertTrue(VectorFunctions.distance_between(p, params) < 0.001);

    }

}