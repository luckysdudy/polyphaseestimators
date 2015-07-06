package poly;

import poly.PolynomialPhaseSignal;
import junit.framework.TestCase;
import org.mckilliam.distributions.Gaussian;
import pubsim.VectorFunctions;

/**
 *
 * @author Robby
 */
public class PolynomialPhaseSignalTest extends TestCase {
    
    public PolynomialPhaseSignalTest(String testName) {
        super(testName);
    }            

    /**
     * Test of generateReceivedSignal method, of class PolynomialPhaseSignal.
     */
    public void testGenerateReceivedSignal() {
        System.out.println("generateReceivedSignal");
        
        int n = 4;
        double[] params = {1, 0.1, 0.2};
        
        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(n, new Gaussian(0, 0.0), params);
        
        //these are taken from matlab
        double[] expreal = {-0.3090, 1.0000, 0.8090, -0.8090};
        double[] expimag = {0.9511, -0.0000, 0.5878, -0.5878};
        
        siggen.generateReceivedSignal();
        
        System.out.println("real = " + VectorFunctions.print(siggen.real()));
        System.out.println("imag = " + VectorFunctions.print(siggen.imag()));
        
        assertEquals(VectorFunctions.distance_between(expreal, siggen.real())<0.001, true);
        assertEquals(VectorFunctions.distance_between(expimag, siggen.imag())<0.001, true);
        
    }

        /**
     * Test of generateReceivedSignal method, of class PolynomialPhaseSignal.
     */
    public void testGenerateRandomParameters() {
        System.out.println("testGenerateRandomParameters");
        int m = 2;
        PolynomialPhaseSignal.RandomParameterGenerator pgen
                = new PolynomialPhaseSignal.RandomParameterGenerator(m);

        double[] p = pgen.generateParameters();
        assertEquals(3, p.length);
        System.out.println(VectorFunctions.print(p));

    }

}
