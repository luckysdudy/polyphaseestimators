package poly;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import org.mckilliam.distributions.Gaussian;
import pubsim.Complex;

/**
 *
 * @author Robby McKilliam
 */
public class PeriodogramFrequencyEstimatorTest {

    public PeriodogramFrequencyEstimatorTest() {
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
     * Test of estimateFreq method, of class simulator.fes.PeriodogramEstimator.
     */
    @Test
    public void testEstimateFreq() {
        System.out.println("estimateFreq");

        final int N = 64;
        int iters = 100;
        double f = 0.1;
        double p = 0.4;
        
        PeriodogramFrequencyEstimator instance = new PeriodogramFrequencyEstimator(N);

        Gaussian noise = new Gaussian(0, 0.0001);
        double[] signalreal = new double[N];
        double[] signalimag = new double[N];
        
        for(int i=0; i < iters; i++ ){
            for(int n = 0; n < N; n++){
                signalreal[n] = Math.cos(n*2*Math.PI*f + p) + noise.noise();
                signalimag[n] = Math.sin(n*2*Math.PI*f + p) + noise.noise();
            }
            double result = instance.estimateFreq(signalreal, signalimag);
            System.out.print(result);
            assertEquals(true, Math.abs(result - f) < 0.02);
        }

    }

}