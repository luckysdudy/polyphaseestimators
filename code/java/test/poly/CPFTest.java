/*
 */
package poly;

import org.junit.After;
import org.junit.AfterClass;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import pubsim.Complex;
import static pubsim.VectorFunctions.print;
import org.mckilliam.distributions.Gaussian;

/**
 *
 * @author Robby McKilliam
 */
public class CPFTest {
    
    public CPFTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }


    /**
     * Test of z method, of class CPF.
     */
    @Test
    public void testZ() {
        System.out.println("test z");
        double tol = 1e-8;
        int n = 3;
        CPF instance = new CPF(n);
        double[] real = {1.0,2.0,3.0};
        double[] imag = {1.0,2.0,3.0};
        instance.estimate(real, imag);
        assertTrue( (instance.z(-2).minus(new Complex(0,0))).abs() < tol );
        assertTrue( (instance.z(-1).minus(new Complex(1,1))).abs() < tol );
        assertTrue( (instance.z(0).minus(new Complex(2,2))).abs() < tol );
        assertTrue( (instance.z(1).minus(new Complex(3,3))).abs() < tol );
        assertTrue( (instance.z(2).minus(new Complex(0,0))).abs() < tol );
    }
    
     /**
     * Test of estimate method, of class DPTEstimator.
     */
    @Test
    public void testCPmax() {
        System.out.println("testHighestOrderParameter");
        
        double tol = 1e-5;
        
        int N = 257;
        int m = 3;
        double[] oparams = {1.0, Math.PI/8, 0.005, 0.00001};
        CPF inst = new CPF(N,10*N);
        double[] params = inst.transformToStandardBasis(oparams);
        System.out.println(print(params));
        
        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(N, new Gaussian(0, 0), params);
        siggen.generateReceivedSignal();
        inst.estimate(siggen.real(), siggen.imag());
        
        //check that maxCP(0) is twice the quadratic parameter
        assertEquals(inst.maxCP(0)/2, oparams[2], tol); 

    }
    
    /**
     * Test estimate of quadratic and cubic parameters
     */
    @Test
    public void testQuadraticAndCubic() {
        System.out.println("test quadratic and cubic");
        
        double tol = 1e-5;
        
        int N = 211;
        int m = 3;
        double[] oparams = {1.0, Math.PI/8, 0.005, 0.00001};
        CPF inst = new CPF(N,10*N);
        double[] params = inst.transformToStandardBasis(oparams);
        
        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(N, new Gaussian(0, 0), params);
        siggen.generateReceivedSignal();
        
        double[] phat = inst.estimate(siggen.real(), siggen.imag());
        
        assertEquals(params[3], phat[3], tol);
        assertEquals(params[2], phat[2], tol);

    }

}
