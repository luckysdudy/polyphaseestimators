package poly;

import org.junit.*;
import static org.junit.Assert.*;
import pubsim.VectorFunctions;
import org.mckilliam.distributions.circular.WrappedUniform;

/**
 *
 * @author mckillrg
 */
public class CircularNoisePolynomialPhaseSignalTest {

    public CircularNoisePolynomialPhaseSignalTest() {
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
     * Test of generateReceivedSignal method, of class
     * CircularNoisePolynomialPhaseSignal.
     */
    @Test
    public void testGenerateReceivedSignal() {
        System.out.println("generateReceivedSignal");
        int N = 5;
        double[] p = {0.0, 0.0, 0.0};
        
        CircularNoisePolynomialPhaseSignal instance = new CircularNoisePolynomialPhaseSignal(N, new WrappedUniform(0, 0.01), p);

        instance.generateReceivedSignal();
        System.out.println(VectorFunctions.print(instance.real()));
        System.out.println(VectorFunctions.print(instance.imag()));
        
        instance.generateReceivedSignal();
        System.out.println(VectorFunctions.print(instance.real()));
        System.out.println(VectorFunctions.print(instance.imag()));
    }
}
