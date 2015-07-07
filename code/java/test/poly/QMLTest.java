package poly;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;
import pubsim.Complex;

/**
 *
 * @author Robby McKilliam
 */
public class QMLTest {
    
    public QMLTest() {
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
    
    public static double diff_complex(Complex x, Complex y){
        return (x.minus(y)).abs();
    }
        
    /**
     * Test of stft method, of class QML.
     */
    @Test
    public void testStft() {
        double tol = 1e-7; 
        System.out.println("stft");
        double h = 6.0;
        Complex[] x = {new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0)};
        assertTrue(diff_complex(new Complex(5,0), QML.stft(2.0, 0.0, h, x)) < tol);
        assertTrue(diff_complex(new Complex(0,0), QML.stft(2.0, 1.0/5, h, x)) < tol);
        assertTrue(diff_complex(new Complex(3,0), QML.stft(-0.1, 0.0, h, x)) < tol);
        assertTrue(diff_complex(new Complex(4,0), QML.stft(3.1, 0.0, h, x)) < tol);
    }
    
    /**
     * Test max of stft method, of class QML.
     */
    @Test
    public void testmaxStft() {
        double tol = 1e-7; 
        System.out.println("max stft");
        double h = 6.0;
        Complex[] x = {new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0)};
        assertTrue( Math.abs(QML.max_stft(2.0, h, x) - 0.0) < tol);
      }

    /**
     * Writes stft of a polynomial phase signal to file for the purpose of ploting (with gnuplot or otherwise)
     */
    @Test
    public void testWriteSTFTtoFile() {
        System.out.println("Write stft to file");
        
        Complex[] x = {new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0), new Complex(1,0)};
        assertTrue( Math.abs(QML.max_stft(2.0, h, x) - 0.0) < tol);
      }
    
}
