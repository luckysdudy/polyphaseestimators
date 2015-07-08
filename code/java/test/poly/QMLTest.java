package poly;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.logging.Level;
import java.util.logging.Logger;
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
        int N = 256;
        double h =12;
        double[] b = {1.0/10, -2.0/10, -1.0/200, 1.0/50000};
        //double noisevar = 1.0/2.0; //SNR = 1
        double noisevar = 0.0; //SNR is infinite
        PolynomialPhaseSignal siggen = new PolynomialPhaseSignal(N, new Gaussian(0,noisevar), b);
        siggen.generateReceivedSignal();
        Complex[] y = new Complex[N];
        for(int n = 0; n < N; n++) y[n] = new Complex(siggen.real()[n], siggen.imag()[n]);
        try { //write STFT to a binary format suitable for gnuplot
            FileOutputStream fos = new FileOutputStream("testdata/stftimageplotdata");
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            writeFloatToStream(bos,N);
            for(int t = 0; t < N; t++) writeFloatToStream(bos,t);
            double fstep = 1.0/h/4.0; //4 times oversampling in frequency dimension
            for(double f = 0.5; f > -0.5; f -= fstep) {
                writeFloatToStream(bos,(float)f);
                for(int t = 0; t < N; t++) writeFloatToStream(bos,(float)QML.stft(t, f, h, y).abs());
            }
            bos.close();
        } catch (IOException ex) {
            fail("Failed to open file");
        }
        
        try { //write maximiser of STFT to file
            BufferedWriter file = new BufferedWriter(new FileWriter(new File("testdata/stftmax")));
             for(int t = 0; t < N; t++) file.write(t + "\t" + QML.max_stft(t, h, y) + "\n");
             file.close();
        } catch (IOException ex) {
            fail("Failed to open file");
        }
        
      }

    //Write in little endian format (stupid java)
    public static void writeFloatToStream(BufferedOutputStream bos, float f) throws IOException{
        ByteBuffer bbuf = ByteBuffer.allocate(4);
        bbuf.order(ByteOrder.LITTLE_ENDIAN);
        bbuf.putFloat(f);
        bos.write(bbuf.array());
    }
    
}
