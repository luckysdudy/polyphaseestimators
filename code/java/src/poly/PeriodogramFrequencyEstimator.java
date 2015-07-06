package poly;

import pubsim.Complex;
import pubsim.FFT;

/**
 * Periodogram estimator of a single frequency.  Uses JTransforms FFT
 * @author Robby McKilliam
 */
public class PeriodogramFrequencyEstimator {

    protected final int oversampled,  N;
    protected final FFT fft;
    protected final Complex[] sig;

    /**Max number of iterations for the Newton step */
    protected static final int MAX_ITER = 15;

    /**Step variable for the Newton step */
    protected static final double EPSILON = 1e-10;

    /** oversampling defaults to 4 as per Rife and Borstyn */
    public PeriodogramFrequencyEstimator(int N) {
        this(N,4);
    }

    /** Constructor that sets the number of samples to be taken of
     * the periodogram.
     */
    public PeriodogramFrequencyEstimator(int N, int oversampled) {
        this.oversampled = oversampled;
        this.N = N;
        sig = new Complex[oversampled * N];
        fft = new FFT(oversampled * N);
    }
    
    public double estimateFreq(double[] real, double[] imag) {
        for (int i = 0; i < N; i++) {
            sig[i] = new Complex(real[i], imag[i]);
        }
        for (int i = N; i < sig.length; i++) {
            sig[i] = new Complex(0.0, 0.0);
        }

        Complex[] ft = fft.forward(sig);

        //FFT output is backwards, flip it
        double maxp = 0;
        double fhat = 0.0;
        double f = 0.0;
        double fstep = 1.0 / ft.length;
        for (int i = 0; i < ft.length; i++) {
            double p = ft[i].abs2();
            if (p > maxp) {
                maxp = p;
                fhat = f;
            }
            f+=fstep;
        }

        //Newton Raphson
        int numIter = 0;
        f = fhat;
        double lastf = f - 2 * EPSILON, lastp = 0;
        while (Math.abs(f - lastf) > EPSILON && numIter <= MAX_ITER) {
            double p = 0, pd = 0, pdd = 0;
            double sumur = 0, sumui = 0, sumvr = 0, sumvi = 0,
                    sumwr = 0, sumwi = 0;
            for (int i = 0; i < N; i++) {
                double cosf = Math.cos(-2 * Math.PI * f * i);
                double sinf = Math.sin(-2 * Math.PI * f * i);
                double ur = real[i] * cosf - imag[i] * sinf;
                double ui = imag[i] * cosf + real[i] * sinf;
                double vr = 2 * Math.PI * i * ui;
                double vi = -2 * Math.PI * i * ur;
                double wr = 2 * Math.PI * i * vi;
                double wi = -2 * Math.PI * i * vr;
                sumur += ur;
                sumui += ui;
                sumvr += vr;
                sumvi += vi;
                sumwr += wr;
                sumwi += wi;
            }
            p = sumur * sumur + sumui * sumui;
            if (p < lastp) //I am not sure this is necessary, Vaughan did it for period estimation.
            {
                f = (f + lastf) / 2;
            } else {
                lastf = f;
                lastp = p;
                if (p > maxp) {
                    maxp = p;
                    fhat = f;
                }
                pd = 2 * (sumvr * sumur + sumvi * sumui);
                pdd = 2 * (sumvr * sumvr + sumwr * sumur + sumvi * sumvi + sumwi * sumui);
                f += pd / Math.abs(pdd);
            }
            numIter++;
        }

        return fhat - Math.round(fhat);
    }
}
