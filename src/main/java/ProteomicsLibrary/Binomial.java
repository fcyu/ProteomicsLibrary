package ProteomicsLibrary;

import java.util.Locale;

public class Binomial {

    private final double[] log10FactorialArray;

    public Binomial(int maxValue) {
        log10FactorialArray = new double[maxValue];
        double temp = 0;
        log10FactorialArray[0] = 0;
        for (int i = 1; i < maxValue; ++i) {
            temp += Math.log10(i);
            log10FactorialArray[i] = temp;
        }
    }

    public double calPValue(int N, int k, double p) throws Exception {
        if (k > N) {
            throw new Exception(String.format(Locale.US, "k = %d is larger than N = %d", k, N));
        }
        if (N < 0 || k < 0) {
            throw new Exception("Either N or k is smaller than 0");
        }
        if (N == 0 && k > 0) {
            throw new Exception("N == 0 but k > 0.");
        }
        if (N > log10FactorialArray.length) {
            throw new IndexOutOfBoundsException(String.format(Locale.US, "N = %d is larger than the maxValue (%d) allowed. Please create a Binomial object with a larger maxValue.", N, log10FactorialArray.length));
        }

        if (N == 0 && k == 0) {
            return 1;
        } else {
            double pValue = 0;
            for (int i = k; i <= N; ++i) {
                pValue += Math.pow(10, calLog10PMF(N, i, p));
            }
            return pValue;
        }
    }

    private double calLog10PMF(int N, int i, double p) {
        return log10FactorialArray[N] - log10FactorialArray[i] - log10FactorialArray[N - i] + i * Math.log10(p) + (N - i) * Math.log10(1 - p);
    }
}
