package ProteomicsLibrary.Types;

import java.util.Locale;

public class Hypergeometric {

    private final double[] log10FactorialArray;

    public Hypergeometric(int maxValue) {
        log10FactorialArray = new double[maxValue];
        double temp = 0;
        log10FactorialArray[0] = 0;
        for (int i = 1; i < maxValue; ++i) {
            temp += Math.log10(i);
            log10FactorialArray[i] = temp;
        }
    }

    public double calPMF(int a, int b, int c, int d) throws Exception {
        // C(a + b, a) * C(c + d, c) / C(a + b + c + d, a + c)
        if (a < 0 || b < 0 || c < 0 || d < 0) {
            throw new Exception("a or b or c or d is smaller than 0.");
        }
        if (a + b + c + d >= log10FactorialArray.length) {
            throw new IndexOutOfBoundsException(String.format(Locale.US, "a + b + c + d = %d is larger than or equal to maxValue (%d). Please create a Hypergeomitric object with a larger maxValue.", a + b + c + d, log10FactorialArray.length));
        }
        return Math.pow(10, log10FactorialArray[a + b] + log10FactorialArray[c + d] + log10FactorialArray[a + c] + log10FactorialArray[b + d] - log10FactorialArray[a + b + c + d] - log10FactorialArray[a] - log10FactorialArray[b] - log10FactorialArray[c] - log10FactorialArray[d]);
    }
}
