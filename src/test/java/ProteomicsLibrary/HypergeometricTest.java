package ProteomicsLibrary;

import ProteomicsLibrary.Types.Hypergeometric;
import org.junit.Test;
import static org.junit.Assert.*;

public class HypergeometricTest {

    @Test
    public void calPMF() throws Exception {
        Hypergeometric hypergeometric = new Hypergeometric(101);
        assertEquals(0.416667, hypergeometric.calPMF(2, 1, 3, 4), 0.00001);
        assertEquals(0.184031, hypergeometric.calPMF(12, 11, 42, 35), 0.00001);
    }

    @Test(expected = Exception.class)
    public void calPMF2() throws Exception {
        Hypergeometric hypergeometric = new Hypergeometric(10);
        assertEquals(0.416667, hypergeometric.calPMF(2, 1, 3, 4), 0.00001);
    }

    @Test(expected = Exception.class)
    public void calPMF3() throws Exception {
        Hypergeometric hypergeometric = new Hypergeometric(9);
        assertEquals(0.416667, hypergeometric.calPMF(2, 1, 3, 4), 0.00001);
    }
}
