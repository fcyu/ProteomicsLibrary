package ProteomicsLibrary;

import org.junit.Test;

import static org.junit.Assert.*;

public class BinomialTest {

    @Test
    public void calProbLargerThanOrEqualTo() throws Exception {
        Binomial binomial = new Binomial(100);
        assertEquals(0.08146, binomial.calProbLargerThanOrEqualTo(5, 2, 0.1), 0.0001);
        assertEquals(0, binomial.calProbLargerThanOrEqualTo(6, 6, 0.01), 0.0001);
        assertEquals(1, binomial.calProbLargerThanOrEqualTo(7, 0, 0.03), 0.0001);
        assertEquals(1, binomial.calProbLargerThanOrEqualTo(0, 0, 0.1), 0.0001);
    }

    @Test(expected = Exception.class)
    public void calProbLargerThanOrEqualTo2() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbLargerThanOrEqualTo(11, 2, 0.01);
        binomial.calProbLargerThanOrEqualTo(3, 5, 0.01);
        binomial.calProbLargerThanOrEqualTo(0, 1, 0.01);
    }
}