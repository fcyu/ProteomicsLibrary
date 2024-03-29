/*
 * Copyright 2018-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
    }

    @Test(expected = Exception.class)
    public void calProbLargerThanOrEqualTo3() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbLargerThanOrEqualTo(10, 5, 0.01);
    }

    @Test(expected = Exception.class)
    public void calProbLargerThanOrEqualTo4() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbLargerThanOrEqualTo(6, 7, 0.01);
    }

    @Test
    public void calProbSmallerThanOrEqualTo() throws Exception {
        Binomial binomial = new Binomial(100);
        assertEquals(0.99144, binomial.calProbSmallerThanOrEqualTo(5, 2, 0.1), 0.0001);
        assertEquals(1, binomial.calProbSmallerThanOrEqualTo(6, 6, 0.01), 0.0001);
        assertEquals(0.80798, binomial.calProbSmallerThanOrEqualTo(7, 0, 0.03), 0.0001);
        assertEquals(1, binomial.calProbSmallerThanOrEqualTo(0, 0, 0.1), 0.0001);
    }

    @Test(expected = Exception.class)
    public void calProbSmallerThanOrEqualTo2() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbSmallerThanOrEqualTo(11, 2, 0.01);
    }

    @Test(expected = Exception.class)
    public void calProbSmallerThanOrEqualTo3() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbSmallerThanOrEqualTo(10, 5, 0.01);
    }

    @Test(expected = Exception.class)
    public void calProbSmallerThanOrEqualTo4() throws Exception{
        Binomial binomial = new Binomial(10);
        binomial.calProbSmallerThanOrEqualTo(6, 8, 0.01);
    }
}