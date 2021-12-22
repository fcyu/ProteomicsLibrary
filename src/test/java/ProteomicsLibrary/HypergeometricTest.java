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
