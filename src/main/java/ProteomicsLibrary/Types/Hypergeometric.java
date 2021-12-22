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
