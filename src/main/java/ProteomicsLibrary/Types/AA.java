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

public class AA implements Comparable<AA> {

    public final char aa;
    public final double ptmDeltaMass;
    private final String toString;
    private final int hashCode;

    public AA(char aa, double ptmDeltaMass) {
        this.aa = aa;
        this.ptmDeltaMass = ptmDeltaMass;

        if (hasMod()) { // If the PTM is smaller than or equal to 0.1, there will be no string for the PTM.
            toString = String.format(Locale.US, "%c(%.3f)", aa, ptmDeltaMass);
        } else {
            toString =  String.valueOf(aa);
        }
        hashCode = toString.hashCode();
    }

    public boolean hasMod() {
        return Math.abs(ptmDeltaMass) > 0.1;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof AA) {
            return ((AA) other).hashCode == hashCode;
        } else {
            return false;
        }
    }

    public AA clone() {
        return new AA(aa, ptmDeltaMass);
    }

    public String toString() {
        return toString;
    }

    public int compareTo(AA other) {
        if (aa < other.aa) {
            return -1;
        } else if (aa > other.aa) {
            return 1;
        } else {
            if (ptmDeltaMass < other.ptmDeltaMass) {
                return -1;
            } else if (ptmDeltaMass > other.ptmDeltaMass) {
                return 1;
            } else {
                return 0;
            }
        }
    }
}
