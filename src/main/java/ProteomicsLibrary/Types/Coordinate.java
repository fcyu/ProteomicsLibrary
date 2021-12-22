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

public class Coordinate implements Comparable<Coordinate> {

    public final int x;
    public final int y;
    private final int hashCode;

    public Coordinate(int x, int y) {
        this.x = x;
        this.y = y;
        String toString = "(" + x + "-" + y + ")";
        hashCode = toString.hashCode();
    }

    public String toString() {
        return "(" + x + "-" + y + ")";
    }

    public int compareTo(Coordinate other) {
        if (x > other.x) {
            return 1;
        } else if (x < other.x) {
            return -1;
        } else {
            if (y > other.y) {
                return 1;
            } else if (y < other.y) {
                return -1;
            } else {
                return 0;
            }
        }
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof Coordinate) {
            Coordinate temp = (Coordinate) other;
            return temp.hashCode() == this.hashCode();
        } else {
            return false;
        }
    }
}
