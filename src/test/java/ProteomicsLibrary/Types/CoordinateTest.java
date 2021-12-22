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

import org.junit.Test;

import static org.junit.Assert.*;

public class CoordinateTest {

    @Test
    public void toStringTest() {
        Coordinate co = new Coordinate(1, 2);
        assertEquals("(1-2)", co.toString());
        co = new Coordinate(-3, 5);
        assertEquals("(-3-5)", co.toString());
    }

    @Test
    public void compareTo() {
        Coordinate co1 = new Coordinate(1, 4);
        Coordinate co2 = new Coordinate(1, 7);
        Coordinate co3 = new Coordinate(2, 3);
        Coordinate co4 = new Coordinate(2, 5);
        assertTrue(co1.compareTo(co2) < 0);
        assertTrue(co1.compareTo(co3) < 0);
        assertTrue(co1.compareTo(co4) < 0);
    }

    @Test
    public void equals() {
        Coordinate co1 = new Coordinate(1, 4);
        Coordinate co2 = new Coordinate(1, 4);
        Coordinate co3 = new Coordinate(1, 3);
        Coordinate co4 = new Coordinate(2, 4);
        Coordinate co5 = new Coordinate(3, 5);
        assertEquals(co1, co2);
        assertNotEquals(co1, co3);
        assertNotEquals(co1, co4);
        assertNotEquals(co1, co5);
    }
}
