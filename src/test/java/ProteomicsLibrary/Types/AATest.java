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


import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class AATest {

    private static AA aa1;
    private static AA aa2;
    private static AA aa3;
    private static AA aa4;

    @Before
    public void setUp() {
        aa1 = new AA('S',79.99);
        aa2 = new AA('S', 0);
        aa3 = new AA('D', -0.9);
        aa4 = new AA('F', 0.01);
    }

    @Test
    public void hasMod() {
        assertTrue(aa1.hasMod());
        assertFalse(aa2.hasMod());
        assertTrue(aa3.hasMod());
        assertFalse(aa4.hasMod());
    }

    @Test
    public void equals() {
        AA aa11 = new AA('S',79.99);
        AA aa12 = new AA('S',79.991);
        AA aa13 = new AA('F', 79.99);
        assertEquals(aa11, aa1);
        assertNotEquals(aa12, aa1);
        assertNotEquals(aa13, aa1);
    }

    @Test
    public void testToString() {
        assertEquals("A", (new AA('A', 0)).toString());
        assertEquals("B(2.100)", (new AA('B', 2.1)).toString());
        assertEquals("C", (new AA('C', 0.01)).toString());
        assertEquals("D(-2.335)", (new AA('D', -2.33462)).toString());
        assertEquals("E(3.222)", (new AA('E', 3.2219)).toString());
        assertEquals("F(6.333)", (new AA('F', 6.3331)).toString());
    }

    @Test
    public void compareTo() {
        assertTrue((new AA('A', 0)).compareTo(new AA('B', 0)) < 0);
        assertTrue((new AA('C', 0)).compareTo(new AA('B', 0)) > 0);
        assertTrue((new AA('A', 0)).compareTo(new AA('A', 2)) < 0);
        assertTrue((new AA('A', 2)).compareTo(new AA('A', 0)) > 0);
        assertTrue((new AA('A', 3)).compareTo(new AA('B', 0)) < 0);
        assertTrue((new AA('A', 0)).compareTo(new AA('A', 0)) == 0);
    }
}
