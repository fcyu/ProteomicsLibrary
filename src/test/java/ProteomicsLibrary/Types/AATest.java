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
        AA aa12 = new AA('S',79.9901);
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
}
