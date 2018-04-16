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
