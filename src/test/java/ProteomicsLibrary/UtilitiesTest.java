package ProteomicsLibrary;

import org.junit.Test;

import java.util.regex.Matcher;

import static org.junit.Assert.*;

public class UtilitiesTest {

    @Test
    public void getScanNum() {
        String s1 = "R2-15.4.4. File:\"R2-15.raw\", NativeID:\"controllerType=0 controllerNumber=1 scan=4\"";
        String s2 = "TITLE=D4.15.15. File:\"\", NativeID:\"sample=1 period=1 cycle=14 experiment=2\"";
        String s3 = "OR20080527_S_mix8FTFTprof-prof_06.00008.00008.3";
        String s4 = "Scan 1643 (rt=24.2110) [Prospector Created]";
        String s5 = "query:11984;rank:1;spectrum:NF1.3099.3099.2_File:\"NF1.raw\",_NativeID:\"controllerType=0_controllerNumber=1_scan=3099\";rt:576.763422;mz:892.941627;charge:2";
        String s6 = "Cmpd 10231, +MS2(324.1625), 23.0eV, 43.8 min #11182";
        assertEquals(4, Utilities.getScanNum(s1));
        assertEquals(15, Utilities.getScanNum(s2));
        assertEquals(8, Utilities.getScanNum(s3));
        assertEquals(1643, Utilities.getScanNum(s4));
        assertEquals(3099, Utilities.getScanNum(s5));
        assertEquals(11182, Utilities.getScanNum(s6));
    }

    @Test(expected = NullPointerException.class)
    public void getScanNum2() {
        assertEquals(3, Utilities.getScanNum("sdsdsd"));
    }

    @Test
    public void testCsvSplitPattern() {
        String s = "1,2,\"3,4,5\",2,\"3,4,\",";
        String[] result = Utilities.csvSplitPattern.split(s);
        String[] groundTruth = new String[]{"1", "2", "\"3,4,5\"", "2", "\"3,4,\""};
        assertArrayEquals(groundTruth, result);
    }

    @Test
    public void testTsvSplitPattern() {
        String s = "1\t2\t\"3\t4,5\"\t2\t\"3,4\t\"\t";
        String[] result = Utilities.tsvSplitPattern.split(s);
        String[] groundTruth = new String[]{"1", "2", "\"3\t4,5\"", "2", "\"3,4\t\""};
        assertArrayEquals(groundTruth, result);
    }
}