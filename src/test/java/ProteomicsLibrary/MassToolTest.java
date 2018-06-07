package ProteomicsLibrary;

import org.junit.Before;
import org.junit.Test;
import ProteomicsLibrary.Types.*;

import java.util.*;
import java.util.regex.Pattern;

import static org.junit.Assert.*;


public class MassToolTest {

    private static Map<Character, Double> fixModMap = new HashMap<>();

    @Before
    public void setUp() {
        fixModMap.put('G', 0d);
        fixModMap.put('A', 0d);
        fixModMap.put('S', 0d);
        fixModMap.put('P', 0d);
        fixModMap.put('V', 0d);
        fixModMap.put('T', 0d);
        fixModMap.put('C', 57.02146);
        fixModMap.put('I', 0d);
        fixModMap.put('L', 0d);
        fixModMap.put('N', 0d);
        fixModMap.put('D', 0d);
        fixModMap.put('Q', 0d);
        fixModMap.put('K', 0d);
        fixModMap.put('E', 0d);
        fixModMap.put('M', 0d);
        fixModMap.put('H', 0d);
        fixModMap.put('F', 0d);
        fixModMap.put('R', 0d);
        fixModMap.put('Y', 0d);
        fixModMap.put('W', 0d);
        fixModMap.put('U', 0d);
        fixModMap.put('O', 0d);
        fixModMap.put('n', 60d);
        fixModMap.put('c', 10d);
    }

    @Test
    public void calResidueMass() {
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(2503.1357421875, massTool.calResidueMass("nGASPVTCILNDQKEMHFRYWc"), 0.001);
    }

    @Test
    public void calResidueMass2() {
        MassTool massTool = new MassTool(1, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(2376.11438353312, massTool.calResidueMass2("nGASPVTCILNDQKEMHFRYWc"), 0.001);
        assertEquals(2377.11438353312, massTool.calResidueMass2("nGASPVT(1.0)CILNDQKEMHFRYWc"), 0.001);
        assertEquals(2377.11438353312, massTool.calResidueMass2("nGASPVTC(1.0)ILNDQKEMHFRYWc"), 0.001);
        assertEquals(2433.13584353312, massTool.calResidueMass2("nGASPVTC(57.02146)ILNDQKEMHFRYWc"), 0.001);
    }

    @Test
    public void mzToBin() {
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(11, massTool.mzToBin(11), 1e-6);
        assertEquals(0, massTool.mzToBin(0), 1e-6);
        assertEquals(-5, massTool.mzToBin(-5), 1e-6);
    }

    @Test
    public void buildChainSet() {
        // 1 missed-cleavage, N-term linkable
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");
        Set<String> result = massTool.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        Set<String> groundTruth = new HashSet<>();
        groundTruth.add("nMRc");
        groundTruth.add("nRc");
        groundTruth.add("nMRGFASSASRc");
        groundTruth.add("nRGFASSASRc");
        groundTruth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        groundTruth.add("nLSKTMDYMRc");
        groundTruth.add("nVFKTYc");

        String[] groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        String[] resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);

        // 2 missed-cleavage, N-term linkable
        massTool = new MassTool(2, fixModMap, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        groundTruth = new HashSet<>();
        groundTruth.add("nMRc");
        groundTruth.add("nRc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nMRGFASSASRc");
        groundTruth.add("nRGFASSASRc");
        groundTruth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        groundTruth.add("nLSKTMDYMRc");
        groundTruth.add("nVFKTYc");
        groundTruth.add("nMRGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nRGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nGFASSASRIATAAAASKPSLNASTSVNPKLSKc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKLSKTMDYMRc");
        groundTruth.add("nLSKTMDYMRIFSVFVVTLWIIRc");
        groundTruth.add("nVDARVFKTYc");

        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);
    }

    @Test
    public void buildPeptideSet() {
        MassTool massTool = new MassTool(2, "KR", "P", true, 1.0005 * 0.5, 0.6, "N14");

        String proteinSequence = "MSDDFKDEDRDDKPSSDKKDF";
        Set<String> result = massTool.buildPeptideSet(proteinSequence);

        Set<String> groundTruth = new HashSet<>();
        groundTruth.add("nMSDDFKc");
        groundTruth.add("nSDDFKc");
        groundTruth.add("nDEDRc");
        groundTruth.add("nDDKPSSDKc");
        groundTruth.add("nKc");
        groundTruth.add("nDFc");
        groundTruth.add("nMSDDFKDEDRc");
        groundTruth.add("nSDDFKDEDRc");
        groundTruth.add("nDEDRDDKPSSDKc");
        groundTruth.add("nDDKPSSDKKc");
        groundTruth.add("nKDFc");
        groundTruth.add("nMSDDFKDEDRDDKPSSDKc");
        groundTruth.add("nSDDFKDEDRDDKPSSDKc");
        groundTruth.add("nDEDRDDKPSSDKKc");
        groundTruth.add("nDDKPSSDKKDFc");

        String[] resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        String[] groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        assertArrayEquals(groundTruthArray, resultArray);

    }

    @Test
    public void seqToAAList() {
        String seq = "nGHUKc";
        AA[] result = MassTool.seqToAAList(seq);
        AA[] groundTruth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 0), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(groundTruth, result);

        seq = "nGH[3.02]UKc";
        result = MassTool.seqToAAList(seq);
        groundTruth = new AA[]{new AA('n', 0), new AA('G', 0), new AA('H', 3.02), new AA('U', 0), new AA('K', 0), new AA('c', 0)};
        assertArrayEquals(groundTruth, result);
    }

    @Test
    public void isAA() {
        assertTrue(MassTool.isAA('A'));
        assertFalse(MassTool.isAA('B'));
        assertFalse(MassTool.isAA('-'));
        assertFalse(MassTool.isAA('Z'));
        assertTrue(MassTool.isAA('K'));
    }

    @Test
    public void containsNonAAAndNC() {
        assertFalse(MassTool.containsNonAAAndNC("ACDS"));
        assertTrue(MassTool.containsNonAAAndNC(".SDS"));
        assertTrue(MassTool.containsNonAAAndNC("XSD"));
        assertTrue(MassTool.containsNonAAAndNC("*"));
        assertFalse(MassTool.containsNonAAAndNC("nSDEc"));
    }

    @Test
    public void binToMz() {
        MassTool massTool = new MassTool(1, "KR", "P", true, 0.5, 0.6, "N14");
        assertEquals(1.4, massTool.binToMz(2), 1e-6);
        assertEquals(-0.6, massTool.binToMz(0), 1e-6);
        assertEquals(0.4, massTool.binToMz(1), 1e-6);
    }

    @Test
    public void getLabelling() {
        MassTool massTool = new MassTool(1, "KR", "P", true, 0.5, 0.6, "N14");
        assertEquals("N14", massTool.getLabelling());
        massTool = new MassTool(1, "KR", "P", true, 0.5, 0.6, "N15");
        assertEquals("N15", massTool.getLabelling());
    }

    @Test
    public void getDigestSitePattern() {
        MassTool massTool = new MassTool(1, "KR", "P", true, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[KR](?![P])").toString(), massTool.getDigestSitePattern().toString());
        massTool = new MassTool(1, "KR", "P", false, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("(?<![P])[KR]").toString(), massTool.getDigestSitePattern().toString());
        massTool = new MassTool(1, "K", "-", true, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[K]").toString(), massTool.getDigestSitePattern().toString());
    }
    @Test
    public void aaListToSeq() {
        AA[] aaArray = new AA[]{new AA('A', 0), new AA('B', 0.1), new AA('C', 90), new AA('D', -9.8787), new AA('E', 90.0978)};
        assertEquals("ABC(90.000)D(-9.879)E(90.098)", MassTool.aaListToSeq(aaArray));
    }

    @Test
    public void unifyPeptide() {
        assertEquals("nSDFSDSc", MassTool.unifyPeptide("A.SDFSDS.S"));
        assertEquals("nS(11.320)DFSDSc", MassTool.unifyPeptide("-.S[11.320233]DFSDS.S"));
        assertEquals("nSDFSDSc", MassTool.unifyPeptide("SDFSDS.S"));
        assertEquals("nSDF(-23.231)SDSc", MassTool.unifyPeptide("nSDF(-23.231)SDS.S"));
        assertEquals("nSDF(-23.231)SDSc", MassTool.unifyPeptide("nSDF[-23.231]SDS.-"));
    }

    @Test
    public void L2I() {
        assertEquals("nSDEIISDc", MassTool.L2I("nSDEILSDc"));
        assertEquals("nSDEIISDc", MassTool.L2I("nSDELLSDc"));
    }
}