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
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(2503.1357421875, massTool.calResidueMass("nGASPVTCILNDQKEMHFRYWc"), 0.001);
    }

    @Test
    public void calResidueMass2() {
        MassTool massTool = new MassTool(1, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(2376.11438353312, massTool.calResidueMass2("nGASPVTCILNDQKEMHFRYWc"), 0.001);
        assertEquals(2377.11438353312, massTool.calResidueMass2("nGASPVT(1.0)CILNDQKEMHFRYWc"), 0.001);
        assertEquals(2377.11438353312, massTool.calResidueMass2("nGASPVTC(1.0)ILNDQKEMHFRYWc"), 0.001);
        assertEquals(2433.13584353312, massTool.calResidueMass2("nGASPVTC(57.02146)ILNDQKEMHFRYWc"), 0.001);
    }

    @Test
    public void mzToBin() {
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        assertEquals(11, massTool.mzToBin(11), 1e-6);
        assertEquals(0, massTool.mzToBin(0), 1e-6);
        assertEquals(-5, massTool.mzToBin(-5), 1e-6);
    }

    @Test
    public void buildChainSet() {
        MassTool massTool = new MassTool(1, fixModMap, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
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

        massTool = new MassTool(1, fixModMap, "D", "-", false, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        groundTruth = new HashSet<>();
        groundTruth.add("nMRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMc");
        groundTruth.add("nRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMc");
        groundTruth.add("nMRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVc");
        groundTruth.add("nRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVc");
        groundTruth.add("nDYMRIFSVFVVTLWIIRVDARVFKTYc");
        groundTruth.add("nDARVFKTYc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);

        massTool = new MassTool(1, fixModMap, "D", "-", false, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildChainSet("MRGFACSSASRIATAAAASKPSLNCASTSVNPKLSKTMDYMRICFSVFVVTLWIIRVDARVFKTY", (short) 2);
        groundTruth = new HashSet<>();
        groundTruth.add("nMRGFACSSASRIATAAAASKPSLNCASTSVNPKLSKTMc");
        groundTruth.add("nRGFACSSASRIATAAAASKPSLNCASTSVNPKLSKTMc");
        groundTruth.add("nMRGFACSSASRIATAAAASKPSLNCASTSVNPKLSKTMDYMRICFSVFVVTLWIIRVc");
        groundTruth.add("nRGFACSSASRIATAAAASKPSLNCASTSVNPKLSKTMDYMRICFSVFVVTLWIIRVc");
        groundTruth.add("nDYMRICFSVFVVTLWIIRVc");
        groundTruth.add("nDYMRICFSVFVVTLWIIRVDARVFKTYc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);

        massTool = new MassTool(2, fixModMap, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
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

        massTool = new MassTool(1, fixModMap, "KR", "P", true, "FYWL", "-", true, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        groundTruth = new HashSet<>();
        groundTruth.add("nRc");
        groundTruth.add("nMRc");
        groundTruth.add("nRGFc");
        groundTruth.add("nKTYc");
        groundTruth.add("nMRGFc");
        groundTruth.add("nSKTMDYc");
        groundTruth.add("nNASTSVNPKLc");
        groundTruth.add("nIATAAAASKPSLc");
        groundTruth.add("nASSASRIATAAAASKPSLc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);

        massTool = new MassTool(1, fixModMap, "KR", "P", true, "D", "-", false, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildChainSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY", (short) 1);
        groundTruth = new HashSet<>();
        groundTruth.add("nRc");
        groundTruth.add("nMRc");
        groundTruth.add("nLSKTMc");
        groundTruth.add("nVFKTYc");
        groundTruth.add("nRGFASSASRc");
        groundTruth.add("nMRGFASSASRc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        groundTruth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);
    }

    @Test
    public void buildPeptideSet() {
        MassTool massTool = new MassTool(2, "KR", "P", true, null, null, null, 1.0005 * 0.5, 0.6, "N14");
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

        massTool = new MassTool(2, "D", "-", false, null, null, null, 1.0005 * 0.5, 0.6, "N14");
        proteinSequence = "MSDDFKDEDRPDKPSSDKKDF";
        result = massTool.buildPeptideSet(proteinSequence);
        groundTruth = new HashSet<>();
        groundTruth.add("nMSc");
        groundTruth.add("nMSDc");
        groundTruth.add("nMSDDFKc");
        groundTruth.add("nDc");
        groundTruth.add("nDDFKc");
        groundTruth.add("nDDFKDEc");
        groundTruth.add("nDFKc");
        groundTruth.add("nDFKDEc");
        groundTruth.add("nDFKDEDRPc");
        groundTruth.add("nDEc");
        groundTruth.add("nDEDRPc");
        groundTruth.add("nDEDRPDKPSSc");
        groundTruth.add("nDRPc");
        groundTruth.add("nDRPDKPSSc");
        groundTruth.add("nDRPDKPSSDKKc");
        groundTruth.add("nDKPSSc");
        groundTruth.add("nDKPSSDKKc");
        groundTruth.add("nDKPSSDKKDFc");
        groundTruth.add("nDKKc");
        groundTruth.add("nDKKDFc");
        groundTruth.add("nDFc");
        groundTruth.add("nSc");
        groundTruth.add("nSDc");
        groundTruth.add("nSDDFKc");
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        assertArrayEquals(groundTruthArray, resultArray);

        massTool = new MassTool(1, fixModMap, "KR", "P", true, "FYWL", "-", true, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildPeptideSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY");
        groundTruth = new HashSet<>();
        groundTruth.add("nLc");
        groundTruth.add("nKc");
        groundTruth.add("nWc");
        groundTruth.add("nRc");
        groundTruth.add("nGFc");
        groundTruth.add("nSKc");
        groundTruth.add("nVFc");
        groundTruth.add("nIFc");
        groundTruth.add("nTYc");
        groundTruth.add("nMRc");
        groundTruth.add("nLSKc");
        groundTruth.add("nRGFc");
        groundTruth.add("nSVFc");
        groundTruth.add("nVFKc");
        groundTruth.add("nIIRc");
        groundTruth.add("nKTYc");
        groundTruth.add("nVVTLc");
        groundTruth.add("nVDARc");
        groundTruth.add("nMRGFc");
        groundTruth.add("nTMDYc");
        groundTruth.add("nMRIFc");
        groundTruth.add("nASSASRc");
        groundTruth.add("nWIIRc");
        groundTruth.add("nIFSVFc");
        groundTruth.add("nVVTLWc");
        groundTruth.add("nVDARVFc");
        groundTruth.add("nSKTMDYc");
        groundTruth.add("nSVFVVTLc");
        groundTruth.add("nGFASSASRc");
        groundTruth.add("nTMDYMRc");
        groundTruth.add("nIIRVDARc");
        groundTruth.add("nNASTSVNPKc");
        groundTruth.add("nNASTSVNPKLc");
        groundTruth.add("nIATAAAASKPSLc");
        groundTruth.add("nASSASRIATAAAASKPSLc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
        assertArrayEquals(groundTruthArray, resultArray);

        massTool = new MassTool(1, fixModMap, "KR", "P", true, "D", "-", false, 1.0005 * 0.5, 0.6, "N14");
        result = massTool.buildPeptideSet("MRGFASSASRIATAAAASKPSLNASTSVNPKLSKTMDYMRIFSVFVVTLWIIRVDARVFKTY");
        groundTruth = new HashSet<>();
        groundTruth.add("nVc");
        groundTruth.add("nRc");
        groundTruth.add("nTMc");
        groundTruth.add("nTYc");
        groundTruth.add("nMRc");
        groundTruth.add("nLSKc");
        groundTruth.add("nDARc");
        groundTruth.add("nVFKc");
        groundTruth.add("nVDARc");
        groundTruth.add("nLSKTMc");
        groundTruth.add("nDYMRc");
        groundTruth.add("nVFKTYc");
        groundTruth.add("nDARVFKc");
        groundTruth.add("nGFASSASRc");
        groundTruth.add("nTMDYMRc");
        groundTruth.add("nRGFASSASRc");
        groundTruth.add("nMRGFASSASRc");
        groundTruth.add("nIFSVFVVTLWIIRc");
        groundTruth.add("nIFSVFVVTLWIIRVc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKc");
        groundTruth.add("nDYMRIFSVFVVTLWIIRc");
        groundTruth.add("nIATAAAASKPSLNASTSVNPKLSKc");
        groundTruth.add("nGFASSASRIATAAAASKPSLNASTSVNPKc");
        groundTruthArray = groundTruth.toArray(new String[0]);
        Arrays.sort(groundTruthArray);
        resultArray = result.toArray(new String[0]);
        Arrays.sort(resultArray);
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
        MassTool massTool = new MassTool(1, "KR", "P", true, null, null, null, 0.5, 0.6, "N14");
        assertEquals(1.4, massTool.binToMz(2), 1e-6);
        assertEquals(-0.6, massTool.binToMz(0), 1e-6);
        assertEquals(0.4, massTool.binToMz(1), 1e-6);
    }

    @Test
    public void getLabelling() {
        MassTool massTool = new MassTool(1, "KR", "P", true, null, null, null, 0.5, 0.6, "N14");
        assertEquals("N14", massTool.getLabelling());
        massTool = new MassTool(1, "KR", "P", true, null, null, null, 0.5, 0.6, "N15");
        assertEquals("N15", massTool.getLabelling());
    }

    @Test
    public void getDigestSitePattern() {
        MassTool massTool = new MassTool(1, "KR", "P", true, null, null, null, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[KR](?![P])").toString(), massTool.getDigestSitePattern1().toString());
        assertNull(massTool.getDigestSitePattern2());
        massTool = new MassTool(1, "KR", "P", false, null, null, null, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("(?<![P])[KR]").toString(), massTool.getDigestSitePattern1().toString());
        assertNull(massTool.getDigestSitePattern2());
        massTool = new MassTool(1, "K", "-", true, null, null, null, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[K]").toString(), massTool.getDigestSitePattern1().toString());
        assertNull(massTool.getDigestSitePattern2());
        massTool = new MassTool(1, "K", "-", false, null, null, null, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[K]").toString(), massTool.getDigestSitePattern1().toString());
        assertNull(massTool.getDigestSitePattern2());
        massTool = new MassTool(1, "KR", "P", true, "DE", "P", true, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[KR](?![P])").toString(), massTool.getDigestSitePattern1().toString());
        assertEquals(Pattern.compile("[DE](?![P])").toString(), massTool.getDigestSitePattern2().toString());
        massTool = new MassTool(1, "KR", "-", true, "D", "-", false, 0.5, 0.6, "N14");
        assertEquals(Pattern.compile("[KR]").toString(), massTool.getDigestSitePattern1().toString());
        assertEquals(Pattern.compile("[D]").toString(), massTool.getDigestSitePattern2().toString());
    }
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