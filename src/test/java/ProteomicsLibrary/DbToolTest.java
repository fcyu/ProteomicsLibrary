package ProteomicsLibrary;

import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.junit.Before;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class DbToolTest {
    private static DbTool db_tool_obj;

    @Before
    public void setUp() throws Exception {
        db_tool_obj = new DbTool(Thread.currentThread().getContextClassLoader().getResource("test.fasta").getPath(), "Others");
    }

    @Test
    public void getProteinSequenceMap() {
        Map<String, String> pro_seq_map = db_tool_obj.getProteinSequenceMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1 test protein one", "ASRIATAAAASKPSLNKF");
        for (String k : pro_seq_map.keySet()) {
            assertEquals(ground_truth.get(k), pro_seq_map.get(k));
        }
    }

    @Test
    public void getProteinAnnotateMap() {
        Map<String, String> pro_annotate_map = db_tool_obj.getProteinAnnotateMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1 test protein one", "Pro1 test protein one");
        for (String k : pro_annotate_map.keySet()) {
            assertEquals(ground_truth.get(k), pro_annotate_map.get(k));
        }
    }

    @Test
    public void shuffleSeq() {
        String sequence = "MJUGHSDKSSDSDKPSDSRSDK";
        String result = DbTool.shuffleSeq(sequence, "KR", "P", true);
        String groundTruth = "MUJHGDSKSSSDKDSPSDRDSK";
        assertEquals(groundTruth, result);

        sequence = "JUGHSDKSSDSDKPSDSRSDKASD";
        result = DbTool.shuffleSeq(sequence, "KR", "P", true);
        groundTruth = "UJHGDSKSSSDKDSPSDRDSKSAD";
        assertEquals(groundTruth, result);

        sequence = "MJUGHSDPKSSDSDKPSDSRSDK";
        result = DbTool.shuffleSeq(sequence, "KR", "P", false);
        groundTruth = "MUJHGDSKPSSSDDKSPSDRDSK";
        assertEquals(groundTruth, result);
    }

    @Test
    public void shuffleSeqFY() {
        String sequence = "MJUGHSDKSSDSDKPSDSRSDK";
        System.out.println(String.format(Locale.US, "Test shuffleSeqFY(%s, %s, %s, %s):", sequence, "KR", "P", "true"));
        String decoySequence = DbTool.shuffleSeqFY(sequence, "KR", "P", true);
        System.out.println("Target:\t" + sequence);
        System.out.println("Decoy:\t" + decoySequence);

        sequence = "JUGHSDKSSDSDKPSDSRSDKASD";
        System.out.println(String.format(Locale.US, "Test shuffleSeqFY(%s, %s, %s, %s):", sequence, "KR", "P", "true"));
        decoySequence = DbTool.shuffleSeqFY(sequence, "KR", "P", true);
        System.out.println("Target:\t" + sequence);
        System.out.println("Decoy:\t" + decoySequence);

        sequence = "MJUGHSDPKSSDSDKPSDSRSDK";
        System.out.println(String.format(Locale.US, "Test shuffleSeqFY(%s, %s, %s, %s):", sequence, "KR", "P", "false"));
        decoySequence = DbTool.shuffleSeqFY(sequence, "KR", "P", false);
        System.out.println("Target:\t" + sequence);
        System.out.println("Decoy:\t" + decoySequence);
    }

    @Test
    public void getLeftRightFlank() {
        try {
            String peptide = "nSDSDD(229.22)SDSDSKc";
            Multimap<String, String> peptideProteinMap = TreeMultimap.create();
            peptideProteinMap.put("nSDSDD(229.22)SDSDSKc", "protein1");
            peptideProteinMap.put("nSDSDD(229.22)SDSDSKc", "protein2");
            peptideProteinMap.put("nSDSDD(229.22)SDSDSKc", "protein3");
            Map<String, String> proteinSequenceMap = new HashMap<>();
            proteinSequenceMap.put("protein1", "SDSDSKSDSDDSDSDSKPFDSDS");
            proteinSequenceMap.put("protein2", "SDSDSASDSDDSDSDSKADS");
            proteinSequenceMap.put("protein3", "SDSDSRSDSDDSDSDSKSDSDS");
            Character[] leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, proteinSequenceMap, "KR", "P", true);
            assertEquals(new Character('R'), leftRightFlank[0]);
            assertEquals( new Character('S'), leftRightFlank[1]);

            peptide = "nRSDSDD(229.22)SDSDSc";
            peptideProteinMap.clear();
            peptideProteinMap.put("nRSDSDD(229.22)SDSDSc", "protein1");
            peptideProteinMap.put("nRSDSDD(229.22)SDSDSc", "protein2");
            peptideProteinMap.put("nRSDSDD(229.22)SDSDSc", "protein3");
            proteinSequenceMap.clear();
            proteinSequenceMap.put("protein1", "PRSDSDDSDSDSRDSDS");
            proteinSequenceMap.put("protein2", "MRSDSDDSDSDSKDSDS");
            proteinSequenceMap.put("protein3", "HRSDSDDSDSDSRDSDS");
            leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, proteinSequenceMap, "KR", "P", false);
            assertEquals(new Character('-'), leftRightFlank[0]);
            assertEquals(new Character('K'), leftRightFlank[1]);
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

    @Test
    public void findPeptideLocation() {
        Set<Integer> result = DbTool.findPeptideLocation("MAGSYFCDSKCKLRCSKAGLADRCLKYCGICCEECKCVPSGTYGNKHECPCYRDKKNSKGKSKCP*", "nM(-2.946)A(26.016)GS(57.021)YFCDSKc", "KR", "P");
        Set<Integer> groundTruth = new HashSet<>();
        groundTruth.add(0);
        assertEquals(groundTruth.size(), result.size());
        Integer[] resultArray = result.toArray(new Integer[0]);
        Arrays.sort(resultArray);
        Integer[] groundTruthArray = groundTruth.toArray(new Integer[0]);
        Arrays.sort(groundTruthArray);
        for (int i = 0; i < resultArray.length; ++i) {
            assertEquals(groundTruthArray[i], resultArray[i]);
        }
    }

    @Test
    public void reduceProteinIdSet() {
        Set<String> inputSet = new HashSet<>();
        inputSet.add("A.1");
        inputSet.add("A.2");
        inputSet.add("B.2");
        Set<String> results = DbTool.reduceProteinIdSet(inputSet, "TAIR");
        Set<String> groundTruth = new HashSet<>();
        groundTruth.add("A.1");
        groundTruth.add("B.2");
        assertTrue(groundTruth.containsAll(results));
        assertTrue(results.containsAll(groundTruth));
        results = DbTool.reduceProteinIdSet(inputSet,"UniProt");
        assertFalse(groundTruth.containsAll(results));
        assertTrue(results.containsAll(groundTruth));
    }

    @Test
    public void getPtmFreePeptide() {
        String peptide = "n(34.063)IINEPTAAAIAYGLDKK(-34.063)c";
        String groundTruth = "nIINEPTAAAIAYGLDKKc";
        assertEquals(groundTruth, DbTool.getPtmFreePeptide(peptide));
    }

    @Test
    public void getCLPtmFreePeptide() {
        String peptide = "n(34.063)KYGSGGAc-1-n(34.063)IINEPTAAAIAYGLDKK(-34.063)c-16";
        String result = DbTool.getCLPtmFreePeptide(peptide);
        String groundTruth = "nKYGSGGAc-1-nIINEPTAAAIAYGLDKKc-16";
        assertEquals(groundTruth, result);

        peptide = "n[34.063]KYGSGGAc-1-n[34.063]IINEPTAAAIAYGLDKK[-34.063]c-16";
        result = DbTool.getCLPtmFreePeptide(peptide);
        groundTruth = "nKYGSGGAc-1-nIINEPTAAAIAYGLDKKc-16";
        assertEquals(groundTruth, result);
    }

    @Test
    public void getCLSequenceOnly() {
        String peptide = "n(34.063)KYGSGGAc-1-n(34.063)IINEPTAAAIAYGLDKK(-34.063)c-16";
        String result = DbTool.getCLSequenceOnly(peptide);
        String groundTruth = "KYGSGGA-1-IINEPTAAAIAYGLDKK-16";
        assertEquals(groundTruth, result);

        peptide = "n[34.063]KYGSGGAc-1-n[34.063]IINEPTAAAIAYGLDKK[-34.063]c-16";
        result = DbTool.getCLSequenceOnly(peptide);
        groundTruth = "KYGSGGA-1-IINEPTAAAIAYGLDKK-16";
        assertEquals(groundTruth, result);
    }
}