package ProteomicsLibrary;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.assertEquals;

public class DbToolTest {
    private static DbTool db_tool_obj;

    @BeforeClass
    public static void setUp() throws Exception {
        db_tool_obj = new DbTool(Thread.currentThread().getContextClassLoader().getResource("test.fasta").getPath(), "Others");
    }

    @Test
    public void returnSeqMap() {
        Map<String, String> pro_seq_map = db_tool_obj.getProteinSequenceMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1 test protein one", "ASRIATAAAASKPSLNKF");
        for (String k : pro_seq_map.keySet()) {
            assertEquals(pro_seq_map.get(k), ground_truth.get(k));
        }
    }

    @Test
    public void findPeptideLocation() {
        Set<Integer> result = DbTool.findPeptideLocation("MAGSYFCDSKCKLRCSKAGLADRCLKYCGICCEECKCVPSGTYGNKHECPCYRDKKNSKGKSKCP*", "nM(-2.946)A(26.016)GS(57.021)YFCDSKc", "KR", "P");
        Set<Integer> groundTruth = new HashSet<>();
        groundTruth.add(0);
        assertEquals(result.size(), groundTruth.size());
        Integer[] resultArray = result.toArray(new Integer[result.size()]);
        Arrays.sort(resultArray);
        Integer[] groundTruthArray = groundTruth.toArray(new Integer[groundTruth.size()]);
        Arrays.sort(groundTruthArray);
        for (int i = 0; i < resultArray.length; ++i) {
            assertEquals(resultArray[i], groundTruthArray[i]);
        }
    }

    @Test
    public void returnAnnotateMap() {
        Map<String, String> pro_annotate_map = db_tool_obj.getProteinAnnotateMap();
        Map<String, String> ground_truth = new HashMap<>();
        ground_truth.put("Pro1 test protein one", "Pro1 test protein one");
        for (String k : pro_annotate_map.keySet()) {
            assertEquals(pro_annotate_map.get(k), ground_truth.get(k));
        }
    }

    @Test
    public void shuffleSeq() {
        String sequence = "MJUGHSDKSSDSDKPSDSRSDK";
        String result = DbTool.shuffleSeq(sequence, "KR", "P", true);
        String groundTruth = "MUJHGDSKSSSDKDSPSDRDSK";
        assertEquals(result, groundTruth);

        sequence = "JUGHSDKSSDSDKPSDSRSDKASD";
        result = DbTool.shuffleSeq(sequence, "KR", "P", true);
        groundTruth = "UJHGDSKSSSDKDSPSDRDSKSAD";
        assertEquals(result, groundTruth);

        sequence = "MJUGHSDPKSSDSDKPSDSRSDK";
        result = DbTool.shuffleSeq(sequence, "KR", "P", false);
        groundTruth = "MUJHGDSKPSSSDDKSPSDRDSK";
        assertEquals(result, groundTruth);
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
            assertEquals(leftRightFlank[0], new Character('R'));
            assertEquals(leftRightFlank[1], new Character('S'));

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
            assertEquals(leftRightFlank[0], new Character('-'));
            assertEquals(leftRightFlank[1], new Character('K'));
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(1);
        }
    }

    @Test
    public void getCLPtmFreePeptide() {
        String peptide = "n(34.063)KYGSGGAc-1-n(34.063)IINEPTAAAIAYGLDKK(-34.063)c-16";
        String result = DbTool.getCLPtmFreePeptide(peptide, "()");
        String groundTruth = "nKYGSGGAc-1-nIINEPTAAAIAYGLDKKc-16";
        assertEquals(result, groundTruth);

        peptide = "n[34.063]KYGSGGAc-1-n[34.063]IINEPTAAAIAYGLDKK[-34.063]c-16";
        result = DbTool.getCLPtmFreePeptide(peptide, "[]");
        groundTruth = "nKYGSGGAc-1-nIINEPTAAAIAYGLDKKc-16";
        assertEquals(result, groundTruth);
    }

    @Test
    public void getCLSequenceOnly() {
        String peptide = "n(34.063)KYGSGGAc-1-n(34.063)IINEPTAAAIAYGLDKK(-34.063)c-16";
        String result = DbTool.getCLSequenceOnly(peptide, "()");
        String groundTruth = "KYGSGGA-1-IINEPTAAAIAYGLDKK-16";
        assertEquals(result, groundTruth);

        peptide = "n[34.063]KYGSGGAc-1-n[34.063]IINEPTAAAIAYGLDKK[-34.063]c-16";
        result = DbTool.getCLSequenceOnly(peptide, "[]");
        groundTruth = "KYGSGGA-1-IINEPTAAAIAYGLDKK-16";
        assertEquals(result, groundTruth);
    }
}