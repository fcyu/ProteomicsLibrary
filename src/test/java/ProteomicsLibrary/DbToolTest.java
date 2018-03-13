package ProteomicsLibrary;

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
}