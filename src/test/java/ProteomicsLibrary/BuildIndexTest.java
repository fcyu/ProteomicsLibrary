package ProteomicsLibrary;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.junit.Before;
import org.junit.Test;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;

public class BuildIndexTest {

    private static Map<String, String> proteinSequenceMap;

    @Before
    public void setUp() {
        proteinSequenceMap = new HashMap<>();
        proteinSequenceMap.put("pro1", "SSDSDSKSDSRSDSDSKPSDS");
        proteinSequenceMap.put("pro2", "SDSKKSDSRDSSK");
    }

    @Test
    public void getTargetPeptideProteinMap() {
        BuildIndex buildIndex = new BuildIndex(proteinSequenceMap, "KR", "P", true, null, null, null, 1);
        Multimap<String, String> results = buildIndex.getTargetPeptideProteinMap();
        Multimap<String, String> groundTruth = HashMultimap.create();
        groundTruth.put("nSSDSDSKc", "pro1");
        groundTruth.put("nSDSRc", "pro1");
        groundTruth.put("nSDSDSKPSDSc", "pro1");
        groundTruth.put("nSSDSDSKSDSRc", "pro1");
        groundTruth.put("nSDSRSDSDSKPSDSc", "pro1");
        groundTruth.put("nSDSKc", "pro2");
        groundTruth.put("nKc", "pro2");
        groundTruth.put("nSDSRc", "pro2");
        groundTruth.put("nDSSKc", "pro2");
        groundTruth.put("nSDSKKc", "pro2");
        groundTruth.put("nKSDSRc", "pro2");
        groundTruth.put("nSDSRDSSKc", "pro2");
        assertEquals(groundTruth.size(), results.size());
        for (String peptide : groundTruth.keySet()) {
            Collection<String> resultSet = results.get(peptide);
            Collection<String> groundTruthSet = groundTruth.get(peptide);
            assertTrue(resultSet.containsAll(groundTruthSet));
            assertTrue(groundTruthSet.containsAll(resultSet));
        }
    }
}