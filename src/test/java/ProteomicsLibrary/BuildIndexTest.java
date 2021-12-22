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
        proteinSequenceMap.put("pro3", "MABC");
        proteinSequenceMap.put("pro4", "MMDEF");
        proteinSequenceMap.put("pro5", "MMMGHI");
        proteinSequenceMap.put("pro6", "MDMJPL");
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
        groundTruth.put("nMABCc", "pro3");
        groundTruth.put("nABCc", "pro3");
        groundTruth.put("nMMDEFc", "pro4");
        groundTruth.put("nMDEFc", "pro4");
        groundTruth.put("nDEFc", "pro4");
        groundTruth.put("nMMMGHIc", "pro5");
        groundTruth.put("nMMGHIc", "pro5");
        groundTruth.put("nMGHIc", "pro5");
        groundTruth.put("nGHIc", "pro5");
        groundTruth.put("nMDMJPLc", "pro6");
        groundTruth.put("nDMJPLc", "pro6");
        assertEquals(groundTruth.size(), results.size());
        for (String peptide : groundTruth.keySet()) {
            Collection<String> resultSet = results.get(peptide);
            Collection<String> groundTruthSet = groundTruth.get(peptide);
            assertTrue(resultSet.containsAll(groundTruthSet));
            assertTrue(groundTruthSet.containsAll(resultSet));
        }
    }
}