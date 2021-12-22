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

import org.junit.Before;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;

public class IsotopeDistributionTest {

    private static IsotopeDistribution isotopeDistribution;

    @Before
    public void setUp() {
        MassTool massTool = new MassTool(2, "KR", "P", true, null, null, null, 0.02, 0, "N14");;
        isotopeDistribution = new IsotopeDistribution(massTool.getElementTable(), 0, "N14");
    }

    @Test
    public void getElementMapFromMonoMass() {
        Map<String, Integer> results = isotopeDistribution.getElementMapFromMonoMass(1000);
        Map<String, Integer> groundTruth = new HashMap<>();
        groundTruth.put("C", 44);
        groundTruth.put("H", 95);
        groundTruth.put("N", 12);
        groundTruth.put("O", 13);
        assertEquals(groundTruth.size(), results.size());
        for (String element : groundTruth.keySet()) {
            assertEquals(groundTruth.get(element), results.get(element));
        }
        results = isotopeDistribution.getElementMapFromMonoMass(1567);
        groundTruth = new HashMap<>();
        groundTruth.put("C", 70);
        groundTruth.put("H", 92);
        groundTruth.put("N", 19);
        groundTruth.put("O", 21);
        groundTruth.put("S", 1);
        assertEquals(groundTruth.size(), results.size());
        for (String element : groundTruth.keySet()) {
            assertEquals(groundTruth.get(element), results.get(element));
        }
    }

    @Test
    public void calculateIsotopePeaks() {
        Map<String, Integer> formulaMap = new HashMap<>();
        formulaMap.put("C", 44);
        formulaMap.put("H", 95);
        formulaMap.put("N", 12);
        formulaMap.put("O", 13);
        List<IsotopeDistribution.Peak> results = isotopeDistribution.calculateIsotopePeaks(formulaMap);
        List<IsotopeDistribution.Peak> groundTruth = new LinkedList<>();
        groundTruth.add(new IsotopeDistribution.Peak(999.71361, 0.5602));
        groundTruth.add(new IsotopeDistribution.Peak(1000.71654, 0.3113));
        groundTruth.add(new IsotopeDistribution.Peak(1001.71923, 0.0998));
        groundTruth.add(new IsotopeDistribution.Peak(1002.72182, 0.0234));
        groundTruth.add(new IsotopeDistribution.Peak(1003.72436, 0.0044));
        groundTruth.add(new IsotopeDistribution.Peak(1004.72686, 0.0007));
        groundTruth.add(new IsotopeDistribution.Peak(1005.72934, 0.0001));
        groundTruth.add(new IsotopeDistribution.Peak(1006.73197, 0));
        groundTruth.add(new IsotopeDistribution.Peak(1007.73424, 0));
        assertEquals(groundTruth.size(), results.size());
        for (int i = 0; i < groundTruth.size(); ++i) {
            assertEquals(groundTruth.get(i).mass, results.get(i).mass, 0.03);
            assertEquals(groundTruth.get(i).realArea, results.get(i).realArea, 1e-4);
        }

        formulaMap.clear();
        formulaMap.put("C", 70);
        formulaMap.put("H", 92);
        formulaMap.put("N", 19);
        formulaMap.put("O", 21);
        formulaMap.put("S", 1);
        results = isotopeDistribution.calculateIsotopePeaks(formulaMap);
        groundTruth.clear();
        groundTruth.add(new IsotopeDistribution.Peak(1566.64304, 0.3812));
        groundTruth.add(new IsotopeDistribution.Peak(1567.64591, 0.3366));
        groundTruth.add(new IsotopeDistribution.Peak(1568.64771, 0.1801));
        groundTruth.add(new IsotopeDistribution.Peak(1569.64929, 0.0714));
        groundTruth.add(new IsotopeDistribution.Peak(1570.65086, 0.0228));
        groundTruth.add(new IsotopeDistribution.Peak(1571.65250, 0.0061));
        groundTruth.add(new IsotopeDistribution.Peak(1572.65422, 0.0014));
        groundTruth.add(new IsotopeDistribution.Peak(1573.65602, 0.0003));
        groundTruth.add(new IsotopeDistribution.Peak(1574.65793, 0.0001));
        groundTruth.add(new IsotopeDistribution.Peak(1575.65982, 0));
        groundTruth.add(new IsotopeDistribution.Peak(1576.66083, 0));
        assertEquals(groundTruth.size(), results.size());
        for (int i = 0; i < groundTruth.size(); ++i) {
            assertEquals(groundTruth.get(i).mass, results.get(i).mass, 0.03);
            assertEquals(groundTruth.get(i).realArea, results.get(i).realArea, 1e-4);
        }
    }

    @Test
    public void getIsotopeCorrectionNum() throws Exception {
        TreeMap<Double, Double> parentPeakList = new TreeMap<>();
        parentPeakList.put(999.96626, 29.68);
        parentPeakList.put(1000.46771, 33.31);
        parentPeakList.put(1000.96878, 21.47);
        parentPeakList.put(1001.46973, 10.10);
        parentPeakList.put(1001.97063, 3.80);
        IsotopeDistribution.Entry entry = isotopeDistribution.getIsotopeCorrectionNum(999.96626, 10, 1, 2, parentPeakList);
        assertEquals(0, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(1000.46771, 10, 1, 2, parentPeakList);
        assertEquals(-1, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(1000.96878, 10, 1, 2, parentPeakList);
        assertEquals(-2, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);

        parentPeakList.clear();
        parentPeakList.put(986.85505, 16.89);
        parentPeakList.put(987.18935, 28.11);
        parentPeakList.put(987.52352, 25.34);
        parentPeakList.put(987.85762, 16.21);
        parentPeakList.put(988.19167, 8.19);
        entry = isotopeDistribution.getIsotopeCorrectionNum(986.85505, 10, 1, 3, parentPeakList);
        assertEquals(0, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(987.18935, 10, 1, 3, parentPeakList);
        assertEquals(-1, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(987.52352, 10, 1, 3, parentPeakList);
        assertEquals(-2, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);

        parentPeakList.clear();
        parentPeakList.put(455.97731, 32d);
        parentPeakList.put(456.22803, 33d);
        parentPeakList.put(456.47854, 20d);
        parentPeakList.put(456.72898, 9d);
        parentPeakList.put(456.97941, 3d);
        entry = isotopeDistribution.getIsotopeCorrectionNum(455.97731, 10, 1, 4, parentPeakList);
        assertEquals(0, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(456.22803, 10, 1, 4, parentPeakList);
        assertEquals(0, entry.isotopeCorrectionNum); // the performance is not good enough to get the correct result.
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
        entry = isotopeDistribution.getIsotopeCorrectionNum(456.47854, 10, 1, 4, parentPeakList);
        assertEquals(-2, entry.isotopeCorrectionNum);
        assertTrue(entry.pearsonCorrelationCoefficient > 0.7);
    }
}