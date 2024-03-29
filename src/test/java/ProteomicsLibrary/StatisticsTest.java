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

import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;

public class StatisticsTest {

    @Test
    public void calMedian() throws Exception {
        List<Double> input = new LinkedList<>();
        input.add(1.1);
        input.add(2.3);
        input.add(5.6);
        input.add(0.1);
        input.add(1.2);
        assertEquals(1.2, Statistics.calMedian(input), 1e-6);

        input.clear();
        input.add(1.1);
        input.add(2.3);
        input.add(5.6);
        input.add(0.1);
        input.add(1.1);
        assertEquals(1.1, Statistics.calMedian(input), 1e-6);

        input.clear();
        input.add(2.3);
        input.add(5.6);
        input.add(0.1);
        input.add(1.1);
        assertEquals(1.7, Statistics.calMedian(input), 1e-6);
    }

    @Test(expected = Exception.class)
    public void calMedian2() throws Exception {
        double temp = Statistics.calMedian(new LinkedList<>());
    }

    @Test
    public void calMean() throws Exception {
        double[] input = new double[]{1.1, 1.2, 3.4};
        assertEquals(1.9, Statistics.calMean(input), 1e-6);
        input = new double[]{1.1, 1.2, 0};
        assertEquals(0.766667, Statistics.calMean(input), 0.0001);
    }

    @Test(expected = Exception.class)
    public void calMean2() throws Exception {
        double temp = Statistics.calMean(new double[0]);
    }

    @Test
    public void calMean1() throws Exception {
        List<Double> input = new ArrayList<>();
        input.add(1.1);
        input.add(1.2);
        input.add(3.4);
        assertEquals(1.9, Statistics.calMean(input), 1e-6);

        input.clear();
        input.add(1.1);
        input.add(1.2);
        input.add(0d);
        assertEquals(0.766667, Statistics.calMean(input), 0.0001);
    }

    @Test(expected = Exception.class)
    public void calMean3() throws Exception {
        double temp = Statistics.calMean(new LinkedList<>());
    }

    @Test
    public void calSd() throws Exception {
        List<Double> input = new ArrayList<>();
        input.add(1.1);
        input.add(1.2);
        input.add(3.4);
        assertEquals(1.061445555, Statistics.calSd(input, 1.9), 1e-6);
    }

    @Test(expected = Exception.class)
    public void calSd2() throws Exception {
        List<Double> input = new ArrayList<>();
        input.add(1.1);
        input.add(1.2);
        assertEquals(0.05, Statistics.calSd(input, 1.15), 0.0001);
    }

    @Test
    public void getPercentile() throws Exception {
        List<Double> input = new ArrayList<>();
        input.add(3.0);
        input.add(6.0);
        input.add(7.5);
        input.add(8.4);
        input.add(9.2);
        input.add(10.0);
        input.add(12.8);
        input.add(20.3);
        input.add(12.0);
        assertEquals(12, Statistics.getPercentile(input, 75), 0.0001);
        assertEquals(12.8, Statistics.getPercentile(input, 80), 0.0001);
        assertEquals(6, Statistics.getPercentile(input, 12), 0.0001);
        assertEquals(7.5, Statistics.getPercentile(input, 25), 0.0001);
    }

    @Test(expected = Exception.class)
    public void getPercentile2() throws Exception {
        double temp = Statistics.getPercentile(new LinkedList<>(), 75);
    }

    @Test
    public void tTestTwoSides() throws Exception {
        assertEquals(0.0571, Statistics.tTestTwoSides(2, 2.9, 0, 10), 1e-4);
        assertEquals(0.3688, Statistics.tTestTwoSides(-1, 6, 0, 30), 1e-4);
    }

    @Test(expected = Exception.class)
    public void tTestTwoSides2() throws Exception {
        assertEquals(0.5080, Statistics.tTestTwoSides(2, 2.9, 0, 2), 1e-4);
    }

    @Test
    public void BHFDR() {
        List<Double> pValueList = new ArrayList<>();
        pValueList.add(0d);
        pValueList.add(0.01);
        pValueList.add(0.03);
        pValueList.add(0.005);
        pValueList.add(0.2);
        pValueList.add(0.9);
        pValueList.add(0.0003);
        Map<Double, Double> results = Statistics.BHFDR(pValueList);
        Map<Double, Double> groundTruth = new HashMap<>();
        groundTruth.put(0d, 0d);
        groundTruth.put(0.01, 0.0175);
        groundTruth.put(0.03, 0.042);
        groundTruth.put(0.005, 0.01166667);
        groundTruth.put(0.2, 0.2333333);
        groundTruth.put(0.9, 0.9);
        groundTruth.put(0.0003, 0.00105);
        assertEquals(groundTruth.size(), results.size());
        for (double pValue : groundTruth.keySet()) {
            assertEquals(groundTruth.get(pValue), results.get(pValue), 1e-4);
        }
    }

    @Test
    public void scaleAndCalPearsonCorrelationCoefficient() throws Exception {
        double[] input1 = new double[]{1, 2, 3, 4, 5, 6};
        double[] input2 = new double[]{2, 4, 6, 8, 10, 12};
        assertEquals(1, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
        input1 = new double[]{1, 2, 3, 4, 5, 6};
        input2 = new double[]{-2, -4, -6, -8, -10, -12};
        assertEquals(-1, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
        input1 = new double[]{1, 2, 3, 4, 5, 6};
        input2 = new double[]{2, 3, 5, 6, 7, 8};
        assertEquals(0.9922, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
        input1 = new double[]{1, 2, 3, 4, 5, 6};
        input2 = new double[]{1, 3, 16, 2, 4, 9};
        assertEquals(0.2716, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
        input1 = new double[]{1, 2, 3, 4, 5, 6};
        input2 = new double[]{1, 1, 1, 1, 1, 1};
        assertEquals(0, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
        input1 = new double[]{1, 1, 1, 2, 1, 1};
        input2 = new double[]{1, 2, 3, 4, 5, 6};
        assertEquals(0.1309, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-4);
    }

    @Test(expected = Exception.class)
    public void scaleAndCalPearsonCorrelationCoefficient2() throws Exception {
        double[] input1 = new double[]{1, 2, 3, 4, 5, 6};
        double[] input2 = new double[]{2, 4, 6, 8, 10};
        assertEquals(1, Statistics.calPearsonCorrelationCoefficient(input1, input2), 1e-6);
    }

    @Test
    public void linearInterpolation() throws Exception {
        TreeMap<Integer, Double> input = new TreeMap<>();
        input.put(1, 2.3);
        input.put(3, 4.5);
        input.put(7, 2.1);
        double[] results = Statistics.linearInterpolation(input, 0, 8);
        double[] groundTruth = new double[]{1.1999999999999997, 2.3, 3.4, 4.5, 3.9, 3.3, 2.7, 2.1, 1.5};
        assertEquals(groundTruth.length, results.length);
        for (int i = 0; i < groundTruth.length; ++i) {
            assertEquals(groundTruth[i], results[i], 1e-4);
        }

        input.clear();
        input.put(2, 1d);
        input.put(4, -3d);
        input.put(7, 4d);
        results = Statistics.linearInterpolation(input, -1, 9);
        groundTruth = new double[]{7, 5, 3, 1, -1, -3, -0.6666666666666665, 1.666666666666667, 4, 6.333333333333334, 8.666666666666666};
        assertEquals(groundTruth.length, results.length);
        for (int i = 0; i < groundTruth.length; ++i) {
            assertEquals(groundTruth[i], results[i], 1e-4);
        }
    }

    @Test(expected = Exception.class)
    public void linearInterpolation2() throws Exception {
        TreeMap<Integer, Double> input = new TreeMap<>();
        input.put(1, 2.3);
        double[] results = Statistics.linearInterpolation(input, 0, 8);
        double[] groundTruth = new double[]{2.3};
        assertEquals(groundTruth.length, results.length);
        for (int i = 0; i < groundTruth.length; ++i) {
            assertEquals(groundTruth[i], results[i], 1e-4);
        }
    }

    @Test
    public void fisherExactTest() throws Exception { // The ground truch is from R fisher.test()
        assertEquals(0.6573, Statistics.fisherExactTest(2, 3, 4, 5, -1), 0.0001);
        assertEquals(1, Statistics.fisherExactTest(2, 3, 4, 5, 0), 0.0001);
        assertEquals(0.7622, Statistics.fisherExactTest(4, 5, 2, 3, -1), 0.0001);
        assertEquals(1, Statistics.fisherExactTest(4, 5, 2, 3, 0), 0.0001);
        assertEquals(0.0005241, Statistics.fisherExactTest(10, 3, 4, 20, 1), 0.0001);
        assertEquals(0.0008452, Statistics.fisherExactTest(10, 3, 4, 20, 0), 0.0001);
        assertEquals(1, Statistics.fisherExactTest(3, 10, 20, 4, 1), 0.0001);
        assertEquals(0.0008452, Statistics.fisherExactTest(3, 10, 20, 4, 0), 0.0001);
        assertEquals(0.7219, Statistics.fisherExactTest(20, 30, 40, 50, 0), 0.0001);
        assertEquals(0.7408, Statistics.fisherExactTest(53, 3405294, 34, 2034187, 0), 0.0001);
    }

    @Test(expected = Exception.class)
    public void fisherExactTest2() throws Exception {
        assertEquals(0.6573, Statistics.fisherExactTest(2, 3, -4, 5, -1), 0.0001);
    }

    @Test
    public void contingencyTableChiSquareTest() {
        assertEquals(0.6106, Statistics.contingencyTableChiSquareTest(20, 30, 40, 50), 0.0001);
        assertEquals(0.6106, Statistics.contingencyTableChiSquareTest(40, 50, 20, 30), 0.0001);
        assertEquals(0.0019, Statistics.contingencyTableChiSquareTest(10, 25, 28, 16), 0.0001);
    }
}