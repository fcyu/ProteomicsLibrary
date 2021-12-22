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

package ProteomicsLibrary.Types;

import org.junit.Before;
import org.junit.Test;

import java.util.*;

import static org.junit.Assert.*;

public class SparseVectorTest {

    private static SparseVector sparseVector;

    @Before
    public void setUp() {
        Map<Integer, Double> sparseMap = new HashMap<>();
        sparseMap.put(1, 2d);
        sparseMap.put(4, 3d);
        sparseMap.put(6, 1d);
        sparseMap.put(10, 5d);
        sparseMap.put(11, -3d);
        sparseMap.put(14, 1d);
        sparseVector = new SparseVector(sparseMap);
    }

    @Test
    public void add() {
        sparseVector.add(1, 3);
        sparseVector.add(2, 1);
        sparseVector.add(6, -3);
        sparseVector.add(10, -5);
        assertEquals(5, sparseVector.get(1), 1e-6);
        assertEquals(1, sparseVector.get(2), 1e-6);
        assertEquals(3, sparseVector.get(4), 1e-6);
        assertEquals(-2, sparseVector.get(6), 1e-6);
        assertEquals(0, sparseVector.get(10), 1e-6);
    }

    @Test
    public void put() {
        sparseVector.put(2, 4);
        sparseVector.put(4, 6);
        assertEquals(4, sparseVector.get(2), 1e-6);
        assertEquals(6, sparseVector.get(4), 1e-6);
        assertNotEquals(3, sparseVector.get(4), 1e-6);
    }

    @Test
    public void get() {
        assertEquals(2, sparseVector.get(1), 1e-6);
        assertEquals(0, sparseVector.get(2), 1e-6);
    }

    @Test
    public void idxSet() {
        Integer[] groundTruth = new Integer[]{1, 4, 6, 10, 11, 14};
        Integer[] results = sparseVector.idxSet().toArray(new Integer[0]);
        Arrays.sort(results);
        assertArrayEquals(groundTruth, results);
    }

    @Test
    public void getValues() {
        Double[] groundTruth = new Double[]{2d, 3d, 1d, 5d, -3d, 1d};
        Arrays.sort(groundTruth);
        Double[] results = sparseVector.getValues();
        Arrays.sort(results);
        assertArrayEquals(groundTruth, results);
    }

    @Test
    public void getMaxValue() {
        assertEquals(5, sparseVector.getMaxValue(), 1e-6);
    }

    @Test
    public void getMinValue() {
        assertEquals(-3, sparseVector.getMinValue(), 1e-6);
    }

    @Test
    public void norm2square() {
        assertEquals(49, sparseVector.norm2square(), 1e-6);
    }

    @Test
    public void dot() {
        Map<Integer, Double> sparseMap = new HashMap<>();
        sparseMap.put(1, -3d);
        sparseMap.put(6, 0.5d);
        sparseMap.put(12, -1d);
        sparseMap.put(14, 4.5d);
        SparseVector sparseVector2 = new SparseVector(sparseMap);
        assertEquals(-1, sparseVector.dot(sparseVector2), 1e-6);
    }

    @Test
    public void isEmpty() {
        assertFalse(sparseVector.isEmpty());
    }

    @Test
    public void getVectorMap() {
        Map<Integer, Double> groundTruth = new HashMap<>();
        groundTruth.put(1, 2d);
        groundTruth.put(4, 3d);
        groundTruth.put(6, 1d);
        groundTruth.put(10, 5d);
        groundTruth.put(11, -3d);
        groundTruth.put(14, 1d);
        Map<Integer, Double> results = sparseVector.getVectorMap();
        assertEquals(groundTruth.size(), results.size());
        for (int i : groundTruth.keySet()) {
            assertEquals(groundTruth.get(i), results.get(i), 1e-6);
        }
    }

    @Test
    public void getNonzeroIdx() {
        Integer[] groundTrugh = new Integer[]{1, 4, 6, 10, 11, 14};
        Integer[] results = sparseVector.getNonzeroIdx().toArray(new Integer[0]);
        Arrays.sort(results);
        assertArrayEquals(groundTrugh, results);
    }

    @Test
    public void isNonzero() {
        assertTrue(sparseVector.isNonzero(1));
        assertFalse(sparseVector.isNonzero(2));
    }
}