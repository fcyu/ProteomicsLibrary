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

public class SparseBooleanVectorTest {

    private static SparseBooleanVector sparseBooleanVector;

    @Before
    public void setUp() {
        Set<Integer> sparseSet = new HashSet<>();
        sparseSet.add(1);
        sparseSet.add(3);
        sparseSet.add(5);
        sparseSet.add(10);
        sparseSet.add(11);
        sparseBooleanVector = new SparseBooleanVector(sparseSet);
    }

    @Test
    public void norm2square() {
        assertEquals(5, sparseBooleanVector.norm2square(), 1e-6);
    }

    @Test
    public void dot() {
        Set<Integer> sparseSet = new HashSet<>();
        sparseSet.add(1);
        sparseSet.add(4);
        sparseSet.add(5);
        sparseSet.add(9);
        sparseSet.add(15);
        SparseBooleanVector sparseBooleanVector2 = new SparseBooleanVector(sparseSet);
        assertEquals(2, sparseBooleanVector.dot(sparseBooleanVector2), 1e-6);
        Map<Integer, Double> sparseMap = new HashMap<>();
        sparseMap.put(1, 2d);
        sparseMap.put(4, 3d);
        sparseMap.put(6, 1d);
        sparseMap.put(10, 5d);
        sparseMap.put(11, -3d);
        SparseVector sparseVector = new SparseVector(sparseMap);
        assertEquals(4, sparseBooleanVector.dot(sparseVector), 1e-6);
    }

    @Test
    public void fastDot() {
        Set<Integer> sparseSet = new HashSet<>();
        sparseSet.add(1);
        sparseSet.add(4);
        sparseSet.add(5);
        sparseSet.add(9);
        sparseSet.add(15);
        SparseBooleanVector sparseBooleanVector2 = new SparseBooleanVector(sparseSet);
        Map<Integer, Double> sparseMap = new HashMap<>();
        sparseMap.put(1, 2d);
        sparseMap.put(4, 3d);
        sparseMap.put(6, 1d);
        sparseMap.put(10, 5d);
        sparseMap.put(11, -3d);
        SparseVector sparseVector = new SparseVector(sparseMap);
        System.out.println("Checkng if the sparseBooleanVector2 being changed:");
        System.out.println(String.format(Locale.US, "Original: %s", sparseBooleanVector2.toString()));
        assertEquals(5, sparseBooleanVector2.fastDot(sparseVector), 1e-6);
        System.out.println(String.format(Locale.US, "New: %s\n", sparseBooleanVector2.toString())); // visually check if the sparseBooleanVector2 being changed.
    }

    @Test
    public void deepCopy() {
        SparseBooleanVector sparseBooleanVector2 = sparseBooleanVector.deepCopy();
        System.out.println("Checkng if the sparseBooleanVector being changed:");
        System.out.println(String.format(Locale.US, "Original: %s", sparseBooleanVector.toString()));
        sparseBooleanVector2.delete(1);
        assertNotEquals(sparseBooleanVector, sparseBooleanVector2);
        System.out.println(String.format(Locale.US, "New: %s", sparseBooleanVector.toString()));
    }

    @Test
    public void isZero() {
        assertFalse(sparseBooleanVector.isZero(1));
        assertTrue(sparseBooleanVector.isZero(2));
        assertFalse(sparseBooleanVector.isZero(10));
    }

    @Test
    public void delete() {
        SparseBooleanVector sparseBooleanVector2 = sparseBooleanVector.deepCopy();
        System.out.println(String.format(Locale.US, "After deleting 3, sparseBooleanVector = %s.", sparseBooleanVector2.toString()));
        sparseBooleanVector2.delete(3);
        assertTrue(sparseBooleanVector2.isZero(3));
        System.out.println(String.format(Locale.US, "After deleting 3, sparseBooleanVector = %s.\n", sparseBooleanVector2.toString()));
    }

    @Test
    public void getNonZeroNum() {
        assertEquals(5, sparseBooleanVector.getNonZeroNum());
    }

    @Test
    public void getNonZeroIdxes() {
        Integer[] results = sparseBooleanVector.getNonZeroIdxes();
        Arrays.sort(results);
        Integer[] groundTruth = new Integer[]{1, 3, 5, 10, 11};
        assertArrayEquals(groundTruth, results);
    }
}