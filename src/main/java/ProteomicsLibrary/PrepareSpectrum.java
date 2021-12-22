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

import ProteomicsLibrary.Types.SparseVector;

import java.util.*;

public class PrepareSpectrum {

    private static final double defaultIntensity = 1; // DO NOT change. Otherwise, change the whole project accordingly.
    private static final int xcorrOffset = 75;
    private static final double removePrecursorPeakTolerance = 1.5; // this equals the isolation window.

    private final MassTool massTool;

    public PrepareSpectrum(MassTool massTool) {
        this.massTool = massTool;
    }

    public TreeMap<Double, Double> preSpectrumTopNStyle (Map<Double, Double> inputPL, double precursorMass, int precursorCharge, double minClear, double maxClear, int topN) {
        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge, minClear, maxClear);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return topNStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass))), topN);
        }
    }

    public TreeMap<Double, Double> preSpectrumTopNStyle (Map<Double, Double> inputPL, double precursorMass, int topN) {
        TreeMap<Double, Double> outputPL = new TreeMap<>(inputPL);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new TreeMap<>();
        } else {
            return topNStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass))), topN);
        }
    }

    public SparseVector preSpectrumCometStyle (Map<Double, Double> inputPL, double precursorMass, int precursorCharge, double minClear, double maxClear, boolean flankingPeak) {
        TreeMap<Double, Double> outputPL = removeCertainPeaks(inputPL, precursorMass, precursorCharge, minClear, maxClear);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new SparseVector();
        } else {
            return prepareXCorr(cometStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass)))), flankingPeak);
        }
    }

    public SparseVector preSpectrumCometStyle (Map<Double, Double> inputPL, double precursorMass, boolean flankingPeak) {
        TreeMap<Double, Double> outputPL = new TreeMap<>(inputPL);
        if (outputPL.subMap(0d, precursorMass).isEmpty()) {
            return new SparseVector();
        } else {
            return prepareXCorr(cometStyleNormalization(sqrtPL(new TreeMap<>(outputPL.subMap(0d, precursorMass)))), flankingPeak);
        }
    }

    public SparseVector prepareXCorr(TreeMap<Double, Double> plMap, boolean flankingPeaks) {
        if (plMap.isEmpty()) {
            return new SparseVector();
        } else {
            // digitize peak list
            double[] plArray = new double[massTool.mzToBin(plMap.lastKey()) + 1];
            for (double mz : plMap.keySet()) {
                if (Math.abs(plMap.get(mz)) > 1e-6) {
                    int idx = massTool.mzToBin(mz);
                    plArray[idx] = Math.max(plMap.get(mz), plArray[idx]);
                }
            }
            return prepareXCorr(plArray, flankingPeaks);
        }
    }

    public SparseVector prepareXCorr(double[] plArray, boolean flankingPeaks) {
        SparseVector xcorrPL = new SparseVector();
        int offsetRange = 2 * xcorrOffset + 1;
        double factor = 1 / (double) (offsetRange - 1); // caution: 1/150 rather than 1/151
        double mySum = 0;
        for (int i = 0; i < xcorrOffset; ++i) {
            mySum += plArray[i];
        }

        double[] tempArray = new double[plArray.length];
        for (int i = xcorrOffset; i < plArray.length + xcorrOffset; ++i) {
            if (i < plArray.length) {
                mySum += plArray[i];
            }
            if (i >= offsetRange) {
                mySum -= plArray[i - offsetRange];
            }
            tempArray[i - xcorrOffset] = (mySum - plArray[i - xcorrOffset]) * factor;
        }

        for (int i = 1; i < plArray.length; ++i) {
            double temp = plArray[i] - tempArray[i];
            if (flankingPeaks) {
                temp += (plArray[i - 1] - tempArray[i - 1]) * 0.5;
                if (i + 1 < plArray.length) {
                    temp += (plArray[i + 1] - tempArray[i + 1]) * 0.5;
                }
            }
            if (Math.abs(temp) > 1e-6) {
                xcorrPL.put(i, temp);
            }
        }

        return xcorrPL;
    }

    public SparseVector digitizePL(TreeMap<Double, Double> plMap) {
        SparseVector digitizedPL = new SparseVector();
        for (double mz : plMap.keySet()) {
            int idx = massTool.mzToBin(mz);
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                digitizedPL.put(idx, Math.max(plMap.get(mz), digitizedPL.get(idx)));
            }
        }
        return digitizedPL;
    }

    public static TreeMap<Double, Double> topNStyleNormalization(TreeMap<Double, Double> inputPL, int localTopN) {
        if (inputPL.isEmpty()) {
            return new TreeMap<>();
        } else {
            // select top N in each 100 Da
            TreeMap<Double, Double> outputPL = new TreeMap<>();
            double minMz = inputPL.firstKey();
            double maxMz = inputPL.lastKey();
            double leftMz = minMz;
            while (leftMz < maxMz) {
                // find the max intensity in each window
                double rightMz = Math.min(leftMz + 100, maxMz);
                NavigableMap<Double, Double> subMap;
                if (rightMz < maxMz) {
                    subMap = inputPL.subMap(leftMz, true, rightMz, false);
                } else {
                    subMap = inputPL.subMap(leftMz, true, rightMz, true);
                }
                if (!subMap.isEmpty()) {
                    Double[] intensityArray = subMap.values().toArray(new Double[0]);
                    Arrays.sort(intensityArray, Comparator.reverseOrder());
                    double temp1 = defaultIntensity / intensityArray[0];
                    double temp2 = subMap.size() > localTopN ? intensityArray[localTopN] : 0;
                    for (double mz : subMap.keySet()) {
                        if (subMap.get(mz) > temp2) {
                            outputPL.put(mz, subMap.get(mz) * temp1);
                        }
                    }
                }
                leftMz = rightMz;
            }

            return outputPL;
        }
    }

    private double[] cometStyleNormalization(TreeMap<Double, Double> plMap) {
        // digitize peak list
        double[] plArray = new double[massTool.mzToBin(plMap.lastKey()) + 1];
        for (double mz : plMap.keySet()) {
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                int idx = massTool.mzToBin(mz);
                plArray[idx] = Math.max(plMap.get(mz), plArray[idx]);
            }
        }

        // find the max intensity
        double maxIntensity = 0;
        for (double intensity : plArray) {
            if (intensity > maxIntensity) {
                maxIntensity = intensity;
            }
        }

        // normalize the peaks
        double[] normalizedPlArray = new double[plArray.length];
        int windowSize = (plArray.length / 10) + 1;
        for (int i = 0; i < 10; ++i) {
            // find the max intensity in each window
            double maxWindowIntensity = 0;
            for (int j = 0; j < windowSize; ++j) {
                int idx = i * windowSize + j;
                if (idx < plArray.length) {
                    if (plArray[idx] > maxWindowIntensity) {
                        maxWindowIntensity = plArray[idx];
                    }
                }
            }

            if (maxWindowIntensity > 0) {
                double temp1 = 50 / maxWindowIntensity;
                double temp2 = 0.05 * maxIntensity; // caution: Xolik does not have this
                for (int j = 0; j < windowSize; ++j) {
                    int idx = i * windowSize + j;
                    if (idx < plArray.length) {
                        if (plArray[idx] > temp2) {
                            normalizedPlArray[idx] = plArray[idx] * temp1;
                        }
                    }
                }
            }
        }

        return normalizedPlArray;
    }

    private static TreeMap<Double, Double> removeCertainPeaks(Map<Double, Double> peakMap, double precursorMass, int precursorCharge, double minClear, double maxClear) {
        TreeMap<Double, Double> mzIntensityMap = new TreeMap<>();
        double precursorMz = precursorMass / precursorCharge + MassTool.PROTON;
        for (double mz : peakMap.keySet()) {
            if (((mz < minClear) || (mz > maxClear)) && (mz > 50)) {
                if ((peakMap.get(mz) > 1e-6) && (Math.abs(peakMap.get(mz) - precursorMz) > removePrecursorPeakTolerance)) {
                    mzIntensityMap.put(mz, peakMap.get(mz));
                }
            }
        }

        return mzIntensityMap;
    }

    private static TreeMap<Double, Double> sqrtPL(TreeMap<Double, Double> plMap) {
        // sqrt the intensity and find the highest intensity.
        TreeMap<Double, Double> sqrtPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            if (plMap.get(mz) > 1e-6) {
                sqrtPlMap.put(mz, Math.sqrt(plMap.get(mz)));
            }
        }
        return sqrtPlMap;
    }

    private static TreeMap<Double, Double> deNoise(TreeMap<Double, Double> plMap) {
        // denoise
        TreeMap<Double, Double> denoisedPlMap = new TreeMap<>();
        double minMz = plMap.firstKey();
        double maxMz = plMap.lastKey();
        double windowSize = (plMap.lastKey() - plMap.firstKey()) * 0.1 + 1;
        for (int i = 0; i < 10; ++i) {
            double leftMz = Math.min(minMz + i * windowSize, maxMz);
            double rightMz = Math.min(leftMz + windowSize, maxMz);
            NavigableMap<Double, Double> subPlMap;
            if (rightMz < maxMz) {
                subPlMap = plMap.subMap(leftMz, true, rightMz, false);
            } else {
                subPlMap = plMap.subMap(leftMz, true, rightMz, true);
            }

            if (subPlMap.size() > 9) {
                double noiseIntensity = estimateNoiseIntensity(subPlMap);
                for (double mz : subPlMap.keySet()) {
                    if (subPlMap.get(mz) > noiseIntensity) {
                        denoisedPlMap.put(mz, subPlMap.get(mz));
                    }
                }
            } else {
                for (double mz : subPlMap.keySet()) {
                    denoisedPlMap.put(mz, subPlMap.get(mz));
                }
            }
        }

        return denoisedPlMap;
    }

    private static double estimateNoiseIntensity(Map<Double, Double> pl) {
        Set<Double> intensitySet = new HashSet<>();
        for (double intensity : pl.values()) {
            intensitySet.add(intensity);
        }
        Double[] uniqueIntensityVector = intensitySet.toArray(new Double[0]);
        Arrays.sort(uniqueIntensityVector);
        double[] cum = new double[uniqueIntensityVector.length];
        for (int i = 0; i < uniqueIntensityVector.length; ++i) {
            for (double intensity : pl.values()) {
                if (intensity <= uniqueIntensityVector[i]) {
                    ++cum[i];
                }
            }
        }
        double[][] diff = new double[2][uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff[0][i] = cum[i + 1] - cum[i];
            diff[1][i] = uniqueIntensityVector[i + 1] - uniqueIntensityVector[i];
        }
        double[] diff2 = new double[uniqueIntensityVector.length - 1];
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            diff2[i] = diff[0][i] / (diff[1][i] + 1e-6);
        }
        double maxValue = 0;
        int maxIdx = 0;
        for (int i = 0; i < uniqueIntensityVector.length - 1; ++i) {
            if (diff2[i] > maxValue) {
                maxValue = diff2[i];
                maxIdx = i;
            }
        }

        return uniqueIntensityVector[maxIdx];
    }
}
