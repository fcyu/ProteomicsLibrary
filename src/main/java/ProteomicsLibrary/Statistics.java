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

import ProteomicsLibrary.Types.Hypergeometric;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;

import java.util.*;

public class Statistics {

    private static final NormalDistribution normalDistribution = new NormalDistribution(0, 1);

    public static double calMedian(Collection<Double> inputList) throws Exception {
        if (inputList.isEmpty()) {
            throw new Exception("There is no element in the input list.");
        } else {
            Double[] inputArray = inputList.toArray(new Double[0]);
            Arrays.sort(inputArray);
            if (inputArray.length % 2 == 0) {
                return (inputArray[inputArray.length / 2] + inputArray[(inputArray.length / 2) - 1]) / 2;
            } else {
                return inputArray[inputList.size() / 2];
            }
        }
    }

    public static double calMean(double[] input) throws Exception {
        if (input.length > 0) {
            double mean = 0;
            for (double v : input) {
                mean += v;
            }
            return mean / input.length;
        } else {
            throw new Exception("There is no element in the input array.");
        }
    }

    public static double calMean(Collection<Double> inputList) throws Exception {
        if (inputList.isEmpty()) {
            throw new Exception("There is no element in the input list.");
        } else {
            double mean = 0;
            for (double v : inputList) {
                mean += v;
            }
            return mean / inputList.size();
        }
    }

    public static double calSd(Collection<Double> inputList, double mean) throws Exception {
        if (inputList.size() > 2) {
            double sd = 0;
            for (double v : inputList) {
                sd += Math.pow(v - mean, 2);
            }
            return Math.sqrt(sd / inputList.size());
        } else {
            throw new Exception(String.format(Locale.US, "The number of elements (%d) is smaller than 3.", inputList.size()));
        }
    }

    public static double getPercentile(Collection<Double> inputList, double percentile) throws Exception {
        if (inputList.isEmpty()) {
            throw new Exception("There is no element in the input list.");
        } else {
            Double[] inputArray = inputList.toArray(new Double[0]);
            Arrays.sort(inputArray);
            int index = (int) Math.ceil((percentile / 100) * inputArray.length);
            return inputArray[index - 1];
        }
    }

    public static double calMaxQuantProteinPValue(double r, double rMinus1, double r0, double r1) { // todo: add a unit test
        double z;
        if (r >= r0) {
            z = (r - r0) / (r1 - r0);
        } else {
            z = (r0 - r) / (r0 - rMinus1);
        }
        return 1 - normalDistribution.cumulativeProbability(z);
    }

    public static double tTestTwoSides(double mean, double sd, double mu, int num) throws Exception { // two-sided t-test
        if (num > 2) {
            double t = (mean - mu) * Math.sqrt(num) / sd;
            TDistribution tDistribution = new TDistribution(num - 1);
            if (t >= 0) {
                return 2 * (1 - tDistribution.cumulativeProbability(t));
            } else {
                return 2 * tDistribution.cumulativeProbability(t);
            }
        } else {
            throw new Exception(String.format(Locale.US, "The number of elements (%d) is smaller than 3.", num));
        }
    }

    public static Map<Double, Double> BHFDR(Collection<Double> pValueList) {
        if (pValueList.isEmpty()) {
            return new HashMap<>();
        } else if (pValueList.size() == 1) {
            Map<Double, Double> pValueFDRMap = new HashMap<>();
            pValueFDRMap.put(pValueList.iterator().next(), pValueList.iterator().next());
            return pValueFDRMap;
        } else {
            Double[] pValueArray = pValueList.toArray(new Double[0]);
            Arrays.sort(pValueArray);

            double[] tempArray = new double[pValueArray.length];
            for (int i = 0; i < pValueArray.length; ++i) {
                tempArray[i] = pValueArray[i] * pValueArray.length / (i + 1);
            }

            double[] fdrArray = new double[tempArray.length];
            double lastFdr = tempArray[tempArray.length - 1];
            fdrArray[fdrArray.length - 1] = lastFdr;
            for (int i = tempArray.length - 2; i >= 0; --i) {
                if (tempArray[i] >= lastFdr) {
                    fdrArray[i] = lastFdr;
                } else {
                    fdrArray[i] = tempArray[i];
                    lastFdr = tempArray[i];
                }
            }

            Map<Double, Double> pValueFDRMap = new HashMap<>();
            for (int i = 0; i < pValueArray.length; ++i) {
                pValueFDRMap.put(pValueArray[i], fdrArray[i]);
            }

            return pValueFDRMap;
        }
    }

    public static double calPearsonCorrelationCoefficient(double[] input1, double[] input2) throws Exception{
        if (input1.length != input2.length) {
            throw new Exception("Two vectors' lengths are different.");
        }

        // calculate Pearson correlation coefficient.
        double mean1 = calMean(input1);
        double mean2 = calMean(input2);
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        for (int i = 0; i < input1.length; ++i) {
            double c1 = input1[i] - mean1;
            double c2 = input2[i] - mean2;
            temp1 += c1 * c2;
            temp2 += Math.pow(c1, 2);
            temp3 += Math.pow(c2, 2);
        }
        return (temp1 == 0 || temp2 == 0) ? 0 : temp1 / (Math.sqrt(temp2 * temp3));
    }

    public static double[] linearInterpolation(TreeMap<Integer, Double> input, int minX, int maxX) throws Exception {
        if (input.size() > 1) {
            double[] output = new double[maxX - minX + 1];
            for (int i = 0; i < output.length; ++i) {
                Integer start = input.floorKey(minX + i);
                Integer end = input.ceilingKey(minX + i);
                if (start != null && end != null) {
                    if (start.equals(end)) {
                        output[i] = input.get(start);
                    } else {
                        output[i] = input.get(start) + (i + minX - start) * (input.get(end) - input.get(start)) / (end - start);
                    }
                } else if (start == null && end != null) {
                    start = input.ceilingKey(minX + i);
                    end = input.ceilingKey(start + 1);
                    output[i] = input.get(start) - (start - i - minX) * (input.get(end) - input.get(start)) / (end - start);
                } else if (start != null) {
                    output[i] = 0;
                    end = input.floorKey(minX + i);
                    start = input.floorKey(end - 1);
                    output[i] = input.get(end) + (i + minX - end) * (input.get(end) - input.get(start)) / (end - start);
                } else {
                    throw new NullPointerException("Something wrong in the linear interpolation");
                }
            }
            return output;
        } else {
            throw new Exception(String.format(Locale.US, "The input length (%d) is smaller than 2", input.size()));
        }
    }

    public static double fisherExactTest(int a, int b, int c, int d, int alternative) throws Exception {
        // alternative < 0: pCase < pControl; alternative > 0: pCase > pControl; alternative == 0: pCase != pControl
        // a    b // case group
        // c    d // control group
        // The two-sided algorithm are from http://mathworld.wolfram.com/FishersExactTest.html

        if (a < 0 || b < 0 || c < 0 || d < 0) {
            throw new Exception(String.format(Locale.US, "There are counts smaller than 0: a = %d, b = %d, c = %d, d = %d.", a, b, c, d));
        }

        Hypergeometric hypergeometric = new Hypergeometric(a + b + c + d + 1);
        if (alternative > 0) {
            double pSum = hypergeometric.calPMF(a, b, c, d);
            while (b > 0 && c > 0) {
                ++a;
                --b;
                --c;
                ++d;
                pSum += hypergeometric.calPMF(a, b, c, d);
            }
            return pSum;
        } else if (alternative < 0) {
            double pSum = hypergeometric.calPMF(a, b, c, d);
            while (a > 0 && d > 0) {
                --a;
                ++b;
                ++c;
                --d;
                pSum += hypergeometric.calPMF(a, b, c, d);
            }
            return pSum;
        } else {
            int aOriginal = a;
            int bOriginal = b;
            int cOriginal = c;
            int dOriginal = d;

            double pBase = hypergeometric.calPMF(aOriginal, bOriginal, cOriginal, dOriginal);
            double pSum = pBase;

            while (b > 0 && c > 0) {
                ++a;
                --b;
                --c;
                ++d;
                double pTemp = hypergeometric.calPMF(a, b, c, d);
                if (pTemp <= pBase) {
                    pSum += pTemp;
                }
            }

            a = aOriginal;
            b = bOriginal;
            c = cOriginal;
            d = dOriginal;
            while (a > 0 && d > 0) {
                --a;
                ++b;
                ++c;
                --d;
                double pTemp = hypergeometric.calPMF(a, b, c, d);
                if (pTemp <= pBase) {
                    pSum += pTemp;
                }
            }
            return Math.min(1, pSum);
        }
    }

    public static double contingencyTableChiSquareTest(int a, int b, int c, int d) {
        // a    b
        // c    d
        if (a < 5 || b < 5 || c < 5 || d < 5) {
            System.err.println(String.format(Locale.US, "One of the entry is smaller than 5 (a = %d, b = %d, c = %d, d = %d). Try Fisher's exact test instead.", a, b, c, d));
        }

        double total = a + b + c + d;
        double aE = (a + b) * (a + c) / total;
        double bE = (a + b) * (b + d) / total;
        double cE = (c + d) * (a + c) / total;
        double dE = (c + d) * (b + d) / total;
        double chiSquare = (a - aE) * (a - aE) / aE + (b - bE) * (b - bE) / bE + (c - cE) * (c - cE) / cE + (d - dE) * (d - dE) / dE;
        ChiSquaredDistribution chiSquaredDistribution = new ChiSquaredDistribution(1);
        return 1 - chiSquaredDistribution.cumulativeProbability(chiSquare);
    }
}
