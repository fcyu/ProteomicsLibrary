package ProteomicsLibrary;

import org.bytedeco.javacpp.gsl;

import java.util.*;

public class Statistics {

    public static double calMedian(Collection<Double> inputList) {
        if (inputList.isEmpty()) {
            return 0;
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

    public static double calMean(double[] input) {
        if (input.length > 0) {
            double mean = 0;
            for (double v : input) {
                mean += v;
            }
            return mean / input.length;
        } else {
            return 0;
        }
    }

    public static double calMean(Collection<Double> inputList) {
        if (inputList.isEmpty()) {
            return 0;
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

    public static double tTestTwoSides(double mean, double sd, double mu, int num) throws Exception { // two-sided t-test
        if (num > 2) {
            double t = (mean - mu) * Math.sqrt(num) / sd;
            if (t >= 0) {
                return 2 * (1 - gsl.gsl_cdf_tdist_P(t, num - 1));
            } else {
                return 2 * gsl.gsl_cdf_tdist_P(t, num - 1);
            }
        } else {
            throw new Exception(String.format(Locale.US, "The number of elements (%d) is smaller than 3.", num));
        }
    }

    public static Map<Double, Double> BHFDR(Collection<Double> pValueList) {
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
}
