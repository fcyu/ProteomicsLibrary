package ProteomicsLibrary;

import ProteomicsLibrary.Types.Coordinate;

import java.util.*;

public class Score {

    public static double calIonFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> expPL, double ms2Tolerance) {
        int matchedPeakNum = 0;
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        int totalIonNum = ionMatrix[0].length * maxRow;
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : expPL.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++matchedPeakNum;
                        break;
                    }
                }
            }
        }

        return (double) matchedPeakNum / (double) totalIonNum;
    }

    public static double calMatchedHighestIntensityFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> expPL, double ms2Tolerance) {
        int matchedPeakNum = 0;
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        int totalIonNum = ionMatrix[0].length * maxRow;
        Double[] intensityArray = expPL.values().toArray(new Double[0]);
        Arrays.sort(intensityArray, Collections.reverseOrder());
        double intensityT = 0;
        if (totalIonNum < intensityArray.length) {
            intensityT = intensityArray[totalIonNum];
        }
        int matchedHighestPeakNum = 0;
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : expPL.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        if (expPL.get(mz) > intensityT) {
                            ++matchedHighestPeakNum;
                        }
                        ++matchedPeakNum;
                        break;
                    }
                }
            }
        }

        if (matchedPeakNum > 0) {
            return (double) matchedHighestPeakNum / (double) matchedPeakNum;
        } else {
            return 0;
        }
    }

    public static double calExplainedAAFraction(double[][] ionMatrix, int precursorCharge, Map<Double, Double> expPL, double ms2Tolerance) {
        Set<Integer> matchedIdxSet = new HashSet<>(); // 0 if the mz = 0 peak; 1 is the peak generated by b1 ion...
        matchedIdxSet.add(0); // Add N-term automatically.
        matchedIdxSet.add(ionMatrix[0].length); // Add C-term automatically because there is no need to observe the peak so that the last amino acid can be inferred.
        int maxRow = Math.max(2, Math.min(ionMatrix.length, 2 * (precursorCharge - 1)));
        for (int i = 0; i < maxRow; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : expPL.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        if (i % 2 == 0) {
                            matchedIdxSet.add(j + 1);
                        } else {
                            if (j > 0) {
                                matchedIdxSet.add(j);
                            }
                        }
                        break;
                    }
                }
            }
        }

        // calculate explained AA num
        Integer[] matchedIdxArray = matchedIdxSet.toArray(new Integer[0]);
        Arrays.sort(matchedIdxArray);
        int explainedAaNum = 0;
        if (matchedIdxArray.length > 1) {
            for (int i = 0; i < matchedIdxArray.length - 1; ++i) {
                if (matchedIdxArray[i + 1] - matchedIdxArray[i] == 1) {
                    ++explainedAaNum;
                }
            }
        }
        return (double) explainedAaNum / (double) ionMatrix[0].length;
    }

    public static double calBinomialScorePValue(TreeMap<Double, Double> expPL, int topN, Binomial binomial, int localMaxMs2Charge, double[][] ionMatrix, double ms2Tolerance, int peptideLengthWithNC) throws Exception { // calculate a p-value as Andromeda.s
        double finalPValue = 2;
        for (int localTopN = 1; localTopN <= topN; ++localTopN) {
            TreeMap<Double, Double> localPLMap = PrepareSpectrum.topNStyleNormalization(expPL, localTopN);
            double pValue = calBinomialScorePValueSub(localPLMap, localTopN, binomial, localMaxMs2Charge, ionMatrix, ms2Tolerance, peptideLengthWithNC);
            if (pValue < finalPValue) {
                finalPValue = pValue;
            }
        }
        return finalPValue;
    }

    public static double calBinomialScorePValueSub(TreeMap<Double, Double> localPLMap, int localTopN, Binomial binomial, int localMaxMs2Charge, double[][] ionMatrix, double ms2Tolerance, int peptideLengthWithNC) throws Exception {
        int matchedIonNum = getMatchedIonNum(localPLMap, localMaxMs2Charge, ionMatrix, ms2Tolerance);
        return binomial.calProbLargerThanOrEqualTo((peptideLengthWithNC - 2) * 2 * localMaxMs2Charge, matchedIonNum, localTopN * 0.01); // todo: the p is not accurate, but we don't have a perfect solution.
    }

    public static double calAScore(TreeMap<Double, Double> expPL, int topN, Binomial binomial, TreeMap<Coordinate, Double> varPTMMap1, double[][] ionMatrix1, TreeMap<Coordinate, Double> varPTMMap2, double[][] ionMatrix2, double ms2Tolerance, int peptideLengthWithNC) throws Exception {
        double finalAScore = -9999;
        for (int localTopN = 1; localTopN <= topN; ++localTopN) {
            TreeMap<Double, Double> localPLMap = PrepareSpectrum.topNStyleNormalization(expPL, localTopN);
            double aScore = calAScoreSub(localPLMap, localTopN, binomial, varPTMMap1, ionMatrix1, varPTMMap2, ionMatrix2, ms2Tolerance, peptideLengthWithNC);
            if (aScore > finalAScore) {
                finalAScore = aScore;
            }
        }
       return finalAScore;
    }

    public static double calAScoreSub(TreeMap<Double, Double> localPLMap, int localTopN, Binomial binomial, TreeMap<Coordinate, Double> varPTMMap1, double[][] ionMatrix1, TreeMap<Coordinate, Double> varPTMMap2, double[][] ionMatrix2, double ms2Tolerance, int peptideLengthWithNC) throws Exception {
        TreeSet<Integer> totalAffectedBIonSet = new TreeSet<>();
        TreeSet<Integer> totalAffectedYIonSet = new TreeSet<>();
        if (varPTMMap2 == null) {
            getAffectedIonSet(varPTMMap1, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet); // don't delete two most outside ions because they are used to fix the location when there is no second peptide
            Set<String> topMatchedPeakSet = getMatchedIonSet(ionMatrix1, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            return  -10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), topMatchedPeakSet.size(), localTopN * 0.01)); // todo: the p is not accurate, but we don't have a perfect solution.
        } else {
            getAffectedIonSet(varPTMMap1, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet);
            getAffectedIonSet(varPTMMap2, peptideLengthWithNC - 2, totalAffectedBIonSet, totalAffectedYIonSet);

            // delete two most outside ions because they are not effected by the different PTM locations.
            if (!totalAffectedBIonSet.contains(1)) {
                totalAffectedBIonSet.pollFirst();
            }
            if (!totalAffectedBIonSet.contains(peptideLengthWithNC - 3)) {
                totalAffectedBIonSet.pollLast();
            }
            if (!totalAffectedYIonSet.contains(1)) {
                totalAffectedYIonSet.pollFirst();
            }
            if (!totalAffectedYIonSet.contains(peptideLengthWithNC - 3)) {
                totalAffectedYIonSet.pollLast();
            }

            Set<String> topMatchedPeakSet = getMatchedIonSet(ionMatrix1, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            Set<String> secondMatchedPeakSet = getMatchedIonSet(ionMatrix2, localPLMap, ms2Tolerance, totalAffectedBIonSet, totalAffectedYIonSet);
            return -10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), topMatchedPeakSet.size(), localTopN * 0.01)) + 10 * Math.log10(binomial.calProbLargerThanOrEqualTo(totalAffectedBIonSet.size() + totalAffectedYIonSet.size(), secondMatchedPeakSet.size(), localTopN * 0.01)); // todo: the p is not accurate, but we don't have a perfect solution.
        }
    }

    public static int getMatchedIonNum(TreeMap<Double, Double> expPL, int localMaxMs2Charge, double[][] ionMatrix, double ms2Tolerance) {
        int K1 = 0;
        for (int i = 0; i < localMaxMs2Charge * 2; ++i) {
            for (int j = 0; j < ionMatrix[0].length; ++j) {
                for (double mz : expPL.keySet()) {
                    if (Math.abs(mz - ionMatrix[i][j]) <= ms2Tolerance) {
                        ++K1;
                        break;
                    }
                }
            }
        }
        return K1;
    }

    private static void getAffectedIonSet(TreeMap<Coordinate, Double> varPtmMap, int peptideLength, Set<Integer> affectedBIonSet, Set<Integer> affectedYIonSet) {
        for (Coordinate co : varPtmMap.keySet()) {
            if (co.x == 0 || co.x == 1) {
                affectedBIonSet.add(1);
                affectedYIonSet.add(peptideLength - 1);
            } else if (co.x == peptideLength + 1 || co.x == peptideLength) {
                affectedBIonSet.add(peptideLength - 1);
                affectedYIonSet.add(1);
            } else {
                affectedBIonSet.add( co.x - 1);
                affectedBIonSet.add(co.x);
                affectedYIonSet.add( peptideLength - co.x);
                affectedYIonSet.add(peptideLength - co.x + 1);
            }
        }
    }

    private static Set<String> getMatchedIonSet(double[][] ionMatrix, TreeMap<Double, Double> localPLMap, double ms2Tolerance, Set<Integer> affectedBIonSet, Set<Integer> affectedYIonSet) {
        Set<String> matchedIonSet = new HashSet<>();
        for (int ion : affectedBIonSet) {
            double theoMass = ionMatrix[0][ion - 1];
            for (double mz : localPLMap.keySet()) {
                if (Math.abs(mz - theoMass) <= ms2Tolerance) {
                    matchedIonSet.add(String.format(Locale.US, "b%d", ion));
                    break;
                }
            }
        }

        for (int ion : affectedYIonSet) {
            double theoMass = ionMatrix[1][ionMatrix[0].length - ion];
            for (double mz : localPLMap.keySet()) {
                if (Math.abs(mz - theoMass) <= ms2Tolerance) {
                    matchedIonSet.add(String.format(Locale.US, "y%d", ion));
                    break;
                }
            }
        }

        return matchedIonSet;
    }
}
