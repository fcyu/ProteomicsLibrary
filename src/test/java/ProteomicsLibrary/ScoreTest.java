package ProteomicsLibrary;

import ProteomicsLibrary.Types.Coordinate;
import org.junit.BeforeClass;
import org.junit.Test;

import java.lang.reflect.Method;
import java.util.*;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

public class ScoreTest {

    private static final MassTool massTool = new MassTool(2, "KR", "P", true, 0.01, 1, "N14");
    private static final Binomial binomial = new Binomial(100);

    private static double[][] ionMatrix1;
    private static double[][] ionMatrix2;
    private static TreeMap<Coordinate, Double> varPTMMap1 = new TreeMap<>();
    private static TreeMap<Coordinate, Double> varPTMMap2 = new TreeMap<>();
    private static Map<Double, Double> expPL;

    @BeforeClass
    public static void setUp() {
        ionMatrix1 = massTool.buildIonArray("n(144.102)SDALETLGFLN(0.984)HYQMK(144.102)c", 2);
        ionMatrix2 = massTool.buildIonArray("n(144.102)SDALETLGFLNHYQ(0.984)MK(144.102)c", 2);

        varPTMMap1.put(new Coordinate(0, 1), 114.102);
        varPTMMap1.put(new Coordinate(11, 12), 0.984);
        varPTMMap1.put(new Coordinate(16, 17), 114.102);

        varPTMMap2.put(new Coordinate(0, 1), 114.102);
        varPTMMap2.put(new Coordinate(14, 15), 0.984);
        varPTMMap2.put(new Coordinate(16, 17), 114.102);

        expPL = new TreeMap<>();
        expPL.put(102.0549, 3215.8);
        expPL.put(102.0931, 920.5);
        expPL.put(107.0955, 896.5);
        expPL.put(110.0711, 50239.2);
        expPL.put(111.0739, 2143.3);
        expPL.put(114.0933, 1072.8);
        expPL.put(114.1106, 44235.9);
        expPL.put(115.1076, 34350.4);
        expPL.put(116.1109, 40027.1);
        expPL.put(117.1142, 56281.4);
        expPL.put(118.1171, 1265.8);
        expPL.put(120.0806, 48398.2);
        expPL.put(121.0836, 2527.0);
        expPL.put(136.0756, 22550.5);
        expPL.put(136.7710, 1079.1);
        expPL.put(136.8413, 957.1);
        expPL.put(138.0670, 1352.0);
        expPL.put(145.1084, 12089.7);
        expPL.put(146.2474, 925.4);
        expPL.put(148.3009, 954.1);
        expPL.put(162.4581, 1004.8);
        expPL.put(167.0795, 1400.2);
        expPL.put(171.9518, 1047.1);
        expPL.put(177.1022, 4321.1);
        expPL.put(181.0965, 2230.9);
        expPL.put(183.1481, 2386.5);
        expPL.put(187.1443, 2584.9);
        expPL.put(193.1771, 973.3);
        expPL.put(193.9456, 999.9);
        expPL.put(203.0669, 1522.8);
        expPL.put(204.1467, 4476.0);
        expPL.put(205.0984, 3468.5);
        expPL.put(211.1442, 2306.3);
        expPL.put(213.0866, 3319.8);
        expPL.put(215.1395, 3783.9);
        expPL.put(218.0538, 1071.4);
        expPL.put(222.1222, 3027.0);
        expPL.put(225.0987, 7068.6);
        expPL.put(227.2116, 1267.7);
        expPL.put(228.1833, 2891.7);
        expPL.put(231.0976, 1761.4);
        expPL.put(232.1418, 7508.7);
        expPL.put(236.0670, 2916.8);
        expPL.put(243.1342, 3272.2);
        expPL.put(249.1674, 1372.7);
        expPL.put(250.1182, 3844.7);
        expPL.put(253.0922, 10890.2);
        expPL.put(261.1573, 1203.1);
        expPL.put(272.1571, 1415.1);
        expPL.put(273.1354, 6343.8);
        expPL.put(274.1036, 11879.3);
        expPL.put(275.1829, 1893.8);
        expPL.put(291.2153, 30264.9);
        expPL.put(292.1297, 2932.0);
        expPL.put(292.2190, 1987.5);
        expPL.put(294.1812, 1632.1);
        expPL.put(300.1548, 1428.1);
        expPL.put(301.1293, 4879.2);
        expPL.put(302.1308, 1186.8);
        expPL.put(318.1796, 2664.3);
        expPL.put(326.1700, 2073.0);
        expPL.put(330.2252, 1181.4);
        expPL.put(330.6713, 1548.7);
        expPL.put(332.3052, 1129.4);
        expPL.put(335.2074, 1161.8);
        expPL.put(339.2780, 963.2);
        expPL.put(344.1798, 1387.7);
        expPL.put(347.1692, 23545.9);
        expPL.put(348.1698, 4840.3);
        expPL.put(348.2363, 12642.4);
        expPL.put(349.2343, 1524.4);
        expPL.put(351.1683, 1542.0);
        expPL.put(354.0852, 996.2);
        expPL.put(358.1861, 1961.3);
        expPL.put(359.1914, 1356.5);
        expPL.put(360.1999, 5332.0);
        expPL.put(363.2095, 1899.6);
        expPL.put(366.1758, 3871.1);
        expPL.put(367.2193, 973.4);
        expPL.put(372.1883, 4395.4);
        expPL.put(381.0766, 1240.8);
        expPL.put(381.1916, 3470.7);
        expPL.put(383.1915, 1258.6);
        expPL.put(385.1460, 1123.5);
        expPL.put(387.1874, 2391.7);
        expPL.put(390.2119, 3507.3);
        expPL.put(392.5189, 1123.4);
        expPL.put(398.2690, 950.3);
        expPL.put(401.1977, 1179.6);
        expPL.put(403.8663, 887.3);
        expPL.put(411.1767, 3625.1);
        expPL.put(416.1609, 4803.6);
        expPL.put(416.7136, 1821.0);
        expPL.put(417.1629, 1002.9);
        expPL.put(417.2076, 1630.3);
        expPL.put(418.2062, 34635.4);
        expPL.put(419.2094, 3802.9);
        expPL.put(420.2535, 1519.6);
        expPL.put(422.2556, 31139.8);
        expPL.put(423.2605, 5716.1);
        expPL.put(425.7220, 17201.0);
        expPL.put(426.2249, 5599.3);
        expPL.put(426.7368, 1185.5);
        expPL.put(428.7295, 1938.9);
        expPL.put(429.1893, 2288.2);
        expPL.put(437.2044, 2120.3);
        expPL.put(440.2235, 1124.0);
        expPL.put(446.2133, 1709.5);
        expPL.put(448.2876, 1258.4);
        expPL.put(459.2634, 1332.2);
        expPL.put(461.3202, 10869.4);
        expPL.put(462.3272, 1661.4);
        expPL.put(465.7163, 2317.7);
        expPL.put(466.2490, 2080.4);
        expPL.put(474.7217, 1176.9);
        expPL.put(475.0379, 1024.1);
        expPL.put(476.2911, 1745.6);
        expPL.put(483.2363, 8022.7);
        expPL.put(483.7374, 2946.8);
        expPL.put(489.2445, 12040.3);
        expPL.put(496.2640, 1260.0);
        expPL.put(503.2950, 9114.4);
        expPL.put(504.2952, 1353.2);
        expPL.put(504.7388, 2762.8);
        expPL.put(514.2886, 1079.4);
        expPL.put(530.7739, 1740.7);
        expPL.put(531.2903, 32317.4);
        expPL.put(532.2997, 17075.7);
        expPL.put(533.2879, 24674.9);
        expPL.put(534.2914, 4529.2);
        expPL.put(539.3586, 2183.2);
        expPL.put(539.7799, 4110.1);
        expPL.put(540.2771, 1713.2);
        expPL.put(542.2751, 1561.7);
        expPL.put(544.2159, 3660.0);
        expPL.put(546.2710, 1313.0);
        expPL.put(550.3146, 26210.3);
        expPL.put(551.3187, 4754.5);
        expPL.put(556.1799, 1238.8);
        expPL.put(569.2763, 1368.0);
        expPL.put(570.2671, 2539.6);
        expPL.put(572.2813, 8105.6);
        expPL.put(573.9478, 2083.2);
        expPL.put(574.4032, 6607.1);
        expPL.put(575.4091, 2930.4);
        expPL.put(590.2913, 8528.7);
        expPL.put(599.2716, 1478.6);
        expPL.put(613.3110, 2512.3);
        expPL.put(625.3989, 924.7);
        expPL.put(632.3397, 3203.7);
        expPL.put(641.8215, 4868.6);
        expPL.put(642.3206, 3806.7);
        expPL.put(652.0047, 1113.9);
        expPL.put(657.2900, 2223.3);
        expPL.put(659.3213, 1328.0);
        expPL.put(660.3326, 39168.7);
        expPL.put(661.3357, 9677.4);
        expPL.put(678.3456, 3133.3);
        expPL.put(679.3669, 1185.3);
        expPL.put(687.4901, 2620.8);
        expPL.put(691.3425, 3486.1);
        expPL.put(692.3361, 1324.4);
        expPL.put(696.3597, 2352.4);
        expPL.put(698.3649, 1959.4);
        expPL.put(700.3586, 1515.6);
        expPL.put(701.3195, 5317.1);
        expPL.put(702.3240, 1720.7);
        expPL.put(706.3347, 2077.0);
        expPL.put(708.3672, 3211.4);
        expPL.put(712.3572, 1694.6);
        expPL.put(713.3770, 55298.7);
        expPL.put(713.7130, 3837.0);
        expPL.put(714.3800, 16846.6);
        expPL.put(714.6431, 1419.5);
        expPL.put(715.3842, 3714.9);
        expPL.put(717.3498, 4312.5);
        expPL.put(719.3295, 19852.1);
        expPL.put(719.7042, 8663.0);
        expPL.put(720.0415, 4511.7);
        expPL.put(720.3317, 5783.2);
        expPL.put(721.0505, 1142.2);
        expPL.put(721.3361, 1564.1);
        expPL.put(722.3668, 3616.1);
        expPL.put(722.8681, 2830.2);
        expPL.put(725.3420, 1599.1);
        expPL.put(733.3091, 1896.6);
        expPL.put(742.3704, 1361.6);
        expPL.put(742.4747, 993.7);
        expPL.put(743.3707, 26163.1);
        expPL.put(744.3724, 8103.1);
        expPL.put(748.8880, 4101.9);
        expPL.put(749.3947, 2085.6);
        expPL.put(760.3939, 1335.6);
        expPL.put(761.3799, 32627.9);
        expPL.put(762.3816, 9047.3);
        expPL.put(763.3843, 1650.1);
        expPL.put(784.5378, 3597.6);
        expPL.put(785.5396, 2704.4);
        expPL.put(788.3420, 1248.3);
        expPL.put(790.3694, 8595.1);
        expPL.put(791.3788, 4219.6);
        expPL.put(804.3882, 2166.2);
        expPL.put(813.3947, 1498.9);
        expPL.put(813.9140, 1688.8);
        expPL.put(848.4129, 1650.9);
        expPL.put(850.4377, 23946.9);
        expPL.put(851.4390, 9381.5);
        expPL.put(852.4345, 1512.5);
        expPL.put(856.4562, 8774.4);
        expPL.put(857.4503, 3261.6);
        expPL.put(861.3882, 2563.0);
        expPL.put(873.4077, 9829.9);
        expPL.put(874.4785, 6073.5);
        expPL.put(875.4637, 2233.5);
        expPL.put(891.4241, 3879.8);
        expPL.put(892.4159, 3155.7);
        expPL.put(913.4674, 1960.2);
        expPL.put(930.5005, 1789.8);
        expPL.put(931.4854, 30842.3);
        expPL.put(932.4897, 8884.9);
        expPL.put(936.6047, 5614.8);
        expPL.put(937.6043, 4621.5);
        expPL.put(938.6149, 1440.1);
        expPL.put(947.4412, 2848.8);
        expPL.put(948.4778, 1532.9);
        expPL.put(964.4821, 1470.4);
        expPL.put(965.4653, 19734.0);
        expPL.put(966.4650, 8950.5);
        expPL.put(967.4569, 2107.8);
        expPL.put(993.4295, 1490.6);
        expPL.put(1006.4843, 1918.1);
        expPL.put(1029.8074, 1122.9);
        expPL.put(1037.4969, 2121.0);
        expPL.put(1050.5564, 2377.6);
        expPL.put(1060.5251, 2298.2);
        expPL.put(1078.5488, 16658.3);
        expPL.put(1079.5515, 6956.8);
        expPL.put(1162.9882, 1400.6);
        expPL.put(1191.6260, 2661.4);
        expPL.put(1225.6193, 3854.4);
        expPL.put(1226.6359, 1794.1);
        expPL.put(1264.6606, 1705.0);
        expPL.put(1265.6200, 1694.0);
        expPL.put(1282.6382, 12915.6);
        expPL.put(1283.6312, 7722.6);
        expPL.put(1284.6368, 3376.3);
        expPL.put(1306.6533, 7133.9);
        expPL.put(1307.6659, 5163.4);
        expPL.put(1339.3894, 1137.5);
        expPL.put(1395.7095, 2211.6);
        expPL.put(1396.7430, 1727.1);
        expPL.put(1462.6909, 1175.8);
        expPL.put(1496.7952, 1931.9);
        expPL.put(1621.1871, 1271.3);
        expPL.put(1645.7383, 1075.0);
        expPL.put(1791.5017, 1326.0);
        expPL.put(1797.4871, 1214.3);
    }

    @Test
    public void calIonFraction() {
        assertEquals((double) 17 / (double) 32 , Score.calIonFraction(ionMatrix1, 2, expPL, 0.01), 0.001);
        assertEquals((double) 14 / (double) 32 , Score.calIonFraction(ionMatrix2, 2, expPL, 0.01), 0.001);
    }

    @Test
    public void calMatchedHighestIntensityFraction() {
        assertEquals((double) 15 / (double) 17 , Score.calMatchedHighestIntensityFraction(ionMatrix1, 2, expPL, 0.01), 0.001);
        assertEquals((double) 12 / (double) 14 , Score.calMatchedHighestIntensityFraction(ionMatrix2, 2, expPL, 0.01), 0.001);
    }

    @Test
    public void calExplainedAAFraction() {
        assertEquals((double) 16 / (double) 16 , Score.calExplainedAAFraction(ionMatrix1, 2, expPL, 0.01), 0.001);
        assertEquals((double) 12 / (double) 16 , Score.calExplainedAAFraction(ionMatrix2, 2, expPL, 0.01), 0.001);
    }

    @Test
    public void calAScore() throws Exception {
        assertEquals(28.3382, Score.calAScore(new TreeMap<>(expPL), 6, binomial, varPTMMap1, ionMatrix1, varPTMMap2, ionMatrix2, 0.01, 18), 0.001);
        assertEquals(38.0304, Score.calAScore(new TreeMap<>(expPL), 6, binomial, varPTMMap1, ionMatrix1, null, null, 0.01, 18), 0.001);
    }

    @Test
    public void calAScoreSub() throws Exception {
        assertEquals(21.4021, Score.calAScoreSub(new TreeMap<>(expPL), 6, binomial, varPTMMap1, ionMatrix1, varPTMMap2, ionMatrix2, 0.01, 18), 0.001);
        assertEquals(31.2747, Score.calAScoreSub(new TreeMap<>(expPL), 6, binomial, varPTMMap1, ionMatrix1, null, null, 0.01, 18), 0.001);
    }

    @Test
    public void calBinomialScorePValue() throws Exception {
        assertEquals(1.34398E-17, Score.calBinomialScorePValue(new TreeMap<>(expPL), 6, binomial, 1, ionMatrix1, 0.01, 18), 0.001);
        assertEquals(6.37371E-13, Score.calBinomialScorePValue(new TreeMap<>(expPL), 6, binomial, 1, ionMatrix2, 0.01, 18), 0.001);
    }

    @Test
    public void calBinomialScorePValueSub () throws Exception {
        assertEquals(1.343981E-17, Score.calBinomialScorePValueSub(new TreeMap<>(expPL), 6, binomial, 1, ionMatrix1, 0.01, 18), 0.001);
        assertEquals(4.906808E-26, Score.calBinomialScorePValueSub(new TreeMap<>(expPL), 1, binomial, 1, ionMatrix1, 0.01, 18), 0.001);
        assertEquals(1.312552E-9, Score.calBinomialScorePValueSub(new TreeMap<>(expPL), 6, binomial, 1, ionMatrix2, 0.01, 18), 0.001);
        assertEquals(3.982398E-20, Score.calBinomialScorePValueSub(new TreeMap<>(expPL), 1, binomial, 1, ionMatrix2, 0.01, 18), 0.001);

    }

    @Test
    public void getAffectedIonSet() throws Exception {
        Object score = Class.forName("ProteomicsLibrary.Score").newInstance();
        Method method = score.getClass().getDeclaredMethod("getAffectedIonSet", TreeMap.class, int.class, Set.class, Set.class);
        method.setAccessible(true);

        Set<Integer> bIonResult = new TreeSet<>();
        Set<Integer> yIonResult = new TreeSet<>();

        method.invoke(score, varPTMMap1, 16, bIonResult, yIonResult);
        Integer[] bIonGroundTruth = new Integer[]{1, 10, 11, 15};
        Integer[] yIonGroundTruth = new Integer[]{1, 5, 6, 15};
        assertArrayEquals(bIonGroundTruth, bIonResult.toArray(new Integer[0]));
        assertArrayEquals(yIonGroundTruth, yIonResult.toArray(new Integer[0]));

        bIonResult.clear();
        yIonResult.clear();
        method.invoke(score, varPTMMap2, 16, bIonResult, yIonResult);
        bIonGroundTruth = new Integer[]{1, 13, 14, 15};
        yIonGroundTruth = new Integer[]{1, 2, 3, 15};
        assertArrayEquals(bIonGroundTruth, bIonResult.toArray(new Integer[0]));
        assertArrayEquals(yIonGroundTruth, yIonResult.toArray(new Integer[0]));
    }

    @Test
    public void getMatchedIonSet() throws Exception {
        Object score = Class.forName("ProteomicsLibrary.Score").newInstance();
        Method method = score.getClass().getDeclaredMethod("getMatchedIonSet", double[][].class, TreeMap.class, double.class, Set.class, Set.class);
        method.setAccessible(true);

        Set<Integer> affectedBIonSet = new TreeSet<>();
        Set<Integer> affectedYIonSet = new TreeSet<>();
        affectedBIonSet.add(1);
        affectedBIonSet.add(10);
        affectedBIonSet.add(11);
        affectedBIonSet.add(15);
        affectedYIonSet.add(1);
        affectedYIonSet.add(5);
        affectedYIonSet.add(6);
        affectedYIonSet.add(15);
        Set<String> temp = (HashSet<String>) method.invoke(score, ionMatrix1, new TreeMap<>(expPL), 0.01, affectedBIonSet, affectedYIonSet);
        String[] result = temp.toArray(new String[0]);
        Arrays.sort(result);
        String[] groundTruth = new String[]{"b1", "y1", "y5", "y6"};
        assertArrayEquals(groundTruth, result);

        affectedBIonSet.clear();
        affectedYIonSet.clear();
        affectedBIonSet.add(1);
        affectedBIonSet.add(13);
        affectedBIonSet.add(14);
        affectedBIonSet.add(15);
        affectedYIonSet.add(1);
        affectedYIonSet.add(2);
        affectedYIonSet.add(3);
        affectedYIonSet.add(15);
        temp = (HashSet<String>) method.invoke(score, ionMatrix2, new TreeMap<>(expPL), 0.01, affectedBIonSet, affectedYIonSet);
        result = temp.toArray(new String[0]);
        Arrays.sort(result);
        groundTruth = new String[]{"b1", "y1", "y2"};
        assertArrayEquals(groundTruth, result);
    }

    @Test
    public void getMatchedPeakNum() {
        assertEquals(17, Score.getMatchedIonNum(new TreeMap<>(expPL), 1, ionMatrix1, 0.01));
        assertEquals(14, Score.getMatchedIonNum(new TreeMap<>(expPL), 1, ionMatrix2, 0.01));
        assertEquals(22, Score.getMatchedIonNum(new TreeMap<>(expPL), 1, ionMatrix1, 0.5));
        assertEquals(21, Score.getMatchedIonNum(new TreeMap<>(expPL), 1, ionMatrix2, 0.5));
        assertEquals(30, Score.getMatchedIonNum(new TreeMap<>(expPL), 2, ionMatrix1, 0.01));
        assertEquals(25, Score.getMatchedIonNum(new TreeMap<>(expPL), 2, ionMatrix2, 0.01));
    }
}