package ProteomicsLibrary;

import ProteomicsLibrary.Types.AA;
import ProteomicsLibrary.Types.SparseVector;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MassTool {

    private static final Pattern modAAPattern = Pattern.compile("([A-Znc])([(\\[]([0-9.\\-]+)[)\\]])?");
    private static final Pattern leftFlankPattern = Pattern.compile("^[A-Z-]\\.");
    private static final Pattern rightFlankPattern = Pattern.compile("\\.[A-Z-]$");

    public static final double PROTON = 1.00727646688;
    public static final double C13_DIFF = 1.00335483;
    public static final double N1514_DIFF = 15.0001088 - 14.0030732;
    public static final Map<Character, Double> n1514DiffMap = new HashMap<>(); // This map only contains amino acids. The PTM different are considered in the identification (PIPI).
    static {
        n1514DiffMap.put('G', N1514_DIFF);
        n1514DiffMap.put('A', N1514_DIFF);
        n1514DiffMap.put('S', N1514_DIFF);
        n1514DiffMap.put('P', N1514_DIFF);
        n1514DiffMap.put('V', N1514_DIFF);
        n1514DiffMap.put('T', N1514_DIFF);
        n1514DiffMap.put('C', N1514_DIFF);
        n1514DiffMap.put('I', N1514_DIFF);
        n1514DiffMap.put('L', N1514_DIFF);
        n1514DiffMap.put('N', N1514_DIFF * 2);
        n1514DiffMap.put('D', N1514_DIFF);
        n1514DiffMap.put('Q', N1514_DIFF * 2);
        n1514DiffMap.put('K', N1514_DIFF * 2);
        n1514DiffMap.put('E', N1514_DIFF);
        n1514DiffMap.put('M', N1514_DIFF);
        n1514DiffMap.put('H', N1514_DIFF * 3);
        n1514DiffMap.put('F', N1514_DIFF);
        n1514DiffMap.put('R', N1514_DIFF * 4);
        n1514DiffMap.put('Y', N1514_DIFF);
        n1514DiffMap.put('W', N1514_DIFF * 2);
        n1514DiffMap.put('U', N1514_DIFF);
        n1514DiffMap.put('O', N1514_DIFF * 3);
        n1514DiffMap.put('n', 0d);
        n1514DiffMap.put('c', 0d);
    }

    public final double H2O;

    private final Map<String, Double> elementTable = new HashMap<>();
    private final Map<Character, Double> massTable = new HashMap<>(30, 1);
    private final int missedCleavage;
    private final double ms2Tolerance;
    private final double inverse2Ms2Tolerance;
    private final double oneMinusBinOffset;
    private final Pattern digestSitePattern;
    private final boolean cleavageFromCTerm;
    private final Pattern digestSitePatternForLinkSiteChecking; // this is for removing the digest site form cross-linking
    private final String labelling;
    private final Map<Character, Double> fixModMap;

    public MassTool(int missedCleavage, Map<Character, Double> fixModMap, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, double ms2Tolerance, double oneMinusBinOffset, String labelling) {
        this.labelling = labelling;
        this.fixModMap = fixModMap;

        elementTable.put("-", 0d);
        elementTable.put("H", 1.0078246);
        elementTable.put("He", 3.01603);
        elementTable.put("Li", 6.015121);
        elementTable.put("Be", 9.012182);
        elementTable.put("B", 10.012937);
        elementTable.put("C", 12.0000000);
        elementTable.put("N", 14.0030732);
        if (labelling.contentEquals("N15")) {
            elementTable.put("N", 15.0001088);
        }
        elementTable.put("O", 15.9949141);
        elementTable.put("F", 18.9984032);
        elementTable.put("Ne", 19.992435);
        elementTable.put("Na", 22.989767);
        elementTable.put("Mg", 23.985042);
        elementTable.put("Al", 26.981539);
        elementTable.put("Si", 27.976927);
        elementTable.put("P", 30.973762);
        elementTable.put("S", 31.972070);
        elementTable.put("Cl", 34.9688531);
        elementTable.put("Ar", 35.967545);
        elementTable.put("K", 38.963707);
        elementTable.put("Ca", 39.962591);
        elementTable.put("Sc", 44.955910);
        elementTable.put("Ti", 45.952629);
        elementTable.put("V", 49.947161);
        elementTable.put("Cr", 49.946046);
        elementTable.put("Mn", 54.938047);
        elementTable.put("Fe", 53.939612);
        elementTable.put("Co", 58.933198);
        elementTable.put("Ni", 57.935346);
        elementTable.put("Cu", 62.939598);
        elementTable.put("Zn", 63.929145);
        elementTable.put("Ga", 68.925580);
        elementTable.put("Ge", 69.924250);
        elementTable.put("As", 74.921594);
        elementTable.put("Se", 73.922475);
        elementTable.put("Br", 78.918336);
        elementTable.put("Kr", 77.914);
        elementTable.put("Rb", 84.911794);
        elementTable.put("Sr", 83.913430);
        elementTable.put("Y", 88.905849);
        elementTable.put("Zr", 89.904703);
        elementTable.put("Nb", 92.906377);
        elementTable.put("Mo", 91.906808);
        elementTable.put("Tc", 98.0);
        elementTable.put("Ru", 95.907599);
        elementTable.put("Rh", 102.905500);
        elementTable.put("Pd", 101.905634);
        elementTable.put("Ag", 106.905092);
        elementTable.put("Cd", 105.906461);
        elementTable.put("In", 112.904061);
        elementTable.put("Sn", 111.904826);
        elementTable.put("Sb", 120.903821);
        elementTable.put("Te", 119.904048);
        elementTable.put("I", 126.904473);
        elementTable.put("Xe", 123.905894);
        elementTable.put("Cs", 132.905429);
        elementTable.put("Ba", 129.906282);
        elementTable.put("La", 137.90711);
        elementTable.put("Ce", 135.907140);
        elementTable.put("Pr", 140.907647);
        elementTable.put("Nd", 141.907719);
        elementTable.put("Pm", 145.0);
        elementTable.put("Sm", 143.911998);
        elementTable.put("Eu", 150.919847);
        elementTable.put("Gd", 151.919786);
        elementTable.put("Tb", 158.925342);
        elementTable.put("Dy", 155.925277);
        elementTable.put("Ho", 164.930319);
        elementTable.put("Er", 161.928775);
        elementTable.put("Tm", 168.934212);
        elementTable.put("Yb", 167.933894);
        elementTable.put("Lu", 174.940770);
        elementTable.put("Hf", 173.940044);
        elementTable.put("Ta", 179.947462);
        elementTable.put("W", 179.946701);
        elementTable.put("Re", 184.952951);
        elementTable.put("Os", 183.952488);
        elementTable.put("Ir", 190.960584);
        elementTable.put("Pt", 189.959917);
        elementTable.put("Au", 196.966543);
        elementTable.put("Hg", 195.965807);
        elementTable.put("Tl", 202.972320);
        elementTable.put("Pb", 203.973020);
        elementTable.put("Bi", 208.980374);
        elementTable.put("Po", 209.0);
        elementTable.put("At", 210.0);
        elementTable.put("Rn", 222.0);
        elementTable.put("Fr", 223.0);
        elementTable.put("Ra", 226.025);
        // elementTable.put("Ac", 227.028); // conflict with Unimod bricks
        elementTable.put("Th", 232.038054);
        elementTable.put("Pa", 231.0359);
        elementTable.put("U", 234.040946);
        elementTable.put("Np", 237.048);
        elementTable.put("Pu", 244.0);
        elementTable.put("Am", 243.0);
        elementTable.put("Cm", 247.0);
        elementTable.put("Bk", 247.0);
        elementTable.put("Cf", 251.0);
        elementTable.put("Es", 252.0);
        elementTable.put("Fm", 257.0);
        elementTable.put("Md", 258.0);
        elementTable.put("No", 259.0);
        elementTable.put("Lr", 260.0);
        elementTable.put("13C", 13.0033554);
        elementTable.put("15N", 15.0001088);
        elementTable.put("18O", 17.9991616);
        elementTable.put("2H", 2.0141021);
        elementTable.put("dHex", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 10);
        elementTable.put("Hep", elementTable.get("C") * 7 + elementTable.get("O") * 6 + elementTable.get("H") * 12);
        elementTable.put("Hex", elementTable.get("C") * 6 + elementTable.get("O") * 5 + elementTable.get("H") * 10);
        elementTable.put("HexA", elementTable.get("C") * 6 + elementTable.get("O") * 6 + elementTable.get("H") * 8);
        elementTable.put("HexN", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 11 + elementTable.get("N"));
        elementTable.put("HexNAc", elementTable.get("C") * 8 + elementTable.get("O") * 5 + + elementTable.get("N") + elementTable.get("H") * 13);
        elementTable.put("Kdn", elementTable.get("C") * 9 + elementTable.get("H") * 14 + elementTable.get("O") * 8);
        elementTable.put("Kdo", elementTable.get("C") * 8 + elementTable.get("H") * 12 + elementTable.get("O") * 7);
        elementTable.put("NeuAc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 8 + elementTable.get("N"));
        elementTable.put("NeuGc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 9 + elementTable.get("N"));
        elementTable.put("Pent", elementTable.get("C") * 5 + elementTable.get("O") * 4 + elementTable.get("H") * 8);
        elementTable.put("Phos", elementTable.get("O") * 3 + elementTable.get("H") + elementTable.get("P"));
        elementTable.put("Sulf", elementTable.get("S") + elementTable.get("O") * 3);
        elementTable.put("Water", elementTable.get("H") * 2 + elementTable.get("O"));
        elementTable.put("Me", elementTable.get("C") + elementTable.get("H") * 2);
        elementTable.put("Ac", elementTable.get("C") * 2 + elementTable.get("H") * 2 + elementTable.get("O")); // Caution! This is not Actinium

        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        inverse2Ms2Tolerance = 1 / (2 * ms2Tolerance);
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm = cleavageFromCTerm;
        massTable.put('G', (elementTable.get("C") * 2 + elementTable.get("H") * 3 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('G')));
        massTable.put('A', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('A')));
        massTable.put('S', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('S')));
        massTable.put('P', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('P')));
        massTable.put('V', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('V')));
        massTable.put('T', (elementTable.get("C") * 4 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('I')));
        massTable.put('C', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('C')));
        massTable.put('I', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('I')));
        massTable.put('L', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('L')));
        massTable.put('N', (elementTable.get("C") * 4 + elementTable.get("H") * 6 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('N')));
        massTable.put('D', (elementTable.get("C") * 4 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('D')));
        massTable.put('Q', (elementTable.get("C") * 5 + elementTable.get("H") * 8 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('Q')));
        massTable.put('K', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('K')));
        massTable.put('E', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('E')));
        massTable.put('M', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('M')));
        massTable.put('H', (elementTable.get("C") * 6 + elementTable.get("H") * 7 + elementTable.get("N") * 3 + elementTable.get("O") + fixModMap.get('H')));
        massTable.put('F', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('F')));
        massTable.put('R', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 4 + elementTable.get("O") + fixModMap.get('R')));
        massTable.put('Y', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('Y')));
        massTable.put('W', (elementTable.get("C") * 11 + elementTable.get("H") * 10 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('W')));
        massTable.put('U', (elementTable.get("C") * 3 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + elementTable.get("Se") + fixModMap.get('U')));
        massTable.put('O', (elementTable.get("C") * 12 + elementTable.get("H") * 21 + elementTable.get("N") * 3 + elementTable.get("O") * 3 + fixModMap.get('O')));
        massTable.put('n', fixModMap.get('n'));
        massTable.put('c', fixModMap.get('c'));
        massTable.put('#', massTable.get('I')); // for I and L.
        massTable.put('$', (massTable.get('Q') + massTable.get('K')) * 0.5); // for Q and K.
        H2O = elementTable.get("H") * 2 + elementTable.get("O");

        digestSitePattern = getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
        digestSitePatternForLinkSiteChecking = getDigestSitePatternForLinkSiteChecking(cleavageSite, protectionSite, cleavageFromCTerm);
    }

    public MassTool(int missedCleavage, String cleavageSite, String protectionSite, boolean cleavageFromCTerm, double ms2Tolerance, double oneMinusBinOffset, String labelling) {
        this.labelling = labelling;

        // build a default fixModMap;
        this.fixModMap = new HashMap<>();
        fixModMap.put('G', 0d);
        fixModMap.put('A', 0d);
        fixModMap.put('S', 0d);
        fixModMap.put('P', 0d);
        fixModMap.put('V', 0d);
        fixModMap.put('T', 0d);
        fixModMap.put('C', 57.02146);
        fixModMap.put('I', 0d);
        fixModMap.put('L', 0d);
        fixModMap.put('N', 0d);
        fixModMap.put('D', 0d);
        fixModMap.put('Q', 0d);
        fixModMap.put('K', 0d);
        fixModMap.put('E', 0d);
        fixModMap.put('M', 0d);
        fixModMap.put('H', 0d);
        fixModMap.put('F', 0d);
        fixModMap.put('R', 0d);
        fixModMap.put('Y', 0d);
        fixModMap.put('W', 0d);
        fixModMap.put('U', 0d);
        fixModMap.put('O', 0d);
        fixModMap.put('n', 0d);
        fixModMap.put('c', 0d);

        elementTable.put("-", 0d);
        elementTable.put("H", 1.0078246);
        elementTable.put("He", 3.01603);
        elementTable.put("Li", 6.015121);
        elementTable.put("Be", 9.012182);
        elementTable.put("B", 10.012937);
        elementTable.put("C", 12.0000000);
        elementTable.put("N", 14.0030732);
        if (labelling.contentEquals("N15")) {
            elementTable.put("N", 15.0001088);
        }
        elementTable.put("O", 15.9949141);
        elementTable.put("F", 18.9984032);
        elementTable.put("Ne", 19.992435);
        elementTable.put("Na", 22.989767);
        elementTable.put("Mg", 23.985042);
        elementTable.put("Al", 26.981539);
        elementTable.put("Si", 27.976927);
        elementTable.put("P", 30.973762);
        elementTable.put("S", 31.972070);
        elementTable.put("Cl", 34.9688531);
        elementTable.put("Ar", 35.967545);
        elementTable.put("K", 38.963707);
        elementTable.put("Ca", 39.962591);
        elementTable.put("Sc", 44.955910);
        elementTable.put("Ti", 45.952629);
        elementTable.put("V", 49.947161);
        elementTable.put("Cr", 49.946046);
        elementTable.put("Mn", 54.938047);
        elementTable.put("Fe", 53.939612);
        elementTable.put("Co", 58.933198);
        elementTable.put("Ni", 57.935346);
        elementTable.put("Cu", 62.939598);
        elementTable.put("Zn", 63.929145);
        elementTable.put("Ga", 68.925580);
        elementTable.put("Ge", 69.924250);
        elementTable.put("As", 74.921594);
        elementTable.put("Se", 73.922475);
        elementTable.put("Br", 78.918336);
        elementTable.put("Kr", 77.914);
        elementTable.put("Rb", 84.911794);
        elementTable.put("Sr", 83.913430);
        elementTable.put("Y", 88.905849);
        elementTable.put("Zr", 89.904703);
        elementTable.put("Nb", 92.906377);
        elementTable.put("Mo", 91.906808);
        elementTable.put("Tc", 98.0);
        elementTable.put("Ru", 95.907599);
        elementTable.put("Rh", 102.905500);
        elementTable.put("Pd", 101.905634);
        elementTable.put("Ag", 106.905092);
        elementTable.put("Cd", 105.906461);
        elementTable.put("In", 112.904061);
        elementTable.put("Sn", 111.904826);
        elementTable.put("Sb", 120.903821);
        elementTable.put("Te", 119.904048);
        elementTable.put("I", 126.904473);
        elementTable.put("Xe", 123.905894);
        elementTable.put("Cs", 132.905429);
        elementTable.put("Ba", 129.906282);
        elementTable.put("La", 137.90711);
        elementTable.put("Ce", 135.907140);
        elementTable.put("Pr", 140.907647);
        elementTable.put("Nd", 141.907719);
        elementTable.put("Pm", 145.0);
        elementTable.put("Sm", 143.911998);
        elementTable.put("Eu", 150.919847);
        elementTable.put("Gd", 151.919786);
        elementTable.put("Tb", 158.925342);
        elementTable.put("Dy", 155.925277);
        elementTable.put("Ho", 164.930319);
        elementTable.put("Er", 161.928775);
        elementTable.put("Tm", 168.934212);
        elementTable.put("Yb", 167.933894);
        elementTable.put("Lu", 174.940770);
        elementTable.put("Hf", 173.940044);
        elementTable.put("Ta", 179.947462);
        elementTable.put("W", 179.946701);
        elementTable.put("Re", 184.952951);
        elementTable.put("Os", 183.952488);
        elementTable.put("Ir", 190.960584);
        elementTable.put("Pt", 189.959917);
        elementTable.put("Au", 196.966543);
        elementTable.put("Hg", 195.965807);
        elementTable.put("Tl", 202.972320);
        elementTable.put("Pb", 203.973020);
        elementTable.put("Bi", 208.980374);
        elementTable.put("Po", 209.0);
        elementTable.put("At", 210.0);
        elementTable.put("Rn", 222.0);
        elementTable.put("Fr", 223.0);
        elementTable.put("Ra", 226.025);
        // elementTable.put("Ac", 227.028); // conflict with Unimod bricks
        elementTable.put("Th", 232.038054);
        elementTable.put("Pa", 231.0359);
        elementTable.put("U", 234.040946);
        elementTable.put("Np", 237.048);
        elementTable.put("Pu", 244.0);
        elementTable.put("Am", 243.0);
        elementTable.put("Cm", 247.0);
        elementTable.put("Bk", 247.0);
        elementTable.put("Cf", 251.0);
        elementTable.put("Es", 252.0);
        elementTable.put("Fm", 257.0);
        elementTable.put("Md", 258.0);
        elementTable.put("No", 259.0);
        elementTable.put("Lr", 260.0);
        elementTable.put("13C", 13.0033554);
        elementTable.put("15N", 15.0001088);
        elementTable.put("18O", 17.9991616);
        elementTable.put("2H", 2.0141021);
        elementTable.put("dHex", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 10);
        elementTable.put("Hep", elementTable.get("C") * 7 + elementTable.get("O") * 6 + elementTable.get("H") * 12);
        elementTable.put("Hex", elementTable.get("C") * 6 + elementTable.get("O") * 5 + elementTable.get("H") * 10);
        elementTable.put("HexA", elementTable.get("C") * 6 + elementTable.get("O") * 6 + elementTable.get("H") * 8);
        elementTable.put("HexN", elementTable.get("C") * 6 + elementTable.get("O") * 4 + elementTable.get("H") * 11 + elementTable.get("N"));
        elementTable.put("HexNAc", elementTable.get("C") * 8 + elementTable.get("O") * 5 + + elementTable.get("N") + elementTable.get("H") * 13);
        elementTable.put("Kdn", elementTable.get("C") * 9 + elementTable.get("H") * 14 + elementTable.get("O") * 8);
        elementTable.put("Kdo", elementTable.get("C") * 8 + elementTable.get("H") * 12 + elementTable.get("O") * 7);
        elementTable.put("NeuAc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 8 + elementTable.get("N"));
        elementTable.put("NeuGc", elementTable.get("C") * 11 + elementTable.get("H") * 17 + elementTable.get("O") * 9 + elementTable.get("N"));
        elementTable.put("Pent", elementTable.get("C") * 5 + elementTable.get("O") * 4 + elementTable.get("H") * 8);
        elementTable.put("Phos", elementTable.get("O") * 3 + elementTable.get("H") + elementTable.get("P"));
        elementTable.put("Sulf", elementTable.get("S") + elementTable.get("O") * 3);
        elementTable.put("Water", elementTable.get("H") * 2 + elementTable.get("O"));
        elementTable.put("Me", elementTable.get("C") + elementTable.get("H") * 2);
        elementTable.put("Ac", elementTable.get("C") * 2 + elementTable.get("H") * 2 + elementTable.get("O")); // Caution! This is not Actinium

        this.missedCleavage = missedCleavage;
        this.ms2Tolerance = ms2Tolerance;
        inverse2Ms2Tolerance = 1 / (2 * ms2Tolerance);
        this.oneMinusBinOffset = oneMinusBinOffset;
        this.cleavageFromCTerm = cleavageFromCTerm;
        massTable.put('G', (elementTable.get("C") * 2 + elementTable.get("H") * 3 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('G')));
        massTable.put('A', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('A')));
        massTable.put('S', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('S')));
        massTable.put('P', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('P')));
        massTable.put('V', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('V')));
        massTable.put('T', (elementTable.get("C") * 4 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('I')));
        massTable.put('C', (elementTable.get("C") * 3 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('C')));
        massTable.put('I', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('I')));
        massTable.put('L', (elementTable.get("C") * 6 + elementTable.get("H") * 11 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('L')));
        massTable.put('N', (elementTable.get("C") * 4 + elementTable.get("H") * 6 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('N')));
        massTable.put('D', (elementTable.get("C") * 4 + elementTable.get("H") * 5 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('D')));
        massTable.put('Q', (elementTable.get("C") * 5 + elementTable.get("H") * 8 + elementTable.get("N") * 2 + elementTable.get("O") * 2 + fixModMap.get('Q')));
        massTable.put('K', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('K')));
        massTable.put('E', (elementTable.get("C") * 5 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 3 + fixModMap.get('E')));
        massTable.put('M', (elementTable.get("C") * 5 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + elementTable.get("S") + fixModMap.get('M')));
        massTable.put('H', (elementTable.get("C") * 6 + elementTable.get("H") * 7 + elementTable.get("N") * 3 + elementTable.get("O") + fixModMap.get('H')));
        massTable.put('F', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") + fixModMap.get('F')));
        massTable.put('R', (elementTable.get("C") * 6 + elementTable.get("H") * 12 + elementTable.get("N") * 4 + elementTable.get("O") + fixModMap.get('R')));
        massTable.put('Y', (elementTable.get("C") * 9 + elementTable.get("H") * 9 + elementTable.get("N") + elementTable.get("O") * 2 + fixModMap.get('Y')));
        massTable.put('W', (elementTable.get("C") * 11 + elementTable.get("H") * 10 + elementTable.get("N") * 2 + elementTable.get("O") + fixModMap.get('W')));
        massTable.put('U', (elementTable.get("C") * 3 + elementTable.get("H") * 7 + elementTable.get("N") + elementTable.get("O") * 2 + elementTable.get("Se") + fixModMap.get('U')));
        massTable.put('O', (elementTable.get("C") * 12 + elementTable.get("H") * 21 + elementTable.get("N") * 3 + elementTable.get("O") * 3 + fixModMap.get('O')));
        massTable.put('n', fixModMap.get('n'));
        massTable.put('c', fixModMap.get('c'));
        massTable.put('#', massTable.get('I')); // for I and L.
        massTable.put('$', (massTable.get('Q') + massTable.get('K')) * 0.5); // for Q and K.
        H2O = elementTable.get("H") * 2 + elementTable.get("O");

        digestSitePattern = getDigestSitePattern(cleavageSite, protectionSite, cleavageFromCTerm);
        digestSitePatternForLinkSiteChecking = getDigestSitePatternForLinkSiteChecking(cleavageSite, protectionSite, cleavageFromCTerm);
    }

    public static boolean isAA(char aa) {
        return aa == 'H' || aa == 'I' || aa == 'L' || aa == 'K' || aa == 'M' || aa == 'F' || aa == 'T' || aa == 'W' || aa == 'V' || aa == 'R' || aa == 'C' || aa == 'Q' || aa == 'G' || aa == 'P' || aa == 'Y' || aa == 'A' || aa == 'D' || aa == 'N' || aa == 'E' || aa == 'S' || aa == 'U' || aa == 'O';
    }

    public static boolean containsNonAAAndNC(String sequence) {
        for (char aa : sequence.toCharArray()) {
            if (!(aa == 'H' || aa == 'I' || aa == 'L' || aa == 'K' || aa == 'M' || aa == 'F' || aa == 'T' || aa == 'W' || aa == 'V' || aa == 'R' || aa == 'C' || aa == 'Q' || aa == 'G' || aa == 'P' || aa == 'Y' || aa == 'A' || aa == 'D' || aa == 'N' || aa == 'E' || aa == 'S' || aa == 'U' || aa == 'O' || aa == 'n' || aa == 'c')) {
                return true;
            }
        }
        return false;
    }

    public double calResidueMass(String sequence) { // n and c are also AA. Consider fixed modification automatically
        double totalMass = 0;
        Matcher matcher = modAAPattern.matcher(sequence);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double deltaMass = 0;
            if (matcher.group(3) != null) {
                deltaMass = Double.valueOf(matcher.group(3));
            }
            totalMass += massTable.get(aa) + deltaMass;
        }

        return totalMass;
    }

    public double calResidueMass2(String sequence) { // n and c are also AA. Don't consider fixed modification automatically
        double totalMass = 0;
        Matcher matcher = modAAPattern.matcher(sequence);
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double deltaMass = 0;
            if (matcher.group(3) != null) {
                deltaMass = Double.valueOf(matcher.group(3));
            }
            totalMass += massTable.get(aa) - fixModMap.get(aa) + deltaMass;
        }

        return totalMass;
    }

    public static AA[] seqToAAList(String sequence) { // n and c are also AA.
        Matcher matcher = modAAPattern.matcher(sequence);
        List<AA> temp = new LinkedList<>();
        while (matcher.find()) {
            char aa = matcher.group(1).charAt(0);
            double deltaMass = 0;
            if (matcher.group(3) != null) {
                deltaMass = Double.valueOf(matcher.group(3));
            }
            temp.add(new AA(aa, deltaMass));
        }
        return temp.toArray(new AA[0]);
    }

    public Set<String> buildPeptideSet(String proteinSequence) {
        Map<Integer, List<int[]>> digestRangeMap = digest(proteinSequence);
        Set<String> peptideSeqSet = new HashSet<>();

        for (int i = 0; i <= missedCleavage; ++i) {
            for (int[] digestRange : digestRangeMap.get(i)) {
                String subString = proteinSequence.substring(digestRange[0], digestRange[1]);
                peptideSeqSet.add("n" + subString + "c");
            }
        }

        // consider first "M" situation
        if (proteinSequence.startsWith("M")) {
            String newSequence = proteinSequence.substring(1);
            digestRangeMap = digest(newSequence);
            for (int i = 0; i <= missedCleavage; ++i) {
                if (!digestRangeMap.get(i).isEmpty()) {
                    int[] digestRange = digestRangeMap.get(i).get(0);
                    String subString = newSequence.substring(digestRange[0], digestRange[1]);
                    peptideSeqSet.add("n" + subString + "c");
                }
            }
        }
        return peptideSeqSet;
    }

    public double[][] buildIonArray(String sequence, int maxCharge) { // there are n and c in the sequence
        AA[] aaArray = seqToAAList(sequence);

        double[] inverseChargeArray = new double[maxCharge];
        for (int charge = 1; charge <= maxCharge; ++charge) {
            inverseChargeArray[charge - 1] = (double) 1 / (double) charge;
        }

        double[][] peptideIonArray = new double[2 * maxCharge][aaArray.length - 2];
        // traverse the sequence to get b-ion
        double bIonMass = massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass; // add N-term modification
        for (int i = 1; i < aaArray.length - 2; ++i) {
            bIonMass += massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            for (int charge = 1; charge <= maxCharge; ++charge) {
                peptideIonArray[2 * (charge - 1)][i - 1]  = bIonMass * inverseChargeArray[charge - 1] + PROTON;
            }
        }
        // calculate the last b-ion with C-term modification
        bIonMass +=  massTable.get(aaArray[aaArray.length - 2].aa) + aaArray[aaArray.length - 2].ptmDeltaMass + massTable.get(aaArray[aaArray.length - 1].aa) + aaArray[aaArray.length - 1].ptmDeltaMass;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1)][aaArray.length - 3] = bIonMass * inverseChargeArray[charge - 1] + PROTON;
        }

        // traverse the sequence with reversed order to get y-ion
        // the whole sequence
        double yIonMass = bIonMass + H2O;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1) + 1][0] = yIonMass * inverseChargeArray[charge - 1] + PROTON;
        }
        // delete the first amino acid and N-term modification
        yIonMass -= massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass + massTable.get(aaArray[1].aa) + aaArray[1].ptmDeltaMass;
        for (int charge = 1; charge <= maxCharge; ++charge) {
            peptideIonArray[2 * (charge - 1) + 1][1] = yIonMass * inverseChargeArray[charge - 1] + PROTON;
        }

        // rest of the sequence
        for (int i = 2; i < aaArray.length - 2; ++i) {
            yIonMass -= massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            for (int charge = 1; charge <= maxCharge; ++charge) {
                peptideIonArray[2 * (charge - 1) + 1][i] = yIonMass * inverseChargeArray[charge - 1] + PROTON;
            }
        }

        return peptideIonArray;
    }

    public double buildVectorAndCalXCorr(double[][] ionMatrix, int precursorCharge, SparseVector xcorrPL) {
        int colNum = ionMatrix[0].length;
        int rowNum = Math.min(ionMatrix.length / 2, precursorCharge - 1) * 2;
        if (precursorCharge == 1) {
            rowNum = 2;
        }

        double xcorr = 0;
        for (int i = 0; i < rowNum; ++i) {
            for (int j = 0; j < colNum; ++j) {
                xcorr += xcorrPL.get(mzToBin(ionMatrix[i][j]));
            }
        }

        return xcorr * 0.25;
    }

    public Map<Character, Double> getMassTable() {
        return massTable;
    }

    public Map<String, Double> getElementTable() {
        return elementTable;
    }

    public int mzToBin(double mz) {
        return (int) Math.floor(mz * inverse2Ms2Tolerance + oneMinusBinOffset);
    }

    public double binToMz(int idx) {
        return (idx - oneMinusBinOffset) * 2 * ms2Tolerance;
    }

    public String getLabelling() {
        return labelling;
    }

    public Pattern getDigestSitePattern() {
        return digestSitePattern;
    }

    public static String aaListToSeq(AA[] aaArray) {
        StringBuilder sb = new StringBuilder();
        for (AA aa : aaArray) {
            sb.append(aa.toString());
        }
        return sb.toString();
    }

    public static String unifyPeptide(String peptide) {
        AA[] aaArray = seqToAAList(deleteLeftRightFlankingAddNC(peptide));
        StringBuilder sb = new StringBuilder();
        for (AA aa : aaArray) {
            sb.append(aa.toString());
        }
        return sb.toString();
    }

    public static String L2I(String peptide) {
        return peptide.replaceAll("L", "I");
    }

    static Pattern getDigestSitePattern(String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
        Pattern digestSitePattern;
        if (cleavageFromCTerm) {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
            } else {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "](?![" + protectionSite + "])");
            }
        } else {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern = Pattern.compile("[" + cleavageSite + "]");
            } else {
                digestSitePattern = Pattern.compile("(?<![" + protectionSite + "])" + "[" + cleavageSite + "]");
            }
        }
        return digestSitePattern;
    }

    private static Pattern getDigestSitePatternForLinkSiteChecking(String cleavageSite, String protectionSite, boolean cleavageFromCTerm) {
        Pattern digestSitePattern2;
        if (cleavageFromCTerm) {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern2 = Pattern.compile("[" + cleavageSite + "]$");
            } else {
                digestSitePattern2 = Pattern.compile("[" + cleavageSite + "](?![" + protectionSite + "])$");
            }
        } else {
            if (protectionSite.contentEquals("-")) {
                digestSitePattern2 = Pattern.compile("^[" + cleavageSite + "]");
            } else {
                digestSitePattern2 = Pattern.compile("^(?<![" + protectionSite + "])" + "[" + cleavageSite + "]");
            }
        }
        return digestSitePattern2;
    }

    // Cross-linking part
    public Set<String> buildChainSet(String proteinSequence, short linkerType) {
        Map<Integer, List<int[]>> digestRangeMap = digest(proteinSequence);
        Set<String> chainSequenceSet = new HashSet<>();

        for (int i = 0; i <= missedCleavage; ++i) {
            for (int[] digestRange : digestRangeMap.get(i)) {
                String subString = proteinSequence.substring(digestRange[0], digestRange[1]);
                Matcher tempMatcher = digestSitePatternForLinkSiteChecking1.matcher(subString);
                String tempString = tempMatcher.replaceAll("");
                if (linkerType == 1 && tempString.contains("K")) {
                    chainSequenceSet.add("n" + subString + "c");
                } else if (linkerType == 2 && tempString.contains("C")) {
                    chainSequenceSet.add("n" + subString + "c");
                }

                if (digestRange[1] == proteinSequence.length()) {
                    // This is the end of the protein. No digestion site, so the link-sites in any position including C-term can be linked.
                    if (linkerType == 1 && subString.contains("K")) {
                        chainSequenceSet.add("n" + subString + "c");
                    } else if (linkerType == 2 && subString.contains("C")) {
                        chainSequenceSet.add("n" + subString + "c");
                    }
                }
            }
            if (linkerType == 1) {
                // Add N-term peptide
                if (!digestRangeMap.get(i).isEmpty()) {
                    int[] digestRange = digestRangeMap.get(i).get(0);
                    String subString = proteinSequence.substring(digestRange[0], digestRange[1]);
                    chainSequenceSet.add("n" + subString + "c");
                }
            }
        }

        if (proteinSequence.startsWith("M")) {
            String newSequence = proteinSequence.substring(1);
            digestRangeMap = digest(newSequence);

            for (int i = 0; i <= missedCleavage; ++i) {
                if (!digestRangeMap.get(i).isEmpty()) {
                    int[] digestRange = digestRangeMap.get(i).get(0);
                    String subString = newSequence.substring(digestRange[0], digestRange[1]);
                    Matcher tempMatcher = digestSitePatternForLinkSiteChecking1.matcher(subString);
                    String tempString = tempMatcher.replaceAll("");
                    if (linkerType == 1 && tempString.contains("K")) {
                        chainSequenceSet.add("n" + subString + "c");
                    } else if (linkerType == 2 && tempString.contains("C")) {
                        chainSequenceSet.add("n" + subString + "c");
                    }

                    if (digestRange[1] == newSequence.length()) {
                        // This is the end of the protein. No digestion site, so the link-sites in any position including C-term can be linked.
                        if (linkerType == 1 && subString.contains("K")) {
                            chainSequenceSet.add("n" + subString + "c");
                        } else if (linkerType == 2 && subString.contains("C")) {
                            chainSequenceSet.add("n" + subString + "c");
                        }
                    }

                    if (linkerType == 1) {
                        // Add N-term peptide
                        if (!digestRangeMap.get(i).isEmpty()) {
                            digestRange = digestRangeMap.get(i).get(0);
                            subString = newSequence.substring(digestRange[0], digestRange[1]);
                            chainSequenceSet.add("n" + subString + "c");
                        }
                    }
                }
            }
        }
        return chainSequenceSet;
    }

    public double generateTheoFragmentAndCalXCorr(String sequence, short linkSite, double additionalMass, int precursorCharge, SparseVector xcorrPL) { // there are n and c in the sequence
        linkSite = (short) Math.max(1, linkSite);

        int localMaxCharge = Math.min(6, Math.max(precursorCharge - 1, 1));
        double[] inverseChargeArray = new double[localMaxCharge];
        for (int charge = 1; charge <= localMaxCharge; ++charge) {
            inverseChargeArray[charge - 1] = (double) 1 / (double) charge;
        }

        AA[] aaArray = seqToAAList(sequence);

        double xcorr = 0;

        // traverse the sequence to get b-ion
        double bIonMass = massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass; // add N-term modification
        for (int i = 1; i < aaArray.length - 2; ++i) {
            bIonMass += massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            if (i < linkSite) {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin(bIonMass * inverseCharge + PROTON));
                }
            } else {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin((bIonMass + additionalMass) * inverseCharge + PROTON));
                }
            }
        }
        // calculate the last b-ion with C-term modification
        bIonMass +=  massTable.get(aaArray[aaArray.length - 2].aa) + aaArray[aaArray.length - 2].ptmDeltaMass + massTable.get(aaArray[aaArray.length - 1].aa) + aaArray[aaArray.length - 1].ptmDeltaMass;
        for (double inverseCharge : inverseChargeArray) {
            xcorr += xcorrPL.get(mzToBin((bIonMass + additionalMass) * inverseCharge + PROTON)); // for the fragment containing all amino acids, the additional mass is always included.
        }

        // traverse the sequence with reversed order to get y-ion
        // the whole sequence
        double yIonMass = bIonMass + H2O;
        for (double inverseCharge : inverseChargeArray) {
            xcorr += xcorrPL.get(mzToBin((yIonMass + additionalMass) * inverseCharge + PROTON)); // for the fragment containing all amino acids, the additional mass is always included.
        }
        // delete the first amino acid and N-term modification
        yIonMass -= massTable.get(aaArray[0].aa) + aaArray[0].ptmDeltaMass + massTable.get(aaArray[1].aa) + aaArray[1].ptmDeltaMass;
        if (1 >= linkSite) {
            for (double inverseCharge : inverseChargeArray) {
                xcorr += xcorrPL.get(mzToBin(yIonMass * inverseCharge + PROTON));
            }
        } else {
            for (double inverseCharge : inverseChargeArray) {
                xcorr += xcorrPL.get(mzToBin((yIonMass + additionalMass) * inverseCharge + PROTON));
            }
        }
        // rest of the sequence
        for (int i = 2; i < aaArray.length - 2; ++i) {
            yIonMass -= massTable.get(aaArray[i].aa) + aaArray[i].ptmDeltaMass;
            if (i >= linkSite) { // caution: here, it is different from b-ion
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin(yIonMass * inverseCharge + PROTON));
                }
            } else {
                for (double inverseCharge : inverseChargeArray) {
                    xcorr += xcorrPL.get(mzToBin((yIonMass + additionalMass) * inverseCharge + PROTON));
                }
            }
        }

        return xcorr * 0.005;
    }
    // End of cross-linking part

    private Map<Integer, List<int[]>> digest(String proteinSequence) {
        // Cut a protein
        TreeSet<Integer> cutPointSet = new TreeSet<>();
        int length = proteinSequence.length();
        int idxStart = 0;
        Matcher matchObj = digestSitePattern.matcher(proteinSequence);
        cutPointSet.add(0);
        while (idxStart < length) {
            if (matchObj.find()) {
                int cutPoint;
                if (cleavageFromCTerm) {
                    cutPoint = matchObj.end();
                } else {
                    cutPoint = matchObj.start();
                }
                cutPointSet.add(cutPoint);
                idxStart = cutPoint;
            } else {
                cutPointSet.add(length);
                break;
            }
        }

        Integer[] cutPointArray = cutPointSet.toArray(new Integer[0]);

        // Deal with missed cleavage
        Map<Integer, List<int[]>> digestRangeMap = new HashMap<>();
        for (int time = 0; time <= missedCleavage; ++time) {
            List<int[]> temp = new LinkedList<>();
            int leftPoint;
            int rightPoint;
            for (int i = 0; i + 1 + time < cutPointSet.size(); ++i) {
                leftPoint = cutPointArray[i];
                rightPoint = cutPointArray[i + 1 + time];
                temp.add(new int[]{leftPoint, rightPoint});
            }
            digestRangeMap.put(time, temp);
        }

        return digestRangeMap;
    }

    private static String deleteLeftRightFlankingAddNC(String peptide) {
        peptide = leftFlankPattern.matcher(peptide).replaceAll("");
        peptide = rightFlankPattern.matcher(peptide).replaceAll("");
        if (!peptide.startsWith("n")) {
            peptide = "n" + peptide;
        }
        if (!peptide.endsWith("c")) {
            peptide = peptide + "c";
        }
        return peptide;
    }
}
