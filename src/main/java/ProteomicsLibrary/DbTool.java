package ProteomicsLibrary;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class DbTool {

    private Map<String, String> proteinSequenceMap = new HashMap<>();
    private Map<String, String> proteinAnnotateMap = new HashMap<>();

    public DbTool(String dbName, String databaseType) throws IOException {
        String id = "";
        String annotate;
        StringBuilder sequence = new StringBuilder(99999);

        boolean newPro = true;

        Pattern headerPattern;
        if (databaseType.contentEquals("TAIR")) {
            headerPattern = Pattern.compile("^>([^\\s]+)[\\s|]+(.+)$");
        } else if (databaseType.contentEquals("UniProt") || databaseType.contentEquals("SwissProt")) {
            headerPattern = Pattern.compile("^>[^|]+\\|(.+)\\|(.+)$");
        } else if (databaseType.contentEquals("neXtProt")) {
            headerPattern = Pattern.compile("^>nxp:NX_([^ ]+) (.+)");
        } else if (databaseType.contentEquals("contaminants") || databaseType.contentEquals("ITAG")) {
            headerPattern = Pattern.compile("^>([^ ]+) (.+)$");
        } else if (databaseType.contentEquals("Others")) {
            headerPattern = Pattern.compile("^>(.+)$");
        } else {
            throw new NullPointerException(String.format(Locale.US, "Incorrect database type (%s) in the parameter file.", databaseType));
        }

        BufferedReader dbReader;
        if (databaseType.contentEquals("contaminants")) {
            InputStream inputStream = getClass().getClassLoader().getResourceAsStream("contaminants.fasta");
            dbReader = new BufferedReader(new InputStreamReader(inputStream));
        } else {
            dbReader = new BufferedReader(new FileReader(dbName));
        }
        String line;
        while ((line = dbReader.readLine()) != null) {
            line = line.trim();
            Matcher headMatcher = headerPattern.matcher(line);
            if (headMatcher.matches()) {
                // This line is a header
                if (!newPro) {
                    // This isn't the first protein
                    proteinSequenceMap.put(id, sequence.toString());
                }
                id = headMatcher.group(1).trim();
                if (databaseType.contentEquals("Others")) {
                    annotate = id;
                } else {
                    annotate = headMatcher.group(2).trim();
                }
                proteinAnnotateMap.put(id, annotate);
                newPro = true;
            } else if (!line.isEmpty()) {
                // This line is a body
                if (newPro) {
                    sequence = new StringBuilder(99999);
                    sequence.append(line);
                    newPro = false;
                } else {
                    sequence.append(line);
                }
            }
        }
        dbReader.close();
        // Last protein
        proteinSequenceMap.put(id, sequence.toString());
    }

    public Map<String, String> getProteinSequenceMap() {
        return proteinSequenceMap;
    }

    public Map<String, String> getProteinAnnotateMap() {
        return proteinAnnotateMap;
    }

    public static Set<Integer> findPeptideLocation(String proteinSequence, String peptide, String cutSite, String protectSite) throws NullPointerException {
        peptide = getSequenceOnly(peptide.trim());
        Set<Integer> output = new HashSet<>();
        int idx = proteinSequence.indexOf(peptide);
        while (idx >= 0) {
            if ((idx == 0 || cutSite.contains(proteinSequence.substring(idx - 1, idx)) || (idx == 1 && proteinSequence.charAt(0) == 'M')) && (idx + peptide.length() == proteinSequence.length() || !protectSite.contains(proteinSequence.substring(idx + peptide.length(), idx + peptide.length() + 1)))) { // caution: we only consider cutting from N-term.
                output.add(idx);
            }
            idx = proteinSequence.indexOf(peptide, idx + 1);
        }
        return output;
    }

    public static TreeSet<String> reduceProteinIdSet(Set<String> input, String databaseType) { // this only works for TAIR
        if (input.size() == 1 || !databaseType.contentEquals("TAIR")) {
            return new TreeSet<>(input);
        } else {
            Map<String, Integer> tempMap = new HashMap<>();
            for (String s : input) {
                String[] tempArray = s.split("\\.");
                if (tempMap.containsKey(tempArray[0])) {
                    if (tempMap.get(tempArray[0]) > Integer.valueOf(tempArray[1])) {
                        tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                    }
                } else {
                    tempMap.put(tempArray[0], Integer.valueOf(tempArray[1]));
                }
            }
            TreeSet<String> output = new TreeSet<>();
            for (String s : tempMap.keySet()) {
                output.add(s + "." + tempMap.get(s));
            }
            return output;
        }
    }

    public static String getPtmFreePeptide(String peptide) {
        return peptide.replaceAll("[^A-Znc]+", "");
    }

    public static String getSequenceOnly(String peptide) {
        return peptide.replaceAll("[^A-Z]+", "");
    }
}
