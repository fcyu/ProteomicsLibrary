package ProteomicsLibrary;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReadNextprot {

    private static final Pattern variantPattern = Pattern.compile("\\\\VariantSimple=([^=\\\\ ]+)");
    private static final Pattern variantPattern2 = Pattern.compile("([0-9]+)\\|([A-Z*])");
    private static final Pattern psiModPattern = Pattern.compile("\\\\ModResPsi=([^=\\\\]+)");
    private static final Pattern psiModPattern2 = Pattern.compile("([0-9]+)\\|(MOD:[^|()]+)");
    private static final Pattern modPattern = Pattern.compile("\\\\ModRes=([^=\\\\]+)");
    private static final Pattern modPattern2 = Pattern.compile("([0-9]+)\\|\\|([^|()]+)");

    private Map<String, Entry> idEntryMap = new HashMap<>();

    public ReadNextprot(String peffPath) throws Exception {
        DbTool dbTool = new DbTool(peffPath, "neXtProt");
        Map<String, String> proteinAnnotationMap = dbTool.getProteinAnnotateMap();
        Map<String, String> proteinSequenceMap = dbTool.getProteinSequenceMap();
        for (String id : proteinAnnotationMap.keySet()) {
            Matcher variantMatcher = variantPattern.matcher(proteinAnnotationMap.get(id));
            Matcher psiModMatcher = psiModPattern.matcher(proteinAnnotationMap.get(id));
            Matcher modMatcher = modPattern.matcher(proteinAnnotationMap.get(id));
            Multimap<Integer, Character> locationVariantMap = HashMultimap.create();
            Multimap<Integer, String> locationPsiModMap = HashMultimap.create();
            Multimap<Integer, String> locationModMap = HashMultimap.create();
            if (variantMatcher.find()) {
                String tempStr = variantMatcher.group(1);
                Matcher variantMatcher2 = variantPattern2.matcher(tempStr);
                while (variantMatcher2.find()) {
                    int location = Integer.valueOf(variantMatcher2.group(1).trim()); // stars from 1
                    if (variantMatcher2.group(2).trim().length() != 1) {
                        throw new Exception(String.format(Locale.US, "The substituted length is larger than one in protein %s.", id));
                    }
                    char aa2 = variantMatcher2.group(2).trim().charAt(0);
                    locationVariantMap.put(location, aa2);
                }
            }
            if (psiModMatcher.find()) {
                String tempStr = psiModMatcher.group(1);
                Matcher psiModMatcher2 = psiModPattern2.matcher(tempStr);
                while (psiModMatcher2.find()) {
                    int location = Integer.valueOf(psiModMatcher2.group(1).trim()); // starts from 1
                    String psiMod = psiModMatcher2.group(2).trim();
                    locationPsiModMap.put(location, psiMod);
                }
            }
            if (modMatcher.find()) {
                String tempSet = modMatcher.group(1);
                Matcher modMatcher2 = modPattern2.matcher(tempSet);
                while (modMatcher2.find()) {
                    int location = Integer.valueOf(modMatcher2.group(1).trim());
                    String mod = modMatcher2.group(2).trim();
                    locationModMap.put(location, mod);
                }
            }
            if (!locationVariantMap.isEmpty() || !locationPsiModMap.isEmpty()) {
                if (idEntryMap.containsKey(id)) {
                    throw new Exception(String.format(Locale.US, "There are duplicate protein IDs: %s.", id));
                }
                idEntryMap.put(id, new Entry(id, proteinSequenceMap.get(id), locationVariantMap, locationModMap, locationPsiModMap));
            }
        }
    }

    public Map<String, Entry> getIdEntryMap() {
        return idEntryMap;
    }

    public class Entry {

        public final String id;
        public final String sequence;
        public final Multimap<Integer, Character> locationVariantMap;
        public final Multimap<Integer, String> locationModMap;
        public final Multimap<Integer, String> locationPsiModMap;

        public Entry(String id, String sequence, Multimap<Integer, Character> locationVariantMap, Multimap<Integer, String> locationModMap, Multimap<Integer, String> locationPsiModMap) {
            this.id = id;
            this.sequence = sequence;
            this.locationVariantMap = locationVariantMap;
            this.locationModMap = locationModMap;
            this.locationPsiModMap = locationPsiModMap;
        }
    }
}
