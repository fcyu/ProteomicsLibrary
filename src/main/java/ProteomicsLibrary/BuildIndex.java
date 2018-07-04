package ProteomicsLibrary;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.util.Map;
import java.util.Set;

public class BuildIndex {

    private Multimap<String, String> targetPeptideProteinMap = HashMultimap.create();

    public BuildIndex(Map<String, String> proteinSequenceMap, String cleavageSite1, String protectionSite1, boolean cleavageFromCTerm1, String cleavageSite2, String protectionSite2, Boolean cleavageFromCTerm2, int missedCleavage) {
        MassTool massTool = new MassTool(missedCleavage, cleavageSite1, protectionSite1, cleavageFromCTerm1, cleavageSite2, protectionSite2, cleavageFromCTerm2, 0.02, 1, "N14");
        for (String protein : proteinSequenceMap.keySet()) {
            String proteinSequence = proteinSequenceMap.get(protein);
            Set<String> peptideSet = massTool.buildPeptideSet(proteinSequence);
            if (proteinSequence.startsWith("M")) {
                peptideSet.addAll(massTool.buildPeptideSet(proteinSequence.substring(1)));
            }

            for (String peptide : peptideSet) {
                targetPeptideProteinMap.put(peptide, protein);
            }
        }
    }

    public Multimap<String, String> getTargetPeptideProteinMap() {
        return targetPeptideProteinMap;
    }
}
