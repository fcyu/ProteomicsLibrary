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
            int i = 0;
            while (proteinSequence.charAt(i) == 'M') {
                peptideSet.addAll(massTool.buildPeptideSet(proteinSequence.substring(++i)));
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
