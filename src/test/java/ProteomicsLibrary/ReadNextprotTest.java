package ProteomicsLibrary;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.junit.Before;
import org.junit.Test;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;

public class ReadNextprotTest {

    private static ReadNextprot readNextprot;

    @Before
    public void setUp() throws Exception {
        readNextprot = new ReadNextprot(Thread.currentThread().getContextClassLoader().getResource("neXtProt.test.peff").getPath());
    }

    @Test
    public void getIdEntryMap() {
        Map<String, ReadNextprot.Entry> results = readNextprot.getIdEntryMap();
        assertEquals(1, results.size());
        ReadNextprot.Entry entry = results.values().iterator().next();
        Multimap<Integer, Character> locationVariantMap = HashMultimap.create();
        locationVariantMap.put(3, '*');
        locationVariantMap.put(8, 'L');
        locationVariantMap.put(12, 'F');
        locationVariantMap.put(12, 'P');
        Multimap<Integer, String> locationModMap = HashMultimap.create();
        locationModMap.put(42, "Disulfide");
        locationModMap.put(112, "Disulfide");
        Multimap<Integer, String> locationPsiModMap = HashMultimap.create();
        locationPsiModMap.put(318, "MOD:00115");
        locationPsiModMap.put(325, "MOD:00115");
        assertEquals("A0A075B6H9-1", entry.id);
        assertEquals("MAWTPLLFLTLLLHCTGSLSQLVLTQSPSASASLGASVKLTCTLSSGHSSYAIAWHQQQPEKGPRYLMKLNSDGSHSKGDGIPDRFSGSSSGAERYLTISSLQSEDEADYYCQTWGTGI", entry.sequence);
        assertEquals(locationVariantMap.size(), entry.locationVariantMap.size());
        for (int location : locationVariantMap.keySet()) {
            assertTrue(locationVariantMap.get(location).containsAll(entry.locationVariantMap.get(location)));
            assertTrue(entry.locationVariantMap.get(location).containsAll(locationVariantMap.get(location)));
        }
        assertEquals(locationModMap.size(), entry.locationModMap.size());
        for (int location : locationModMap.keySet()) {
            assertTrue(entry.locationModMap.get(location).containsAll(locationModMap.get(location)));
            assertTrue(locationModMap.get(location).containsAll(entry.locationModMap.get(location)));
        }
        assertEquals(locationPsiModMap.size(), entry.locationPsiModMap.size());
        for (int location : locationPsiModMap.keySet()) {
            assertTrue(entry.locationPsiModMap.get(location).containsAll(locationPsiModMap.get(location)));
            assertTrue(locationPsiModMap.get(location).containsAll(entry.locationPsiModMap.get(location)));
        }
    }
}