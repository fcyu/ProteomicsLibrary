package ProteomicsLibrary.Types;

import java.util.Locale;

public class AA {

    public char aa;
    public double ptmDeltaMass;
    private final int hashCode;

    public AA(char aa, double ptmDeltaMass) {
        this.aa = aa;
        this.ptmDeltaMass = ptmDeltaMass;
        hashCode = String.format(Locale.US, "%c(%.3f)", aa, ptmDeltaMass).hashCode();
    }

    public boolean hasMod() {
        return Math.abs(ptmDeltaMass) > 0.1;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        if (other instanceof AA) {
            return ((AA) other).hashCode == hashCode;
        } else {
            return false;
        }
    }

    public AA clone() {
        return new AA(aa, ptmDeltaMass);
    }

    public String toString() {
        if (hasMod()) {
            return String.format(Locale.US, "%c(%.3f)", aa, ptmDeltaMass);
        } else {
            return String.valueOf(aa);
        }
    }
}
