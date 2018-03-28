package ProteomicsLibrary;

import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Utilities {

    private static final Pattern mgfScanNumPattern1 = Pattern.compile("Scan:([0-9]+) ", Pattern.CASE_INSENSITIVE);
    private static final Pattern mgfScanNumPattern2 = Pattern.compile("scan=([0-9]+)", Pattern.CASE_INSENSITIVE);
    private static final Pattern mgfScanNumPattern3 = Pattern.compile("scan ([0-9]+)", Pattern.CASE_INSENSITIVE);
    private static final Pattern mgfScanNumPattern4 = Pattern.compile(".+ [0-9]+, \\+MS2\\([0-9.]+\\), [0-9. ]+eV, [0-9. ]+min #([0-9]+)$", Pattern.CASE_INSENSITIVE);
    private static final Pattern mgfScanNumPattern5 = Pattern.compile("^[^.]+\\.([0-9]+)\\.[0-9]+\\.[0-9]");
    private static final Pattern mgfScanNumPattern6 = Pattern.compile("\\.([0-9]+)\\.[0-9]+\\.");

    public static final Pattern csvSplitPattern = Pattern.compile(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");
    public static final Pattern tsvSplitPattern = Pattern.compile("\t(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)");

    public static int getScanNum(String s) {
        Matcher matcher1 = mgfScanNumPattern1.matcher(s);
        if (matcher1.find()) {
            return Integer.valueOf(matcher1.group(1));
        } else {
            Matcher matcher2 = mgfScanNumPattern2.matcher(s);
            if (matcher2.find()) {
                return Integer.valueOf(matcher2.group(1));
            } else {
                Matcher matcher3 = mgfScanNumPattern3.matcher(s);
                if (matcher3.find()) {
                    return Integer.valueOf(matcher3.group(1));
                } else {
                    Matcher matcher4 = mgfScanNumPattern4.matcher(s);
                    if (matcher4.find()) {
                        return Integer.valueOf(matcher4.group(1));
                    } else {
                        Matcher matcher5 = mgfScanNumPattern5.matcher(s);
                        if (matcher5.find()) {
                            return Integer.valueOf(matcher5.group(1));
                        } else {
                            Matcher matcher6 = mgfScanNumPattern6.matcher(s);
                            if (matcher6.find()) {
                                return Integer.valueOf(matcher6.group(1));
                            } else {
                                throw new NullPointerException(String.format(Locale.US, "Cannot get scan number from the MGF title %s. Please report your MGF title to the author.", s));
                            }
                        }
                    }
                }
            }
        }
    }
}
