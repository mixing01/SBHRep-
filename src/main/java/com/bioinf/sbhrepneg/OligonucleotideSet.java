package com.bioinf.sbhrepneg;

import java.util.ArrayList;

public class OligonucleotideSet {
    public final int sequenceLength;
    public final int probeLength;
    public final String startOligoValue;

    public final ArrayList<Oligonucleotide> oligonucleotides;

    public OligonucleotideSet(int n, int k, String s) {
        sequenceLength = n;
        probeLength = k;
        startOligoValue = s;
        oligonucleotides = new ArrayList<>();
    }

    public void addOligonucleotide(Oligonucleotide o) {
        oligonucleotides.add(o);
    }
}
