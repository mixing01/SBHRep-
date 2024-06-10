package com.bioinf.sbhrepneg;

public class Oligonucleotide {
    public final String value;
    private int repMin;

    public Oligonucleotide(String value, int repMin) {
        this.value = value;
        this.repMin = repMin;
    }

    public int checkOffset(Oligonucleotide oligo2) {
        String value2 = oligo2.value;
        int length = value.length();
        for(int i = 1; i<length; i++) {
            if(value.substring(i).equals(value2.substring(0,length-i)))
                return i;
        }
        return length;
    }

    public boolean compare(Oligonucleotide oligo2) {
        return value.equals(oligo2.value);
    }

    public void decrementRep() {
        repMin--;
    }

    public boolean checkRep() {
        return repMin > 1;
    }

    public int getRepMin() {
        return repMin;
    }
}
