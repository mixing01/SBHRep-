package com.bioinf.sbhrepneg;

import java.util.ArrayList;

public class Greedy {

    private Greedy(){}
    public static String calculate(OligonucleotideSet initialSet) {
        StringBuilder result = new StringBuilder();
        result.append(initialSet.startOligoValue);
        int probeLength = initialSet.probeLength;
        ArrayList<Oligonucleotide> remainingOligo = new ArrayList<>(initialSet.oligonucleotides);
        int startingOligoIndex = findOligo(remainingOligo, initialSet.startOligoValue);
        Oligonucleotide lastAddedOligo;
        if(startingOligoIndex == -1) {
            lastAddedOligo = new Oligonucleotide(initialSet.startOligoValue, 1);
        }
        else {
            lastAddedOligo = remainingOligo.get(startingOligoIndex);
            removeOrDecrementNucleotide(remainingOligo, lastAddedOligo);
        }
        int targetLength = initialSet.sequenceLength;
        while (result.length() < targetLength) {
            int[] scores = new int[remainingOligo.size()];
            for (int i = 0; i < remainingOligo.size(); i++) {
                ArrayList<Oligonucleotide> tempSet = new ArrayList<>(remainingOligo);
                Oligonucleotide currOligo = remainingOligo.get(i);
                scores[i] += lastAddedOligo.checkOffset(currOligo);
                removeOrDecrementNucleotide(tempSet, currOligo);
                int min = probeLength + 1;
                for (Oligonucleotide tempOligo : tempSet) {
                    int tempScore = currOligo.checkOffset(tempOligo);
                    if (tempScore < min) {
                        min = tempScore;
                    }
                }
                scores[i] += min;
            }
            int minScore = scores[0];
            int minIndex = 0;
            for (int i = 1; i < scores.length; i++) {
                if(scores[i]<minScore) {
                    minScore = scores[i];
                    minIndex = i;
                }
            }
            Oligonucleotide newOligo = remainingOligo.get(minIndex);
            result.append(newOligo.value.substring(probeLength - lastAddedOligo.checkOffset(newOligo)));
            if(result.length() > targetLength)
                return result.substring(0,result.length() - probeLength);
            removeOrDecrementNucleotide(remainingOligo, minIndex);
        }
        return result.toString();
    }

    private static void removeOrDecrementNucleotide(ArrayList<Oligonucleotide> set, Oligonucleotide o) {
        for (int i = 0; i < set.size(); i++) {
            Oligonucleotide temp = set.get(i);
            if (temp.compare(o)) {
                if (temp.checkRep()) {
                    temp.decrementRep();
                    set.set(i, temp);
                } else {
                    set.remove(i);
                }
                break;
            }

        }
    }

    private static void removeOrDecrementNucleotide(ArrayList<Oligonucleotide> set, int i) {
        Oligonucleotide temp = set.get(i);
        if (temp.checkRep()) {
            temp.decrementRep();
            set.set(i, temp);
        } else {
            set.remove(i);
        }
    }

    private static int findOligo(ArrayList<Oligonucleotide> set, String value) {
        for (int i = 0; i < set.size(); i++) {
            if (set.get(i).value.equals(value))
                return i;
        }
        return -1;
    }
}
