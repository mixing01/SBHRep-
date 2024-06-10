package com.bioinf.sbhrepneg;

import org.cicirello.permutations.Permutation;
import org.cicirello.permutations.PermutationIterator;

import java.util.ArrayList;

public class Exact {
    private Exact(){}

    public static String calculate(OligonucleotideSet set) {
        StringBuilder result = new StringBuilder();
        int probeLength = set.probeLength;
        ArrayList<Oligonucleotide> currentSet = new ArrayList<>(set.oligonucleotides);
        ArrayList<Oligonucleotide> multiset = new ArrayList<>();
        for (Oligonucleotide curr : currentSet) {
            for (int j = curr.getRepMin(); j > 0; j--) {
                multiset.add(curr);
            }
        }
        int graphSize = multiset.size();
        int[][] graph = new int[graphSize][graphSize];
        for(int i = 0; i < graphSize; i++) {
            for(int j = 0; j < graphSize; j++) {
                if(i!=j)
                    graph[i][j] = multiset.get(i).checkOffset(multiset.get(j));
                else
                    graph[i][j] = -1;
            }
        }
        int startingIndex = findOligo(multiset, set.startOligoValue);
        Oligonucleotide lastAddedOligo;
        if(startingIndex == -1) {
            lastAddedOligo = new Oligonucleotide(set.startOligoValue, 0);
        }
        else {
            lastAddedOligo = multiset.get(startingIndex);
            multiset.remove(startingIndex);
            graphSize--;
        }
        result.append(lastAddedOligo.value);
        PermutationIterator permutationIterator = new PermutationIterator(graphSize);
        int bestScore;
        Permutation bestPermutation;
        Permutation currPermutation = permutationIterator.next();
        int currScore = 0;
        for(int i = 0; i< graphSize - 1; i++) {
            int currI = currPermutation.get(i);
            int currJ = currPermutation.get(i+1);
            currScore += multiset.get(currI).checkOffset(multiset.get(currJ));
        }
        bestScore = currScore;
        bestPermutation = currPermutation;
        while(permutationIterator.hasNext()) {
            currPermutation = permutationIterator.next();
            currScore = 0;
            for(int i = 0; i< graphSize - 1; i++) {
                int currI = currPermutation.get(i);
                int currJ = currPermutation.get(i+1);
                currScore += multiset.get(currI).checkOffset(multiset.get(currJ));
            }
            if(currScore < bestScore) {
                bestScore = currScore;
                bestPermutation = currPermutation;
            }

        }
        for(int i = 1; i<graphSize; i++) {
            int currIndex = bestPermutation.get(i);
            Oligonucleotide currOligo = multiset.get(currIndex);
            result.append(currOligo.value.substring(probeLength - lastAddedOligo.checkOffset(currOligo)));
            lastAddedOligo = currOligo;
        }
        return result.toString();
    }

    private static int findOligo(ArrayList<Oligonucleotide> set, String value) {
        for (int i = 0; i < set.size(); i++) {
            if (set.get(i).value.equals(value))
                return i;
        }
        return -1;
    }
}
