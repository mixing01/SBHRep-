package com.bioinf.sbhrepneg;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

public class Genetic {
    static Random random = new Random();

    private Genetic(){}

    public static String calculate(OligonucleotideSet set, int popNum, int epochLim, int totalEpochLim, int multiplier, double mutationProb) throws InvalidProbabilityException {
        if(mutationProb < 0 || mutationProb > 1) {
            throw new InvalidProbabilityException("Mutation probability should be between 0 and 1 inclusive");
        }
        ArrayList<Oligonucleotide> initList = set.oligonucleotides;
        ArrayList<String> initListString = new ArrayList<>();

        for(Oligonucleotide o : initList) {
            for(int i = 0; i<o.getRepMin(); i++) {
                initListString.add(o.value);
            }
        }
        HashMap<String, Integer> repNums = new HashMap<>();
        for(Oligonucleotide e : set.oligonucleotides) {
            repNums.put(e.value, e.getRepMin());
        }

        ArrayList<ArrayList<String>> population = new ArrayList<>();
        for(int i = 0; i < popNum; i++) {
            ArrayList<String> newChromosome = generateRandomChromosome(initListString, set.startOligoValue, set.sequenceLength, set.probeLength);
            population.add(newChromosome);
        }
        int initPopSize = population.size();
        int popSizeEven = (initPopSize % 2 == 0) ? initPopSize : initPopSize - 1;
        int bestScore = getFitnessScore(population.get(0),repNums,multiplier,set.sequenceLength, set.probeLength);
        int totalEpoch = 0;
        ArrayList<String> bestChromosome = population.get(0);

        for(int epoch = 0; epoch < epochLim; epoch++) {
            totalEpoch++;
            if(totalEpoch > totalEpochLim)
                break;

            for(ArrayList<String> c : population) {
                int currFitness = getFitnessScore(c,repNums,multiplier,set.sequenceLength, set.probeLength);
                if(currFitness < bestScore) {
                    bestScore = currFitness;
                    bestChromosome = c;
                    epoch = 0;
                }
            }

            ArrayList<ArrayList<String>> children = new ArrayList<>();
            for(int i = 0; i<popSizeEven; i+=2) {
                children.add(crossing(population.get(i),population.get(i+1)));
            }
            population.addAll(children);

            for(int i = 0; i<population.size(); i++) {
                if(random.nextInt(0, (int) (1/mutationProb)) == 0) {
                    population.set(i,mutation(population.get(i)));
                }
            }

            tournament(population,repNums,multiplier,set.sequenceLength, set.probeLength);

            while(population.size() < initPopSize) {
                ArrayList<String> newChromosome = generateRandomChromosome(initListString, set.startOligoValue, set.sequenceLength, set.probeLength);
                population.add(newChromosome);
            }
        }

        return generateResult(bestChromosome);
    }

    private static ArrayList<String> generateRandomChromosome(ArrayList<String> initList, String firstOligo, int length, int probeSize) {
        ArrayList<String> newList = new ArrayList<>();
        ArrayList<String> oldList = new ArrayList<>(initList);
        newList.add(firstOligo);
        oldList.remove(firstOligo);
        while(checkTotalOffset(newList) < length - probeSize && oldList.size() > 0) {
            int randomPlace = random.nextInt(0, oldList.size());
            newList.add(oldList.get(randomPlace));
            oldList.remove(randomPlace);
        }
        return newList;
    }

    private static int getFitnessScore(ArrayList<String> chromosome, HashMap<String, Integer> repNums, int multiplier, int length, int probeLength) {
        int score = 0;
        int offsets = 0;
        HashMap<String, Integer> localRepNums = new HashMap<>();
        for(String val : chromosome) {
            localRepNums.merge(val, 1, Integer::sum);
        }
        for(String key : repNums.keySet()) {
            if(localRepNums.containsKey(key))
                score += Math.abs(localRepNums.get(key) - repNums.get(key))*multiplier;
            else
                score += repNums.get(key)*multiplier;
        }
        for(int i = 0; i < chromosome.size()-1; i++) {
            offsets+=offsetStrings(chromosome.get(i),chromosome.get(i+1));
        }
        if(offsets > length - probeLength) {
            score+=offsets*length*multiplier;
        }
        else if (offsets < length - probeLength) {
            score+=offsets*multiplier;
        }
        else {
            score += offsets;
        }
        return score;
    }

    private static int offsetStrings(String s1, String s2) {
        int length = s1.length();
        for(int i = 1; i<length; i++) {
            if(s1.substring(i).equals(s2.substring(0,length-i)))
                return i;
        }
        return length;
    }

    private static ArrayList<String> crossing(ArrayList<String> c1, ArrayList<String> c2) {
        int splitPoint = random.nextInt(1,c1.size());
        ArrayList<String> newChromosome = new ArrayList<>();
        for(int i = 0; i<splitPoint; i++) {
            newChromosome.add(c1.get(i));
        }
        for(int i = splitPoint; i<c2.size(); i++) {
            newChromosome.add(c2.get(i));
        }
        return newChromosome;
    }

    private static ArrayList<String> mutation(ArrayList<String> c) {
        ArrayList<String> result = new ArrayList<>(c);
        int mutationPoint1 = random.nextInt(1,c.size());
        int mutationPoint2 = random.nextInt(1,c.size());
        String temp = c.get(mutationPoint1);
        result.set(mutationPoint1,c.get(mutationPoint2));
        result.set(mutationPoint2,temp);
        return result;
    }

    private static void tournament(ArrayList<ArrayList<String>> population, HashMap<String, Integer> repNums, int multiplier, int length, int probeLength) {
        Collections.shuffle(population);
        ArrayList<Integer> indicesToRemove = new ArrayList<>();
        int analizedLength = (population.size() % 2 == 0) ? population.size() : population.size() - 1;
        for(int i = 0; i<analizedLength; i+=2) {
            int indexToRemove = (getFitnessScore(population.get(i),repNums, multiplier, length, probeLength) > getFitnessScore(population.get(i+1),repNums, multiplier, length, probeLength)) ? i : i+1;
            indicesToRemove.add(indexToRemove);
        }
        for(int i = indicesToRemove.size()-1; i >= 0; i--) {
            population.remove((int) indicesToRemove.get(i));
        }
    }

    private static String generateResult(ArrayList<String> c) {
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(c.get(0));
        for(int i = 0; i<c.size()-1; i++) {
            stringBuilder.append(getOffset(c.get(i),c.get(i+1)));
        }
        return stringBuilder.toString();
    }

    private static int checkTotalOffset(ArrayList<String> c) {
        int result = 0;
        for(int i = 0; i < c.size()-1; i++) {
            result += offsetStrings(c.get(0),c.get(1));
        }
        return result;
    }

    private static String getOffset(String c1, String c2) {
        int length = c1.length();
        for(int i = 1; i<length; i++) {
            if(c1.substring(i).equals(c2.substring(0,length-i)))
                return c2.substring(length-i,length);
        }
        return c2;
    }
}
