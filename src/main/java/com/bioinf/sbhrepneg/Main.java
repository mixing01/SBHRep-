package com.bioinf.sbhrepneg;

public class Main {
    public static void main(String[] args) throws Exception {
        InstanceReader instanceReader = new InstanceReader();
        OligonucleotideSet oligonucleotideSet = instanceReader.readInstanceFromXML(1000,10,20);
        System.out.println("Start: "+oligonucleotideSet.startOligoValue);
        for(Oligonucleotide o : oligonucleotideSet.oligonucleotides) {
            System.out.println(o.value);
        }
        //System.out.println(Greedy.calculate(oligonucleotideSet));
        //System.out.println(Exact.calculate(oligonucleotideSet));
        System.out.println(Genetic.calculate(oligonucleotideSet, 10, 50,1000,500,0.01).length());
    }
}
