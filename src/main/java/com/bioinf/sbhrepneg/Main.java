package com.bioinf.sbhrepneg;

import java.util.ArrayList;

public class Main {
    public static void main(String[] args) throws Exception {
        InstanceReader instanceReader = new InstanceReader();
       ArrayList<OligonucleotideSet> testSets = new ArrayList<>();
       for(int i = 100; i <= 500; i+=100) {
           testSets.add(instanceReader.readInstanceFromXML(i,8,i/10));
       }
        System.nanoTime();
       long timeStart;
       long timeStop;
       String result;
        for(OligonucleotideSet s : testSets) {
            timeStart = System.nanoTime();
            result = Genetic.calculate(s, 20, 10, 10000, 500, 0.01);
            timeStop = System.nanoTime();
            System.out.println((timeStop-timeStart)+" "+result.length());

            timeStart = System.nanoTime();
            result = Genetic.calculate(s, 20, 50, 10000, 500, 0.01);
            timeStop = System.nanoTime();
            System.out.println((timeStop-timeStart)+" "+result.length());

            timeStart = System.nanoTime();
            result = Genetic.calculate(s, 20, 100, 10000, 500, 0.01);
            timeStop = System.nanoTime();
            System.out.println((timeStop-timeStart)+" "+result.length());

            timeStart = System.nanoTime();
            result = Genetic.calculate(s, 20, 200, 10000, 500, 0.01);
            timeStop = System.nanoTime();
            System.out.println((timeStop-timeStart)+" "+result.length());

            timeStart = System.nanoTime();
            result = Genetic.calculate(s, 20, 300, 10000, 500, 0.01);
            timeStop = System.nanoTime();
            System.out.println((timeStop-timeStart)+" "+result.length());

        }
    }
}
