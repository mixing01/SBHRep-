package com.bioinf.sbhrepneg;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.IOException;
import java.net.URI;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;

public class InstanceReader {
    /**
     * Returns an OligonucleotideSet read from a remote generator.
     * <p>
     * Method sends an HTTP request to www.<!-- -->cs.<!-- -->put.<!-- -->poznan.<!-- -->pl/pwawrzyniak/bio/bio.<!-- -->php.
     * The generated oligonucleotides can be formed into a DNA sequence.
     * The oligonucleotides are supplied with an information about how many times they could be repeated in a sequence.
     * Parameters specify the sequence that the oligonucleotides can be formed into, the length of a
     * single oligonucleotide and the number of negative errors in a sequence.
     *
     * @param n    sequence length
     * @param k    oligonucleotide length
     * @param sqne number of errors
     * @return OligonucleotideSet containing read Oligonucleotides
     * @throws IOException                     If an input or output exception occured
     * @throws TooLongOligonucleotideException If the length of single oligonucleotide is larger than sequence length
     * @throws ErrorNumberException            If the number of errors is not between 0 and n/4
     * @see Oligonucleotide
     * @see OligonucleotideSet
     */
    public OligonucleotideSet readInstanceFromXML(int n, int k, int sqne) throws Exception {
        /*
          [number of matches] -> [intensity signal]
          0  ->  0
          1  -> 1-3
          2  -> 3-5
          3  -> 5-7
          4  -> 6-8
          5  -> 7-8
          6  -> 7-9
          7  -> 8-9
          8  -> 8-9
          9+ ->  9
         */

        // Exception handling
        if (n < 1)
            return null;
        if (k > n)
            throw new TooLongOligonucleotideException("Oligonucleotide length cannot be larger than sequence length");
        if (sqne > n / 4)
            throw new ErrorNumberException("Number of errors cannot exceed n/4");
        if (sqne < 0)
            throw new ErrorNumberException("Number of errors cannot be negative");


        // Creating a URI
        String uri = String.format("https://www.cs.put.poznan.pl/pwawrzyniak/bio/bio.php?n=%d&k=%d&mode=basic&intensity=1&position=0&sqpe=0&sqne=%d&pose=0", n, k, sqne);

        // XML Parsing from URI
        DocumentBuilderFactory documentBuilderFactory = DocumentBuilderFactory.newInstance();
        documentBuilderFactory.setFeature("http://xml.org/sax/features/external-general-entities", false);
        documentBuilderFactory.setFeature("http://xml.org/sax/features/external-parameter-entities", false);
        DocumentBuilder builder = documentBuilderFactory.newDocumentBuilder();
        Document doc = builder.parse(uri);
        Element dna = (Element) doc.getElementsByTagName("dna").item(0);
        Element probe = (Element) dna.getElementsByTagName("probe").item(0);
        NodeList cells = probe.getElementsByTagName("cell");

        // Creating OligonucleotideSet
        String startOligoValue = dna.getAttribute("start");
        OligonucleotideSet oligonucleotideSet = new OligonucleotideSet(n, k, startOligoValue);
        for (int i = 0; i < cells.getLength(); i++) {
            Element cell = (Element) cells.item(i);
            int rep = calculateReps(Integer.parseInt(cell.getAttribute("intensity")));
            oligonucleotideSet.addOligonucleotide(new Oligonucleotide(cell.getChildNodes().item(0).getNodeValue(), rep));
        }
        return oligonucleotideSet;
    }



    private int calculateReps(int intensity) {
        switch (intensity) {
            case 1, 2, 3 -> {
                return 1;
            }
            case 4, 5 -> {
                return 2;
            }
            case 6, 7 -> {
                return 3;
            }
            case 8 -> {
                return 4;
            }
            case 9 -> {
                return 6;
            }
            default -> {
                return 0;
            }
        }
    }

    public String get(String uri) throws Exception {
        HttpClient client = HttpClient.newHttpClient();
        HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(uri))
                .build();

        HttpResponse<String> response =
                client.send(request, HttpResponse.BodyHandlers.ofString());

        return response.body();
    }
}
