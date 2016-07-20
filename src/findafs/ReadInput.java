/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package findafs;

import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 *
 * @author bmumey
 */
public class ReadInput {

    public static Graph readDotFile(String fileName) {

        Graph g = new Graph();
        TreeMap<Integer, TreeSet<Integer>> nodeNeighbors = new TreeMap<Integer, TreeSet<Integer>>();
        TreeMap<Integer, ArrayList<Integer>> nodeStarts = new TreeMap<Integer, ArrayList<Integer>>();
        TreeMap<Integer, Integer> nodeLength = new TreeMap<Integer, Integer>();
        try {
            Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(fileName)));
            Scanner lineScanner = new Scanner(scanner.nextLine());
            lineScanner.useDelimiter(" ");
            lineScanner.next();
            g.name = lineScanner.next();

            g.maxStart = 0;

            int linesRead = 0;
            String prevLine = "xxx";
            while (scanner.hasNextLine()) {
                String nextLine = scanner.nextLine();
                linesRead++;
                if (linesRead % 20000 == 0) {
                    System.out.println("linesRead: " + linesRead);
                }
                while (nextLine.equals(prevLine) && scanner.hasNextLine()) {
                    nextLine = scanner.nextLine();
                    linesRead++;
                    if (linesRead % 20000 == 0) {
                        System.out.println("linesRead: " + linesRead);
                    }
                }
                prevLine = nextLine;

                if (nextLine.contains("label")) { //node
                    lineScanner = new Scanner(nextLine);
                    int node = lineScanner.nextInt();
                    nodeNeighbors.put(node, new TreeSet<Integer>());
                    String label = lineScanner.next();
                    label = label.split("\"")[1];
                    String[] l = label.split(":");
                    String[] starts = l[0].split(",");
                    nodeStarts.put(node, new ArrayList<Integer>());
                    for (String s : starts) {
                        int start = Integer.parseInt(s);
                        nodeStarts.get(node).add(start);
                        g.maxStart = Math.max(g.maxStart, start);
                    }
                    nodeLength.put(node, Integer.parseInt(l[1]));
                    if (node % 10000 == 0) {
                        System.out.println("reading node: " + node);
                    }
                }
                if (nextLine.contains("->")) { //edge
                    lineScanner = new Scanner(nextLine);
                    lineScanner.useDelimiter("->");
                    int tail = Integer.parseInt(lineScanner.next().trim());
                    int head = Integer.parseInt(lineScanner.next().trim());
                    nodeNeighbors.get(tail).add(head);
                    //System.out.println("edge: " + tail + "->" + head);
                }
            }
        } catch (Exception ex) {
            System.err.println(ex);
            ex.printStackTrace();
            System.exit(-1);
        }

        g.numNodes = nodeNeighbors.keySet().size();
        g.neighbor = new int[g.numNodes][];
        for (int i = 0; i < g.neighbor.length; i++) {
            g.neighbor[i] = new int[nodeNeighbors.get(i).size()];
            int j = 0;
            for (Integer jobj : nodeNeighbors.get(i)) {
                g.neighbor[i][j++] = jobj;
            }
        }
        g.starts = new int[g.numNodes][];
        for (int i = 0; i < g.neighbor.length; i++) {
            g.starts[i] = new int[nodeStarts.get(i).size()];
            int j = 0;
            for (Integer jobj : nodeStarts.get(i)) {
                g.starts[i][j++] = jobj;
            }
        }
        g.length = new int[g.numNodes];
        for (int i = 0; i < g.neighbor.length; i++) {
            g.length[i] = nodeLength.get(i);
        }

        return g;
    }

    public static ArrayList<Sequence> readFastaFile(String fileName) {

        ArrayList<Sequence> sequences = new ArrayList<Sequence>();
        try {
            Scanner scanner = new Scanner(new BufferedInputStream(new FileInputStream(fileName)));
            while (scanner.hasNextLine()) {
                String nextLine = scanner.nextLine();
                if (nextLine.contains(">")) { // new sequence
                    String label = nextLine.substring(1, nextLine.length());
                    Sequence nextSeq = new Sequence();
                    nextSeq.label = label.split(" ")[0]; // get what is left of first space
                    //nextSeq.seq = "";
                    nextSeq.length = 0;
                    sequences.add(nextSeq);
                } else { // sequence
                    Sequence lastSeq = sequences.get(sequences.size() - 1);
                    lastSeq.length += nextLine.length();
                    //lastSeq.seq = lastSeq.seq + nextLine.toUpperCase();
                }
            }
        } catch (Exception ex) {
            System.err.println(ex);
            ex.printStackTrace();
            System.exit(-1);
        }

        return sequences;
    }
}
