/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package findafs;

import java.io.FileInputStream;
import java.io.BufferedInputStream;
import java.io.*;
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
        g.maxStart = 0;
        int linesRead = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                linesRead++;
                if (linesRead % 50000 == 0) {
                    System.out.println("linesRead: " + linesRead);
                }

                Scanner lineScanner;
                if (line.contains("label")) { //node
                    lineScanner = new Scanner(line);
                    lineScanner.useDelimiter(" ");
                    //System.out.println("line:" + line);
                    lineScanner.next();
                    int node = Integer.parseInt(lineScanner.next().trim());
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
                    if (node % 50000 == 0) {
                        System.out.println("reading node: " + node);
                    }
                }
                if (line.contains("->")) { //edge
                    lineScanner = new Scanner(line);
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
        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.contains(">")) { // new sequence
                    String label = line.substring(1, line.length());
                    Sequence nextSeq = new Sequence();
                    nextSeq.label = label.split(" ")[0]; // get what is left of first space
                    //nextSeq.seq = "";
                    nextSeq.length = 0;
                    sequences.add(nextSeq);
                } else { // sequence
                    Sequence lastSeq = sequences.get(sequences.size() - 1);
                    lastSeq.length += line.length();
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
