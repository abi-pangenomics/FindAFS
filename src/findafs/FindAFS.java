package findafs;

/**
 *
 * @author bmumey
 */
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.PriorityQueue;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import ilog.concert.*;
import ilog.cplex.*;

public class FindAFS implements Runnable {

    static Graph g;
    static ArrayList<Sequence> sequences;
    static TreeMap<Integer, Integer> startToNode;

    static int[][] paths;
    static PriorityBlockingQueue<AFSNode> frontierQ;
    static int numThreads;
    static int nodesExplored = 0;

    // per thread variables:
    IloCplex cplex;
    int myThreadNum;

    // command line arguments:
    static int K = -1; // k-mer size
    static double eps_r = 0.0; // epsilon_r parameter
    static double mu_i = 0.0; // mu_i parameter
    static double eps_c = 0.0; // epsilon_c parameter
    static String filePrefix = ""; // file name prefix
    static int maxExploreNodes;
    static int maxSolns;
    static double minSup;

    static void readData() {
        g = ReadInput.readDotFile(filePrefix + ".dot");
        sequences = ReadInput.readFastaFile(filePrefix + ".fa");
    }

    static void buildPaths() {
        //System.out.println("g.maxStart: " + g.maxStart);
        startToNode = new TreeMap<Integer, Integer>();
        for (int i = 0; i < g.starts.length; i++) {
            for (int j = 0; j < g.starts[i].length; j++) {
                startToNode.put(g.starts[i][j], i);
            }
        }
        g.starts = null; // no longer need this.

        ArrayList<ArrayList<Integer>> pathsAL = new ArrayList<ArrayList<Integer>>();
        int curStart = 1;
        for (Sequence s : sequences) {
            s.startPos = curStart;
            int lastPos = curStart + s.length - 1;
            //System.out.println("curStart: " + curStart + " lastPos: " + lastPos);
            ArrayList path = new ArrayList<Integer>();
            while (curStart + g.length[startToNode.get(curStart)] - 1 <= lastPos) {
                path.add(startToNode.get(curStart));
                curStart += g.length[startToNode.get(curStart)] - (K - 1);
                //System.out.println(curStart);
            }
            //System.out.println("path found of length: " + path.size());
            pathsAL.add(path);
            curStart = lastPos + 2; // shift curStart to next sequence
        }
        System.out.println("number of paths: " + pathsAL.size());

        paths = new int[pathsAL.size()][];
        for (int i = 0; i < pathsAL.size(); i++) {
            ArrayList<Integer> path = pathsAL.get(i);
            paths[i] = new int[path.size()];
            for (int j = 0; j < path.size(); j++) {
                paths[i][j] = path.get(j);
            }
        }
        pathsAL.clear();

        TreeMap<Integer, TreeSet<Integer>> nodePaths = new TreeMap<Integer, TreeSet<Integer>>();
        for (int i = 0; i < g.numNodes; i++) {
            nodePaths.put(i, new TreeSet<Integer>());
        }

        for (int i = 0; i < paths.length; i++) {
            for (int j = 0; j < paths[i].length; j++) {
                nodePaths.get(paths[i][j]).add(i);
            }
        }
        g.nodePaths = new int[g.numNodes][];
        for (int i = 0; i < g.numNodes; i++) {
            g.nodePaths[i] = new int[nodePaths.get(i).size()];
            int j = 0;
            for (Integer pobj : nodePaths.get(i)) {
                g.nodePaths[i][j++] = pobj;
            }
        }

        //System.out.println("path " + nextNodeToAdd + ": " + Arrays.toString(paths[nextNodeToAdd]));
        for (int i = 0; i < g.numNodes; i++) {
            //System.out.println("node " + nextNodeToAdd + " paths: " + g.nodePaths[nextNodeToAdd]);
        }
    }

    static ArrayList<PathSegment> findMaximalSegments(AFSNode afsNode, TreeSet<Integer> pathSet) {
        ArrayList<PathSegment> segList = new ArrayList<PathSegment>();
        ArrayList<Integer> anchorPath = afsNode.getAnchorPath();

        for (Integer iobj : pathSet) {
            int[] testPath = paths[iobj];
            ArrayList<Integer> matchPos = new ArrayList<Integer>();
            for (int j = 0; j < testPath.length; j++) {
                if (afsNode.contains(testPath[j])) {
                    matchPos.add(j);
                }
            }
            //int[][] mat = new int[matchPos.size()][matchPos.size()];
            int maxK = -1;
            for (int j = 0; j < matchPos.size(); j++) {
                boolean foundNew = false;
                for (int k = Math.max(j, maxK + 1); k < matchPos.size(); k++) {
                    int numMatches = (k - j) + 1;
                    int segLength = (matchPos.get(k) - matchPos.get(j)) + 1;
                    if (numMatches >= (1.0 - eps_r) * anchorPath.size() && // check eps_r and mu_i constraints
                            (segLength - numMatches) <= mu_i) {
                        maxK = k;
                        foundNew = true;
                    }
                }
                if (foundNew) {
                    PathSegment ps = new PathSegment();
                    ps.path = iobj;
                    ps.start = matchPos.get(j);
                    ps.stop = matchPos.get(maxK);
                    ps.support = 1.0;
                    segList.add(ps);
                }
            }
        }
        return segList;
    }

    void epsCFilter(AFSNode node, ArrayList<PathSegment> supportingSegments) {

        // create an ILP
        try {
            IloNumVar[] psv = new IloNumVar[supportingSegments.size()];
            for (int i = 0; i < supportingSegments.size(); i++) {
                psv[i] = cplex.numVar(0, 1.0, "ps(" + i + ")");
            }
            ArrayList<Integer> pa = node.getAnchorPath();

            // overlap constraints:
            int maxJ = -1;
            for (int i = 0; i < supportingSegments.size(); i++) {
                boolean foundNew = false;
                for (int j = Math.max(i, maxJ + 1); j < supportingSegments.size(); j++) {
                    if (supportingSegments.get(i).path == supportingSegments.get(j).path
                            && supportingSegments.get(i).stop >= supportingSegments.get(j).start) {
                        maxJ = j;
                        foundNew = true;
                    }
                }
                if (foundNew) {
                    IloLinearNumExpr expr = cplex.linearNumExpr(0.0);
                    for (int k = i; k <= maxJ; k++) {
                        expr.addTerm(1.0, psv[k]);
                    }
                    cplex.addLe(expr, 1.0);
                }
            }

            // eps_c constraints:
            IloLinearNumExpr[] expr = new IloLinearNumExpr[pa.size()];
            for (int i = 0; i < pa.size(); i++) {
                expr[i] = cplex.linearNumExpr(0.0);
            }
            for (int i = 0; i < supportingSegments.size(); i++) {
                PathSegment ps = supportingSegments.get(i);
                for (int j = ps.start; j <= ps.stop; j++) {
                    int index = pa.indexOf(paths[ps.path][j]);
                    if (index >= 0) {
                        expr[index].addTerm(1.0, psv[i]);
                    }
                }
            }
            for (int i = 0; i < pa.size(); i++) {
                for (int j = 0; j < supportingSegments.size(); j++) {
                    expr[i].addTerm(-(1.0 - eps_c), psv[j]);
                }
                cplex.addGe(expr[i], 0.0);
            }

            // objective:
            IloLinearNumExpr maxExpr = cplex.linearNumExpr();
            for (int i = 0; i < supportingSegments.size(); i++) {
                maxExpr.addTerm(1.0, psv[i]);
            }

            // solveLP the problem
            IloObjective obj = cplex.maximize(maxExpr);
            cplex.add(obj);
            cplex.setOut(null);
            //cplex.exportModel("scd.lp");
            if (cplex.solve()) {
                for (int i = 0; i < supportingSegments.size(); i++) {
                    supportingSegments.get(i).support = cplex.getValue(psv[i]);
                }
            }
            cplex.clearModel();
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }
    }

    ArrayList<PathSegment> computeSupport(AFSNode afsNode) {
        TreeSet<Integer> pathSet = new TreeSet<Integer>();
        for (int i = 0; i < g.nodePaths[afsNode.node].length; i++) {
            pathSet.add(g.nodePaths[afsNode.node][i]);
        }
        if (afsNode.parent != null) {
            for (int p : afsNode.parent.supportingPaths) {
                pathSet.add(p);
            }
        }
        ArrayList<PathSegment> supportingSegments = findMaximalSegments(afsNode, pathSet);
        epsCFilter(afsNode, supportingSegments);

        afsNode.support = (float) 0.0;

        TreeSet<Integer> supPaths = new TreeSet<Integer>();
        for (PathSegment ps : supportingSegments) {
            afsNode.support += ps.support;
            supPaths.add(ps.path);
        }

        afsNode.supportingPaths = new int[supPaths.size()];
        int i = 0;
        for (Integer iobj : supPaths) {
            afsNode.supportingPaths[i++] = iobj;
        }

        return supportingSegments;
    }

    void expand(AFSNode parentNode) {
        for (int i = 0; i < g.neighbor[parentNode.node].length; i++) {
            int neighbor = g.neighbor[parentNode.node][i];
            if (!parentNode.contains(neighbor)) { // only consider simple paths
                AFSNode neighborAFSNode = new AFSNode();
                neighborAFSNode.node = neighbor;
                neighborAFSNode.parent = parentNode;
                computeSupport(neighborAFSNode);

                frontierQ.add(neighborAFSNode);
                //bestQ.add(neighborAFSNode);
            }
        }
    }

    void exploreSolns() { // executed by each thread:
        for (int nextNodeToAdd = 0; nextNodeToAdd < g.numNodes; nextNodeToAdd++) {
            if (nextNodeToAdd % numThreads == myThreadNum) {
                AFSNode newNode = new AFSNode();
                newNode.parent = null;
                newNode.node = nextNodeToAdd;
                computeSupport(newNode);

                //bestQ.add(newNode);
                frontierQ.add(newNode);
                if (nextNodeToAdd % 10000 == 0) {
                    System.out.println("nodes added: " + nextNodeToAdd);
                }
            }
        }

        AFSNode top = frontierQ.poll();
        while (top != null && nodesExplored < maxExploreNodes) {
            expand(top);
            nodesExplored++;
            if (nodesExplored % 10000 == 0) {
                System.out.println("nodes explored: " + nodesExplored);
            }
            top = frontierQ.poll();
        }

    }

    void findBestSolns() {
        PriorityQueue<AFSNode> bestQ = new PriorityQueue<AFSNode>();
        for (AFSNode leafNode : frontierQ) {
            AFSNode curNode = leafNode;
            while (curNode != null && curNode.support < minSup) {
                curNode = curNode.parent;
            }
            if (leafNode != null && leafNode.support >= minSup) {
                bestQ.add(leafNode);
            }
        }
        for (AFSNode bNode : bestQ) {
            AFSNode curNode = bNode.parent;
            while (curNode != null) {
                bestQ.remove(curNode);
                curNode = curNode.parent;
            }
        }

        writeBED(bestQ);
    }

    static int[] findFastaLoc(PathSegment ps) {
        int[] startStop = new int[2];
        int curPos = sequences.get(ps.path).startPos;

        int curIndex = 0;
        while (curIndex != ps.start) {
            curPos += g.length[startToNode.get(curPos)] - (K - 1);
            curIndex++;
        }
        startStop[0] = curPos - sequences.get(ps.path).startPos; // assume fasta seq indices start at 0
        while (curIndex != ps.stop) {
            curPos += g.length[startToNode.get(curPos)] - (K - 1);
            curIndex++;
        }
        startStop[1] = curPos + g.length[startToNode.get(curPos)] - 1 - sequences.get(ps.path).startPos + 1; // last position is excluded in BED format
        return startStop;
    }

    static void printAFS(AFSNode afsNode, ArrayList<PathSegment> supportingSegments) {
        System.out.println("anchor path: " + afsNode.getAnchorPath());
        System.out.println("total support: " + afsNode.support);

        for (PathSegment ps : supportingSegments) {
            System.out.println("from fasta seq: " + sequences.get(ps.path).label);
            System.out.print("subpath: [");
            for (int i = ps.start; i <= ps.stop; i++) {
                if (i > ps.start) {
                    System.out.print(",");
                }
                System.out.print(paths[ps.path][i]);
            }
            System.out.println("] support:" + ps.support);
            System.out.print("fasta location: ");
            int[] startStop = findFastaLoc(ps);
            System.out.println(startStop[0] + "," + startStop[1]);
        }
        System.out.println();
    }

    void writeBED(PriorityQueue<AFSNode> bestQ) {
        String[] colors = {"122,39,25", "92,227,60", "225,70,233", "100,198,222", "232,176,49", "50,39,85", "67,101,33", "222,142,186", "92,119,227", "206,225,151", "227,44,118", "229,66,41", "47,36,24", "225,167,130", "120,132,131", "104,232,178", "158,43,133", "228,228,42", "213,217,213", "118,64,79", "88,155,219", "226,118,222", "146,197,53", "222,100,89", "224,117,41", "160,96,228", "137,89,151", "126,209,119", "145,109,70", "91,176,164", "54,81,103", "164,174,137", "172,166,48", "56,86,143", "210,184,226", "175,123,35", "129,161,88", "158,47,85", "87,231,225", "216,189,112", "49,111,75", "89,137,168", "209,118,134", "33,63,44", "166,128,142", "53,137,55", "80,76,161", "170,124,221", "57,62,13", "176,40,40", "94,179,129", "71,176,51", "223,62,170", "78,25,30", "148,69,172", "122,105,31", "56,33,53", "112,150,40", "239,111,176", "96,55,25", "107,90,87", "164,74,28", "171,198,226", "152,131,176", "166,225,211", "53,121,117", "220,58,86", "86,18,56", "225,197,171", "139,142,217", "216,151,223", "97,229,117", "225,155,85", "31,48,58", "160,146,88", "185,71,129", "164,233,55", "234,171,187", "110,97,125", "177,169,175", "177,104,68", "97,48,122", "237,139,128", "187,96,166", "225,90,127", "97,92,55", "124,35,99", "210,64,194", "154,88,84", "100,63,100", "140,42,54", "105,132,99", "186,227,103", "224,222,81", "191,140,126", "200,230,182", "166,87,123", "72,74,58", "212,222,124", "205,52,136"};
        try {
            String paramString = "-k" + K + "-r" + eps_r + "-i" + mu_i + "-c" + eps_c + "-mn" + maxExploreNodes + "-ms" + maxSolns;
            BufferedWriter out = new BufferedWriter(new FileWriter(filePrefix + paramString + ".bed"));
            int count = 0;
            while (!bestQ.isEmpty() && count < maxSolns) {
                AFSNode top = bestQ.poll();
                ArrayList<PathSegment> supportingSegments = computeSupport(top);
                printAFS(top, supportingSegments);
                for (PathSegment ps : supportingSegments) {
                    if (ps.support > 0.0) {
                        String name = sequences.get(ps.path).label;
                        int[] startStop = findFastaLoc(ps);
                        out.write(
                                name // chrom
                                + "\t" + startStop[0] // chromStart (starts with 0)
                                + "\t" + startStop[1] // chromEnd
                                + "\t" + "afs-" + count // name
                                + "\t" + Math.round(top.support) // score
                                + "\t+" // strand
                                + "\t" + 0 // thickstart
                                + "\t" + 0 // thickend
                                + "\t" + colors[count % colors.length] // itemRGB
                                + "\n");
                    }
                }
                count++;
            }
            out.close();
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }

    }

    public void run() {
        System.out.println("starting thread: " + myThreadNum);
        try {
            cplex = new IloCplex();
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }
        exploreSolns();

        // write results:
        if (myThreadNum == 0) {
            findBestSolns();
        }
    }

    public static void main(String[] args) {

        // parse args:
        if (args.length != 8) {
            System.out.println("Usage: java findAFS K eps_r mu_i eps_c minSup filePrefix maxExploreNodes maxSolns");
            System.out.println(Arrays.toString(args));
            System.exit(0);
        }
        K = Integer.parseInt(args[0]);
        eps_r = Double.parseDouble(args[1]);
        mu_i = Double.parseDouble(args[2]);
        eps_c = Double.parseDouble(args[3]);
        minSup = Double.parseDouble(args[4]);
        filePrefix = args[5];
        maxExploreNodes = Integer.parseInt(args[6]);
        maxSolns = Integer.parseInt(args[7]);

        readData();
        buildPaths();

        frontierQ = new PriorityBlockingQueue<AFSNode>();

        // launch threads:
        numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] worker = new Thread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            FindAFS next = new FindAFS();
            worker[i] = new Thread(next);
            next.myThreadNum = i;
            worker[i].start();
        }
        try {
            for (int i = 0; i < numThreads; i++) {
                worker[i].join();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            System.exit(-1);
        }
    }
}
