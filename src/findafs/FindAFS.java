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
import java.util.concurrent.ConcurrentSkipListSet;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.PriorityQueue;
import java.util.LinkedList;
import java.util.HashMap;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import ilog.concert.*;
import ilog.cplex.*;

public class FindAFS implements Runnable {

    static Graph g;
    static ArrayList<Sequence> sequences;
    static char[] fastaConcat;
    static int fastaConcatLen;
    static TreeMap<Integer, Integer> startToNode;
    static int[][] paths;

    static PriorityBlockingQueue<AFSNode> frontierQ;
    //static PriorityBlockingQueue<AFSNode> bestQ;
    static ConcurrentSkipListSet<AFSNode> afsQ;
    static LinkedBlockingQueue<AFSNode> suffixQ;
    static AFSNode[] startNode;
    static int numThreads;
    static int nodesExpanded = 0;
    static int epsCNodes = 0;

    // per thread variables:
    static Thread[] worker;
    IloCplex cplex;
    int myThreadNum;

    // command line arguments:
    static int K = -1; // j-mer size
    static double eps_r = 0.0; // epsilon_r parameter
    static double mu_i = 0.0; // mu_i parameter
    static double eps_c = 0.0; // epsilon_c parameter
    static String filePrefix = ""; // file name prefix
    static int maxExploreNodes;
    static int maxSolns;
    static double minSup;
    static int maxSegPPath;

    static void readData() {
        g = ReadInput.readDotFile(filePrefix + ".dot");

        startToNode = new TreeMap<Integer, Integer>();
        for (int i = 0; i < g.starts.length; i++) {
            for (int j = 0; j < g.starts[i].length; j++) {
                startToNode.put(g.starts[i][j], i);
            }
            int firstStart = g.starts[i][0];
            g.starts[i] = new int[1]; // only save 1st start
            g.starts[i][0] = firstStart;
        }

        sequences = ReadInput.readFastaFile(filePrefix + ".fa");
    }

    static int[] removeDuplicates(int[] a) {
        int[] b = null;

        if (a.length > 0) {
            ArrayList<Integer> B = new ArrayList<Integer>();
            B.add(a[0]);
            for (int i = 1; i < a.length; i++) {
                if (a[i] != B.get(B.size() - 1)) {
                    B.add(a[i]);
                }
            }
            int i = 0;
            b = new int[B.size()];
            for (Integer pobj : B) {
                b[i++] = pobj;
            }
        }
        return b;
    }

    static void buildPaths() {
        ArrayList<ArrayList<Integer>> pathsAL = new ArrayList<ArrayList<Integer>>();
        int curStart = 1;
        int seqStart = 1;
        int seqEnd;
        int prevStart = 0;
        for (Sequence s : sequences) {
            s.startPos = seqStart;
            s.length = s.seq.length();
            seqEnd = seqStart + s.length - 1;
            curStart = seqStart;
            ArrayList path = new ArrayList<Integer>();
            while (curStart + g.length[startToNode.get(curStart)] - 1 < seqEnd) {
                path.add(startToNode.get(curStart));
                curStart += g.length[startToNode.get(curStart)] - (K - 1);
            }
            pathsAL.add(path);
            seqStart = seqEnd + 2;

            fastaConcatLen += 1 + s.seq.length();
        }
        fastaConcat = new char[fastaConcatLen];
        for (Sequence s : sequences) {
            fastaConcat[s.startPos - 1] = '$';
            for (int i = 0; i < s.length; i++) {
                fastaConcat[s.startPos + i] = s.seq.charAt(i);
            }
            s.seq = null; // no longer needed
        }
        //System.out.println(Arrays.toString(fastaConcat));
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
        pathsAL = null; // can be gc'ed

        System.out.println("finding node paths");

        // find paths for each node:
        int[] len = new int[g.numNodes];
        for (int i = 0; i < g.numNodes; i++) {
            len[i] = 0;
        }
        for (int i = 0; i < paths.length; i++) {
            for (int j = 0; j < paths[i].length; j++) {
                len[paths[i][j]]++;
            }
        }
        g.nodePaths = new int[g.numNodes][];
        for (int i = 0; i < g.numNodes; i++) {
            g.nodePaths[i] = new int[len[i]];
            len[i] = 0;
        }

        for (int i = 0; i < paths.length; i++) {
            for (int j = 0; j < paths[i].length; j++) {
                g.nodePaths[paths[i][j]][len[paths[i][j]]++] = i;
            }
        }
        for (int i = 0; i < g.numNodes; i++) {
            if (g.nodePaths[i] != null) {
                g.nodePaths[i] = removeDuplicates(g.nodePaths[i]);
            }
        }
    }

    static int[][] lcs(int[] pa, int[] seg) { // needs to be updated to only use matching nodes..
        int[][] m = new int[2][];
        m[0] = new int[pa.length]; // m[0][i] == 1 iff pa[i] in LCS
        m[1] = new int[seg.length]; // m[1][j] == 1 iff seg[j] in LCS
        int[][] len = new int[pa.length + 1][seg.length + 1];
        int[][] best = new int[pa.length + 1][seg.length + 1];

        for (int i = 1; i <= pa.length; i++) {
            for (int j = 1; j <= seg.length; j++) {
                if (pa[i - 1] == seg[j - 1]) {
                    len[i][j] = 1 + len[i - 1][j - 1];
                    best[i][j] = 1; // match, go up and left
                } else if (len[i][j - 1] >= len[i - 1][j]) {
                    len[i][j] = len[i][j - 1];
                    best[i][j] = 2; // mismatch, go left       
                } else {
                    len[i][j] = len[i - 1][j];
                    best[i][j] = 3; // mismatch, go up                   
                }
            }
        }
        int curI = pa.length;
        int curJ = seg.length;
        while (curI > 0 && curJ > 0) {
            if (best[curI][curJ] == 1) {
                m[0][curI - 1] = 1;
                m[1][curJ - 1] = 1;
                curI--;
                curJ--;
            } else if (best[curI][curJ] == 2) {
                curJ--;
            } else {
                curI--;
            }
        }
        return m;
    }

    static boolean[][] set(int[] pa, int[] path, ArrayList<Integer> matchPos, int start, int stop) {
        boolean[][] m = new boolean[2][];
        m[0] = new boolean[pa.length]; // m[0][i] == 1 iff pa[i] in LCS
        m[1] = new boolean[stop - start + 1]; // m[1][i] == 1 iff path[matchPos[i]] matched, start <= i <= stop

        HashMap<Integer, LinkedList<Integer>> hs = new HashMap<Integer, LinkedList<Integer>>();
        for (int i = 0; i < pa.length; i++) {
            if (!hs.containsKey(pa[i])) {
                hs.put(pa[i], new LinkedList<Integer>());
                hs.get(pa[i]).add(i);
            }
        }
        for (int i = start; i <= stop; i++) {
            if (hs.containsKey(path[matchPos.get(i)]) && !hs.get(path[matchPos.get(i)]).isEmpty()) {
                m[1][i - start] = true;
                m[0][hs.get(path[matchPos.get(i)]).pop()] = true;
            }
        }
        return m;
    }

    static int matchLength(int[] pa, boolean[] matched) {
        int len = 0;
        int lastMatched = 0;
        int curStart = 1;
        int curEnd;
        for (int i = 0; i < pa.length; i++) {
            curEnd = curStart + g.length[pa[i]] - 1;
            if (matched[i]) {
                len += g.length[pa[i]] - Math.max(lastMatched - curStart + 1, 0);
                lastMatched = curEnd;
            }
            curStart += g.length[pa[i]] - (K - 1);
        }
        return len;
    }

    static int matchLength(int[] path, ArrayList<Integer> matchPos, int start, int stop, boolean[] matched) {
        int len = 0;
        int lastMatched = 0;
        int curStart = 1;
        int curEnd;
        int curMatchPos = 0;
        for (int i = matchPos.get(start); i <= matchPos.get(stop); i++) {
            curEnd = curStart + g.length[path[i]] - 1;
            if (i == matchPos.get(start + curMatchPos)) {
                if (matched[curMatchPos]) {
                    len += g.length[path[i]] - Math.max(lastMatched - curStart + 1, 0);
                    lastMatched = curEnd;
                }
                curMatchPos++;
            }
            curStart += g.length[path[i]] - (K - 1);
        }
        return len;
    }

    static int[] comparePaths(int[] pa, int[] path, ArrayList<Integer> matchPos, int start, int stop) {
        int[] matchLen = new int[2]; //matchLen[0] pa matchLen, matchLen[1] seg matchLen
//        booolean[][] match = lcs(pa, seg);
        boolean[][] match = set(pa, path, matchPos, start, stop);
        matchLen[0] = matchLength(pa, match[0]);
        matchLen[1] = matchLength(path, matchPos, start, stop, match[1]);
        return matchLen;
    }

    static ArrayList<PathSegment> findMaximalSegments(AFSNode afsNode, TreeSet<Integer> pathSet) {
        ArrayList<PathSegment> segList = new ArrayList<PathSegment>();
        int[] anchorPath = afsNode.getAnchorPath();
        int palength = (K - 1);
        for (int i : anchorPath) {
            palength += g.length[i] - (K - 1);
        }
        for (Integer P : pathSet) {
            int[] testPath = paths[P];
            ArrayList<Integer> matchPos = new ArrayList<Integer>();
            for (int j = 0; j < testPath.length; j++) {
                if (afsNode.pathContains(testPath[j])) {
                    matchPos.add(j);
                }
            }
            if (matchPos.isEmpty()) {
                break;
            }
            int maxJ = -1;
            boolean foundNew;
            int curStart = 1;

            int foundCount = 0;
            for (int i = 0; i < matchPos.size(); i++) {
                foundNew = false;
                int prevJ = i;
                int curEnd = curStart + g.length[testPath[matchPos.get(i)]] - 1;
                for (int j = Math.max(i, maxJ + 1); j < matchPos.size(); j++) {

                    for (int k = matchPos.get(prevJ) + 1; k <= matchPos.get(j); k++) {
                        curEnd += g.length[testPath[k]] - (K - 1);
                    }
                    int segLength = (curEnd - curStart + 1);
                    prevJ = j;
                    int[] mlength = comparePaths(anchorPath, testPath, matchPos, i, j);
                    // mlength[0] = pa match length, malength[1] = seg match length
                    if (mlength[0] >= (1.0 - eps_r) * palength && // check eps_r and mu_i constraints
                            (segLength - mlength[1]) <= mu_i * palength) {
                        maxJ = j;
                        foundNew = true;
                    } else if (foundNew) { // stop moving j is a match has been previously found
                        break;
                    }
                }
                if (foundNew) {
                    PathSegment ps = new PathSegment();
                    ps.path = P;
                    ps.start = matchPos.get(i);
                    ps.stop = matchPos.get(maxJ);
                    ps.support = 0.0;
                    segList.add(ps);
                    foundCount++;
                    if (maxSegPPath > 0 && foundCount >= maxSegPPath) {
                        break;
                    }
                }

                if (i < matchPos.size() - 1) {
                    for (int k = matchPos.get(i); k < matchPos.get(i + 1); k++) {
                        curStart += g.length[testPath[k]] - (K - 1);
                    }
                }
            }
        }
        return segList;
    }

    static int getIndex(int[] a, int x) {
        for (int i = 0; i < a.length; i++) {
            if (a[i] == x) {
                return i;
            }
        }
        return -1;
    }

    void epsCFilter(AFSNode node, ArrayList<PathSegment> supportingSegments) {

        // create an ILP
        try {
            IloNumVar[] psv = new IloNumVar[supportingSegments.size()];
            for (int i = 0; i < supportingSegments.size(); i++) {
                psv[i] = cplex.numVar(0, 1.0, "ps(" + i + ")");
            }
            int[] pa = node.getAnchorPath();

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
            IloLinearNumExpr[] expr = new IloLinearNumExpr[pa.length];
            for (int i = 0; i < pa.length; i++) {
                expr[i] = cplex.linearNumExpr(0.0);
            }
            for (int i = 0; i < supportingSegments.size(); i++) {
                PathSegment ps = supportingSegments.get(i);
                for (int j = ps.start; j <= ps.stop; j++) {
                    int index = getIndex(pa, paths[ps.path][j]);
                    if (index >= 0) {
                        expr[index].addTerm(1.0, psv[i]);
                    }
                }
            }
            for (int i = 0; i < pa.length; i++) {
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
        for (PathSegment ps : supportingSegments) {
            afsNode.support += ps.support;
        }
        TreeSet<Integer> supPaths = new TreeSet<Integer>();
        for (PathSegment ps : supportingSegments) {
            if (ps.support > 0) {
                supPaths.add(ps.path);
            }
        }
        afsNode.supportingPaths = new int[supPaths.size()];
        int i = 0;
        for (Integer I : supPaths) {
            afsNode.supportingPaths[i++] = I;
        }
        return supportingSegments;
    }

    void expand(AFSNode parentNode) {
        ArrayList<AFSNode> children = new ArrayList<AFSNode>();
        for (int i = 0; i < g.neighbor[parentNode.node].length; i++) {
            int neighbor = g.neighbor[parentNode.node][i];
            if (g.nodePaths[neighbor] != null
                    && g.nodePaths[neighbor].length >= (1.0 - eps_c) * minSup) {
                // && !parentNode.pathContains(neighbor)) {
                AFSNode newNode = new AFSNode();
                newNode.node = neighbor;
                newNode.parent = parentNode;
                newNode.child = null;
                children.add(newNode);
                computeSupport(newNode);
                frontierQ.add(newNode);
                //bestQ.add(newNode);
            }
        }
        parentNode.child = new AFSNode[children.size()];
        int i = 0;
        for (AFSNode n : children) {
            parentNode.child[i++] = n;
        }
    }

    void exploreSolns() { // executed by each thread:
        for (int nextNodeToAdd = 0; nextNodeToAdd < g.numNodes; nextNodeToAdd++) {
            if (nextNodeToAdd % numThreads == myThreadNum
                    && g.nodePaths[nextNodeToAdd] != null
                    && g.nodePaths[nextNodeToAdd].length >= (1.0 - eps_c) * minSup) {
                AFSNode newNode = new AFSNode();
                newNode.parent = null;
                newNode.node = nextNodeToAdd;
                computeSupport(newNode);
                if (newNode.support >= (1.0 - eps_c) * minSup) {
                    frontierQ.add(newNode);
                    //bestQ.add(newNode);
                    startNode[newNode.node] = newNode;
                    suffixQ.add(newNode);
                    if (nextNodeToAdd % 20000 == 0) {
                        System.out.println("nodes added: " + nextNodeToAdd);
                    }
                }
            }
        }

        AFSNode top;
        while ((top = frontierQ.poll()) != null && nodesExpanded < maxExploreNodes) {
            expand(top);
            nodesExpanded++;
            if (nodesExpanded % 100 == 0) {
                System.out.println("nodes expanded: " + nodesExpanded);
            }
        }
//        while ((top = bestQ.poll()) != null && epsCNodes < 10 * maxSolns) {
//            computeSupport(top, true);
//            epsCNodes++;
//        }
    }

    void findSuffixLinks() {
        // at this point suffixQ contains all level 1 nodes
        System.out.println("finding suffix links");
        AFSNode front;
        while ((front = suffixQ.poll()) != null) {
            if (front.parent == null) {
                front.suffixlink = null;
            } else {
                AFSNode parent = front.parent;
                AFSNode sparent = parent.suffixlink;
                while (front.suffixlink == null) {
                    if (sparent != null && sparent.child != null) {
                        for (AFSNode n : sparent.child) {
                            if (n.node == front.node) {
                                front.suffixlink = n;
                                break;
                            }
                        }
                    }
                    if (sparent == null) {
                        front.suffixlink = startNode[front.node];
                    } else if (front.suffixlink == null) {
                        sparent = sparent.suffixlink;
                    }
                }
            }
            // add front's children to suffixQ:
            if (front.child != null) {
                for (AFSNode n : front.child) {
                    suffixQ.add(n);
                }
            }
        }
    }

    double maxSupport(AFSNode n) {
        if (n.child == null) {
            return n.support;
        }
        double maxChildSup = 0.0;
        for (int i = 0; i < n.child.length; i++) {
            maxChildSup = Math.max(maxChildSup, maxSupport(n.child[i]));
        }
        if (n.support > maxChildSup && n.support >= minSup) {
            afsQ.add(n);
        }
        return Math.max(n.support, maxChildSup);
    }

    void findBestSolns() {
        if (myThreadNum == 0) {
            try {
                for (int i = 1; i < numThreads; i++) { // wait for other threads to finish
                    worker[i].join();
                }
            } catch (Exception ex) {
                ex.printStackTrace();
                System.exit(-1);
            }

            System.out.println("single thread finishing up..");
            for (int i = 0; i < g.numNodes; i++) {
                if (startNode[i] != null) {
                    maxSupport(startNode[i]);
                }
            }

            findSuffixLinks();
            TreeSet<AFSNode> removeSet = new TreeSet<AFSNode>();
            for (AFSNode bNode : afsQ) {
                AFSNode curNode = bNode;
                while (curNode != null) {
                    AFSNode sufNode = curNode.suffixlink;
                    while (sufNode != null) {
                        if (sufNode.support <= bNode.support) {
                            removeSet.add(sufNode);
                        }
                        sufNode = sufNode.suffixlink;
                    }
                    curNode = curNode.parent;
                }
            }
            System.out.println("removeSet size: " + removeSet.size());
            for (AFSNode rmNode : removeSet) {
                afsQ.remove(rmNode);
            }

            // write results:
            outputAFSs();
        }
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
        System.out.println("anchor path: " + Arrays.toString(afsNode.getAnchorPath()));
        System.out.println("total support: " + afsNode.support);

        for (PathSegment ps : supportingSegments) {
            System.out.println("from fasta seq: " + sequences.get(ps.path).label);
            System.out.println("support:" + ps.support);
            System.out.println("length (nodes): " + (ps.stop - ps.start + 1));
            System.out.print("fasta location: ");
            int[] startStop = findFastaLoc(ps);
            System.out.println(startStop[0] + "," + startStop[1]);
        }
        System.out.println();
    }

    int findLen(int[] pa) {
        int c = 0;
        for (int i = 0; i < pa.length; i++) {
            if (i == 0) {
                c += g.length[pa[i]];

            } else {
                c += g.length[pa[i]] - (K - 1);
            }
        }
        return c;
    }

    void outputAFSs() {
        String[] colors = {"122,39,25", "92,227,60", "225,70,233", "100,198,222", "232,176,49", "50,39,85", "67,101,33", "222,142,186", "92,119,227", "206,225,151", "227,44,118", "229,66,41", "47,36,24", "225,167,130", "120,132,131", "104,232,178", "158,43,133", "228,228,42", "213,217,213", "118,64,79", "88,155,219", "226,118,222", "146,197,53", "222,100,89", "224,117,41", "160,96,228", "137,89,151", "126,209,119", "145,109,70", "91,176,164", "54,81,103", "164,174,137", "172,166,48", "56,86,143", "210,184,226", "175,123,35", "129,161,88", "158,47,85", "87,231,225", "216,189,112", "49,111,75", "89,137,168", "209,118,134", "33,63,44", "166,128,142", "53,137,55", "80,76,161", "170,124,221", "57,62,13", "176,40,40", "94,179,129", "71,176,51", "223,62,170", "78,25,30", "148,69,172", "122,105,31", "56,33,53", "112,150,40", "239,111,176", "96,55,25", "107,90,87", "164,74,28", "171,198,226", "152,131,176", "166,225,211", "53,121,117", "220,58,86", "86,18,56", "225,197,171", "139,142,217", "216,151,223", "97,229,117", "225,155,85", "31,48,58", "160,146,88", "185,71,129", "164,233,55", "234,171,187", "110,97,125", "177,169,175", "177,104,68", "97,48,122", "237,139,128", "187,96,166", "225,90,127", "97,92,55", "124,35,99", "210,64,194", "154,88,84", "100,63,100", "140,42,54", "105,132,99", "186,227,103", "224,222,81", "191,140,126", "200,230,182", "166,87,123", "72,74,58", "212,222,124", "205,52,136"};
        try {
            String paramString = "-k" + K + "-r" + eps_r + "-i" + mu_i + "-c" + eps_c + "-sp" + minSup + "-xp" + maxExploreNodes + "-ms" + maxSolns + "-msp" + maxSegPPath;
            BufferedWriter bedOut = new BufferedWriter(new FileWriter(filePrefix + paramString + ".bed"));
            int count = 0;
            AFSNode top;
            while ((top = afsQ.pollFirst()) != null && count < maxSolns) {
                ArrayList<PathSegment> supportingSegments = computeSupport(top);
//                printAFS(top, supportingSegments);
                int[] pa = top.getAnchorPath();
                BufferedWriter afsOut = new BufferedWriter(new FileWriter(filePrefix + paramString + "-afs-" + count + ".fa"));
                afsOut.write(">afs-" + count + " support: " + top.support + " length: " + findLen(pa) + "\n");
                int c = 0;
                for (int i = 0; i < pa.length; i++) {
                    int start = g.starts[pa[i]][0];
                    if (i == 0) {
                        for (int j = 0; j < g.length[pa[i]]; j++) {
                            afsOut.write(fastaConcat[start + j]);
                            if (++c % 80 == 0) {
                                afsOut.write("\n");
                            }
                        }
                    } else {
                        for (int j = (K - 1); j < g.length[pa[i]]; j++) {
                            afsOut.write(fastaConcat[start + j]);
                            if (++c % 80 == 0) {
                                afsOut.write("\n");
                            }
                        }
                    }
                }
                if (c % 80 != 0) {
                    afsOut.write("\n");
                }

                for (PathSegment ps : supportingSegments) {
                    if (ps.support > 0.0) {
                        String name = sequences.get(ps.path).label;
                        int[] startStop = findFastaLoc(ps);
                        bedOut.write(
                                name // chrom
                                + "\t" + startStop[0] // chromStart (starts with 0)
                                + "\t" + startStop[1] // chromEnd
                                + "\t" + "afs-" + count // name
                                + "\t" + Math.round(top.support) // score
                                + "\t+" // strand
                                + "\t" + 0 // thickstart
                                + "\t" + 0 // thickend
                                + "\t" + colors[count % colors.length] // itemRGB
                                + "\t" + (startStop[1] - startStop[0]) // AFS length
                                + "\n");

                        afsOut.write(">" + name + " segment for: afs-" + count + " support: " + ps.support + " length: " + (startStop[1] - startStop[0]) + "\n");
                        c = 0;
                        for (int i = ps.start; i <= ps.stop; i++) {
                            int start = g.starts[paths[ps.path][i]][0];
                            if (i == ps.start) {
                                for (int j = 0; j < g.length[paths[ps.path][i]]; j++) {
                                    afsOut.write(fastaConcat[start + j]);
                                    if (++c % 80 == 0) {
                                        afsOut.write("\n");
                                    }
                                }
                            } else {
                                for (int j = (K - 1); j < g.length[paths[ps.path][i]]; j++) {
                                    afsOut.write(fastaConcat[start + j]);
                                    if (++c % 80 == 0) {
                                        afsOut.write("\n");
                                    }
                                }
                            }
                        }
                        if (c % 80 != 0) {
                            afsOut.write("\n");
                        }
                    }
                }
                count++;
                afsOut.close();
            }
            bedOut.close();

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
        findBestSolns();
    }

    public static void main(String[] args) {
        // parse args:
        if (args.length != 9) {
            System.out.println("Usage: java findAFS K eps_r mu_i eps_c minSup filePrefix maxExploreNodes maxSolns maxSegPPath");
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
        maxSegPPath = Integer.parseInt(args[8]);

        readData();
        buildPaths();

        frontierQ = new PriorityBlockingQueue<AFSNode>();
        afsQ = new ConcurrentSkipListSet<AFSNode>();
        startNode = new AFSNode[g.numNodes];
        suffixQ = new LinkedBlockingQueue<AFSNode>();

        // launch threads:
        numThreads = Runtime.getRuntime().availableProcessors();
        worker = new Thread[numThreads];
        for (int i = 0; i < numThreads; i++) {
            FindAFS next = new FindAFS();
            worker[i] = new Thread(next);
            next.myThreadNum = i;
            worker[i].start();
        }
    }
}
