/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package findafs;

/**
 *
 * @author bmumey
 */
import java.util.Arrays;
import java.util.ArrayList;
import java.util.TreeSet;

public class AFSNode implements Comparable<AFSNode> {

    int node;
    AFSNode parent;
    int[] supportingPaths;
    //ArrayList<PathSegment> supportingSegments;
    float support;
    AFSNode suffixlink;
    AFSNode[] child;
    

    public int compareTo(AFSNode other) {
        return Double.compare(other.support, support);
    }

    public boolean pathContains(int n) {
        if (node == n) {
            return true;
        } else if (parent == null) {
            return false;
        } else {
            return parent.pathContains(n);
        }
    }

    int depth() {
        if (parent == null) {
            return 1;
        } else {
            return 1 + parent.depth();
        }
    }
    
    int[] getAnchorPath() {
        ArrayList<Integer> ap = new ArrayList<Integer>();
        AFSNode curNode = this;
        while (curNode != null) {
            ap.add(node);
            curNode = curNode.parent;
        }
        int[] pa = new int[ap.size()];
        int i = ap.size() - 1;
        for (Integer N : ap) {
            pa[i--] = N;
        }
        return pa;
    }
}
