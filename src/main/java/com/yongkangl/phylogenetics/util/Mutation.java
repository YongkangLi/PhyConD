package com.yongkangl.phylogenetics.util;

import java.util.ArrayList;
import java.util.List;

public class Mutation {
    private final double time;
    private final int site;
    private final int to;

    public Mutation(double time, int site, int state) {
        this.time = time;
        this.site = site;
        this.to = state;
    }

    public double getTime() {
        return time;
    }

    public int getSite() {
        return site;
    }

    public int getTo() {
        return to;
    }

    public static void addAllMutations(List<Mutation> mutations, List<Mutation> toBeAdded) {
        for (int i = 0, j = 0; i < toBeAdded.size(); i++) {
            while (j < mutations.size() && toBeAdded.get(i).getTime() > mutations.get(j).getTime()) {
                j++;
            }
            mutations.add(j, toBeAdded.get(i));
            j++;
        }
    }
    public static List<Mutation> mergeMutations(List<Mutation> mutations1, List<Mutation> mutations2) {
        List<Mutation> mergedMutations = new ArrayList<>(mutations1);
        addAllMutations(mergedMutations, mutations2);
        return mergedMutations;
    }

    public static List<Mutation> mergePaths(List<Mutation>[] paths) {
        int n = paths.length;
        if (n ==0) {
            return new ArrayList<>();
        } else if (n == 1) {
            return paths[0];
        } else if (n == 2) {
            return mergeTwoPaths(paths[0], paths[1]);
        } else {
            ArrayList<Mutation>[] paths1 = new ArrayList[n / 2];
            ArrayList<Mutation>[] paths2 = new ArrayList[n - n / 2];
            System.arraycopy(paths, 0, paths1, 0, n / 2);
            if (n - n / 2 >= 0) {
                System.arraycopy(paths, n / 2, paths2, 0, n - n / 2);
            }
            return mergeTwoPaths(mergePaths(paths1), mergePaths(paths2));
        }
    }

    public static ArrayList<Mutation> mergeTwoPaths(List<Mutation> path1, List<Mutation> path2) {
        ArrayList<Mutation> merged = new ArrayList<>();
        int i1 = 0;
        int i2 = 0;
        while (i1 < path1.size() && i2 < path2.size()) {
            Mutation m1 = path1.get(i1);
            Mutation m2 = path2.get(i2);
            if (m1.getTime() < m2.getTime()) {
                merged.add(m1);
                i1++;
            } else {
                merged.add(m2);
                i2++;
            }
        }
        if (i1 < path1.size()) {
            for (int i = i1; i < path1.size(); i++) {
                merged.add(path1.get(i));
            }
        } else if (i2 < path2.size()) {
            for (int i = i2; i < path2.size(); i++) {
                merged.add(path2.get(i));
            }
        }
        return merged;
    }

    public static String toString(List<Mutation> mutations) {
        StringBuilder time = new StringBuilder();
        StringBuilder site = new StringBuilder();
        StringBuilder to = new StringBuilder();
        time.append("[");
        site.append("[");
        to.append("[");
        for (Mutation mutation : mutations) {
            time.append(mutation.getTime()).append(", ");
            site.append(mutation.getSite()).append(", ");
            to.append(mutation.getTo()).append(", ");
        }
        if (!mutations.isEmpty()) {
            time.delete(time.length() - 2, time.length());
            site.delete(site.length() - 2, site.length());
            to.delete(to.length() - 2, to.length());
        }
        time.append("]");
        site.append("]");
        to.append("]");
        return time + ", " + site + ", " + to;
    }
}
