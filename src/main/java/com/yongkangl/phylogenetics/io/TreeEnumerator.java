package com.yongkangl.phylogenetics.io;

import java.util.*;

public class TreeEnumerator {
    private final int leafCount;

    public TreeEnumerator(int leafCount) {
        if (leafCount <= 0) {
            throw new IllegalArgumentException("Leaf count must be positive.");
        }
        this.leafCount = leafCount;
    }

    /**
     * Generates all unique, labeled, rooted binary trees for the given number of leaves.
     * @return A list of root nodes for all possible tree topologies.
     */
    public List<TreeNode> enumerate() {
        if (this.leafCount == 1) {
            TreeNode leaf = new TreeNode(null);
            leaf.setTipName("tip_1");
            return Collections.singletonList(leaf);
        }

        List<String> leafLabels = new ArrayList<>();
        for (int i = 1; i <= leafCount; i++) {
            leafLabels.add("tip_" + i);
        }
        return buildBinaryTrees(leafLabels);
    }

    private List<TreeNode> buildBinaryTrees(List<String> leaves) {
        // Corrected Base Case: A single leaf is a tree by itself.
        if (leaves.size() == 1) {
            TreeNode leaf = new TreeNode(null);
            leaf.setTipName(leaves.get(0));
            return Collections.singletonList(leaf);
        }

        List<TreeNode> results = new ArrayList<>();
        List<List<List<String>>> partitions = generateBinaryPartitions(leaves);

        for (List<List<String>> partition : partitions) {
            // Each partition has exactly two non-empty subsets
            List<TreeNode> leftSubtrees = buildBinaryTrees(partition.get(0));
            List<TreeNode> rightSubtrees = buildBinaryTrees(partition.get(1));

            // Combine all possible left and right subtrees
            for (TreeNode left : leftSubtrees) {
                for (TreeNode right : rightSubtrees) {
                    TreeNode root = new TreeNode(null);
                    root.addChild(left);
                    root.addChild(right);
                    results.add(root);
                }
            }
        }

        return results;
    }

    /**
     * Generates all unique ways to partition a set of items into two non-empty subsets.
     */
    private List<List<List<String>>> generateBinaryPartitions(List<String> items) {
        List<List<List<String>>> result = new ArrayList<>();
        int n = items.size();
        if (n < 2) return result;

        // Using a HashSet with a canonical signature avoids duplicate partitions
        // e.g., ({A}, {B, C}) is the same partition as ({B, C}, {A})
        Set<String> seenSignatures = new HashSet<>();

        // Iterate through all possible subsets using a bitmask, skipping the empty set and the full set.
        for (int mask = 1; mask < (1 << n) - 1; mask++) {
            List<String> left = new ArrayList<>();
            List<String> right = new ArrayList<>();

            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    left.add(items.get(i));
                } else {
                    right.add(items.get(i));
                }
            }

            String signature = canonicalSignature(left, right);
            if (seenSignatures.add(signature)) {
                result.add(Arrays.asList(left, right));
            }
        }

        return result;
    }

    /**
     * Creates a canonical string representation of a partition to detect duplicates.
     */
    private String canonicalSignature(List<String> a, List<String> b) {
        List<String> list1 = new ArrayList<>(a);
        List<String> list2 = new ArrayList<>(b);
        Collections.sort(list1);
        Collections.sort(list2);

        String sig1 = String.join(",", list1) + "|" + String.join(",", list2);
        String sig2 = String.join(",", list2) + "|" + String.join(",", list1);

        // Return the lexicographically smaller signature
        return sig1.compareTo(sig2) < 0 ? sig1 : sig2;
    }
}