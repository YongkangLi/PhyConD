package com.yongkangl.phylogenetics.io;

import com.yongkangl.phylogenetics.util.Taxon;

public class BeastTreeLineParser {
    private String input;
    private int position;
    private int traitDimension;
    private double lnp;
    private double joint;
    private String auxiliary;
    private int treeNumber;

    public BeastTreeLineParser(String input) {
        this.input = input;
        this.position = 0;
        this.lnp = 0.0;
        this.joint = 0.0;
        this.treeNumber = 0;
        // Initialize partials array as needed, possibly requiring additional parsing logic
    }

    public int getTreeNumber() {
        return treeNumber;
    }

    public double getLnp() {
        return lnp;
    }

    public double getJoint() {
        return joint;
    }

    public String getAuxiliary() {
        return auxiliary;
    }

    public TreeNode parseBeast() {
        parseBeastHeader();
        return parseTree(null);
    }

    public double[][] parsePartials(String partialsStr) {
        // Remove the leading "{{" and trailing "}}" characters
        partialsStr = partialsStr.substring(2, partialsStr.length() - 2);
        // Split the string into rows based on "}, {"
        String[] rows = partialsStr.split("\\}, \\{");
        double[][] partials = new double[rows.length][];

        for (int i = 0; i < rows.length; i++) {
            // Split each row into its elements based on ", "
            String[] elements = rows[i].split(", ");
            partials[i] = new double[elements.length];
            for (int j = 0; j < elements.length; j++) {
                // Parse each element as a double and store it in the array
                partials[i][j] = Double.parseDouble(elements[j]);
            }
            traitDimension = elements.length >> 2;
        }

        return partials;
    }

    private void parseBeastHeader() {
        int start = input.indexOf("tree") + 5; // Skip "tree " part
        int end = input.indexOf(" ", start);
        String treeName = input.substring(start, end);
        try {
            this.treeNumber = Integer.parseInt(treeName.replaceAll("\\D+", ""));
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid tree number format");
        }

        int propertiesStart = input.indexOf("[&", end) + 2; // Skip "[&"
        int propertiesEnd = input.indexOf("]", propertiesStart);
        String properties = input.substring(propertiesStart, propertiesEnd);
        String[] parts = properties.split("=");
        try {
            this.lnp = Double.parseDouble(parts[1].split(",")[0]);
            this.joint = Double.parseDouble(parts[2]);
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number format in tree header");
        }

        this.position = input.indexOf("] ", propertiesEnd + 1) + 2; // Skip "= [&...]" part to start of NEWICK-like tree
        this.auxiliary = input.substring(propertiesEnd + 4, this.position);
    }

    private TreeNode parseTree(TreeNode parent) {
        TreeNode root = new TreeNode(parent);
        while (position < input.length()) {
            char current = input.charAt(position);
            if (current == '(' || current == ',') {
                position++; // Skip '('
                root.addChild(parseTree(root));
            } else if (current == ')') {
                position++; // Skip ')'
                break; // End of this node's children
            } else {
                root.setTipName(parseTipName());
                break;
            }
        }
        parseAbsoluteBranchLength(root);
        return root;
    }

    private TreeNode parsePrecomputedTree(TreeNode parent) {
        TreeNode root = new TreeNode(parent);
        while (position < input.length()) {
            char current = input.charAt(position);
            if (current == '(' || current == ',') {
                position++; // Skip '('
                root.addChild(parseTree(root));
            } else if (current == ')') {
                position++; // Skip ')'
                break; // End of this node's children
            } else {
                root.setTipName(parseTipName());
                break;
            }
        }
        parseCalibratedBranchLength(root);
        return root;
    }

    private int parseNodeNumber() {
        int nodeNumber = -1;
        // Method to extract and handle [NODE_NUMBER] before a node
        if (input.charAt(position) == '[') {
            int end = input.indexOf("]", position);
            String nodeNumberStr = input.substring(position + 1, end);
            try {
                nodeNumber = Integer.parseInt(nodeNumberStr);
                // Use or store nodeNumber as needed for the TreeNode
            } catch (NumberFormatException e) {
                // Handle error or invalid format
            }
            position = end + 1; // Move past the node number
        }
        return nodeNumber;
    }

    private String parseTipName() {
        int next = input.indexOf(":", position);
        String tipName = input.substring(position, next);
        position = next; // Move past the tip number
        return tipName;
    }

    private void parseAbsoluteBranchLength(TreeNode treeNode) {
        if (input.charAt(position) == ':') {
            int nextComma = input.indexOf(",", position);
            int nextBracket = input.indexOf(")", position);
            int next = Integer.compareUnsigned(nextComma, nextBracket) < 0 ? nextComma : nextBracket;
            int nextEqual = input.indexOf("=", position);
            int nextSquareBracket = input.indexOf("]", position);
            double rate = Double.parseDouble(input.substring(nextEqual + 1, nextSquareBracket));
            double branchLength = Double.parseDouble(input.substring(nextSquareBracket + 1, next));
            treeNode.setRate(rate);
            treeNode.setBranchLength(branchLength);
            position = next; // Move past the branch length
        }
    }

    private void parseCalibratedBranchLength(TreeNode treeNode) {
        if (input.charAt(position) == ':') {
            int nextComma = input.indexOf(",", position);
            int nextBracket = input.indexOf(")", position);
            int next = Integer.compareUnsigned(nextComma, nextBracket) < 0 ? nextComma : nextBracket;
            String[] numbers = input.substring(position, next).split("=")[1].split("]");
            double rate = Double.parseDouble(numbers[0]);
            double branchLength = Double.parseDouble(numbers[1]);
            treeNode.setRate(rate);
            treeNode.setBranchLength(branchLength);
            position = next; // Move past the branch length
        }
    }
}