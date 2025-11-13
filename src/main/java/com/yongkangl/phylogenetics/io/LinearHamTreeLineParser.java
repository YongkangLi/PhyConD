package com.yongkangl.phylogenetics.io;

import com.yongkangl.phylogenetics.substitution.SiteModel;
import com.yongkangl.phylogenetics.util.Taxon;

public class LinearHamTreeLineParser {
    private String input;
    private int position;
    private int traitDimension;
    private String auxiliary;

    public LinearHamTreeLineParser(String input) {
        this.input = input;
        this.position = 0;
    }

    public TreeNode parseLinearHam() {
        TreeNode root = parseTree(null);
        if (root.getChild(0).getTipName().equals("naive")) {
            TreeNode naive = root.getChild(0);
            root.setFixedTaxon(naive.getFixedTaxon());
            double branchLength = naive.getBranchLength();
            root.getChild(1).setBranchLength(branchLength + root.getChild(1).getBranchLength());
            root.removeChild(0);
        }
        return root;
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
                String tipName = parseTipName();
                root.setTipName(tipName);
                break;
            }
        }
        Taxon taxon = parseTaxon();
        if (taxon != null) {
            root.setFixedTaxon(taxon);
        }
        root.setRate(1.0);
        root.setBranchLength(parseBranchLength());
        return root;
    }

    private String parseTipName() {
        int nextSemicolon = input.indexOf(":", position);
        int nextBracket = input.indexOf("[", position);
        int next = Integer.compareUnsigned(nextSemicolon, nextBracket) < 0 ? nextSemicolon : nextBracket;
        String tipName = input.substring(position, next);
        position = next;
        return tipName;
    }

    private Taxon parseTaxon() {
        Taxon taxon = null;
        if (input.charAt(position) == '[') {
            int startingQuote = input.indexOf("\"", position);
            int endingQuote = input.indexOf("\"", startingQuote + 1);
            taxon = new Taxon(input.substring(startingQuote + 1, endingQuote));
            position = input.indexOf("]", endingQuote) + 1;
        }
        return taxon;
    }

    private double parseBranchLength() {
        double branchLength = 0.0;
        if (input.charAt(position) == ':') {
            int nextComma = input.indexOf(",", position);
            int nextBracket = input.indexOf(")", position);
            int next = Integer.compareUnsigned(nextComma, nextBracket) < 0 ? nextComma : nextBracket;
            String number = input.substring(position + 1, next);
            branchLength = Double.parseDouble(number);
            position = next; // Move past the branch length
        }
        return branchLength;
    }
}