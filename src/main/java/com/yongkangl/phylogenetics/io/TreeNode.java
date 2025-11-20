package com.yongkangl.phylogenetics.io;

import com.yongkangl.phylogenetics.model.MarkovChain;
import com.yongkangl.phylogenetics.sampler.SplitTaxonSampler;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import jeigen.DenseMatrix;

import java.util.ArrayList;
import java.util.List;

public class TreeNode {
    private String tipName;
    private Taxon fixedTaxon;
    private SplitTaxonSampler taxonSampler;
    private double rate = 0.0;
    private double branchLength;
    private TreeNode parent;
    private List<TreeNode> children;
    public SplitTaxonSampler[] messagesFromChildren;
    private double logWeight = Double.NEGATIVE_INFINITY;
    private MetropolisHastingsPath[] paths;
    private MarkovChain<MetropolisHastingsPath> metropolisHastings;

    public TreeNode(TreeNode parent) {
        this.parent = parent;
        this.children = new ArrayList<>();
    }

    public void pairSequences(FastaParser fastaParser) {
        if (isTip()) {
            fixedTaxon = new Taxon(fastaParser.queryTaxon(tipName));
        } else {
            for (TreeNode child : children) {
                child.pairSequences(fastaParser);
            }
        }
    }

    public void sampleClone(SiteModel siteModel) {
        for (TreeNode child : children) {
            if (!child.isTip()) {
                child.sample(siteModel, fixedTaxon, false);
                child.sampleClone(siteModel);
            }
        }
    }

    public void setUpMessages(SiteModel siteModel, int contextLength) {
        setUpMessagesFromChildren(siteModel, contextLength);
    }


    private SplitTaxonSampler createTaxonSampler (SplitTaxonSampler[] splitTaxonSamplers) {
        assert splitTaxonSamplers.length == 2;
        return new SplitTaxonSampler(splitTaxonSamplers);
    }

    private SplitTaxonSampler createTaxonSampler (SplitTaxonSampler splitTaxonSampler, SiteModel siteModel, double T) {
        DenseMatrix kernel = new DenseMatrix(BlockwiseSiteModel.getKernel(1, new int[]{}, new int[]{}, siteModel));
        DenseMatrix transition = kernel.mul(T).mexp();
        return new SplitTaxonSampler(splitTaxonSampler, transition, true);
    }

    private SplitTaxonSampler createTaxonSampler (SplitTaxonSampler splitTaxonSampler, DenseMatrix transition) {
        return new SplitTaxonSampler(splitTaxonSampler, transition, true);
    }

    public void setUpTaxonSampler() {
        if (taxonSampler == null) {
            taxonSampler = createTaxonSampler(messagesFromChildren);
        }
    }

    public void setUpMessagesFromChildren(SiteModel siteModel, int contextLength) {
        if (getChildCount() > 0) {
            for (int i = 0; i < getChildCount(); i++) {
                getChild(i).setUpMessagesFromChildren(siteModel, contextLength);
            }
            messagesFromChildren = new SplitTaxonSampler[getChildCount()];
            for (int i = 0; i < getChildCount(); i++) {
                TreeNode child = getChild(i);
                messagesFromChildren[i] = createTaxonSampler(child.getTaxonSampler(), siteModel, child.getBranchLength());
            }
            setUpTaxonSampler();
        }
        else {
            messagesFromChildren = new SplitTaxonSampler[1];
            List<Integer> boundaries = new ArrayList<>();
            for (int i = contextLength; i <= fixedTaxon.len() - contextLength; i++) {
                boundaries.add(i);
            }
            taxonSampler = new SplitTaxonSampler(fixedTaxon, contextLength, boundaries);
            messagesFromChildren[0] = taxonSampler;
        }
    }

    public void setUpMessagesFromChildren(DenseMatrix transition) {
        if (getChildCount() > 0) {
            for (int i = 0; i < getChildCount(); i++) {
                getChild(i).setUpMessagesFromChildren(transition);
            }
            messagesFromChildren = new SplitTaxonSampler[getChildCount()];
            for (int i = 0; i < getChildCount(); i++) {
                TreeNode child = getChild(i);
                messagesFromChildren[i] = createTaxonSampler(child.getTaxonSampler(), transition);
            }
            setUpTaxonSampler();
        }
        else {
            messagesFromChildren = new SplitTaxonSampler[1];
            List<Integer> boundaries = new ArrayList<>();
            boundaries.add(2);
            boundaries.add(fixedTaxon.len() - 2);
            taxonSampler = new SplitTaxonSampler(fixedTaxon, 2, boundaries);
            messagesFromChildren[0] = taxonSampler;
        }
    }

    public void calculateExactly(SiteModel siteModel, boolean fixingUCA, Taxon UCA) {
        calculateExactlyUpward(siteModel);
        if (fixingUCA) {
            int contextLength = siteModel.getContextLength();
            List<Integer> boundaries = new ArrayList<>();
            boundaries.add(contextLength);
            boundaries.add(UCA.len() - contextLength);
            taxonSampler = new SplitTaxonSampler(UCA, contextLength, boundaries);
            calculateExactlyDownward(siteModel);
        }
    }

    public void calculateExactlyUpward(SiteModel siteModel) {
        int contextLength = siteModel.getContextLength();
        if (getChildCount() > 0) {
            for (int i = 0; i < getChildCount(); i++) {
                getChild(i).calculateExactlyUpward(siteModel);
            }
            messagesFromChildren = new SplitTaxonSampler[getChildCount()];
            for (int i = 0; i < getChildCount(); i++) {
                TreeNode child = getChild(i);
                messagesFromChildren[i] = new SplitTaxonSampler(child.getTaxonSampler(), siteModel, child.getBranchLength(), true);
            }
            setUpTaxonSampler();
        }
        else {
            messagesFromChildren = new SplitTaxonSampler[1];
            List<Integer> boundaries = new ArrayList<>();
            boundaries.add(contextLength);
            boundaries.add(fixedTaxon.len() - contextLength);
            taxonSampler = new SplitTaxonSampler(fixedTaxon, contextLength, boundaries);
            messagesFromChildren[0] = taxonSampler;
        }
    }

    public void calculateExactlyDownward(SiteModel siteModel) {
        if (parent != null) {
            SplitTaxonSampler[] messages = new SplitTaxonSampler[2];
            messages[0] = new SplitTaxonSampler(parent.taxonSampler, siteModel, branchLength, false);
            messages[1] = taxonSampler;
            taxonSampler = new SplitTaxonSampler(messages);
        }
        for (int i = 0; i < getChildCount(); i++) {
            getChild(i).calculateExactlyDownward(siteModel);
        }
    }

    public double getConditionalLogLikelihood() {
        return taxonSampler.getLogLikelihood();
    }

    public void setFixedTaxon (Taxon taxon) {
        this.fixedTaxon = taxon;
    }

    public void setTipName(String tipName) {
        this.tipName = tipName;
    }

    public void setRate(double rate) {
        this.rate = rate;
    }

    public void setBranchLength(double branchLength) {
        this.branchLength = branchLength;
    }

    public void addChild(TreeNode child) {
        children.add(child);
    }

    public void removeChild(int i) {
        children.remove(i);
    }

    public SplitTaxonSampler getTaxonSampler() {
        return taxonSampler;
    }

    public void setLogWeight(double logWeight) {
        this.logWeight = logWeight;
    }

    public Taxon getFixedTaxon() {
        return fixedTaxon;
    }

    public void sample() {
        fixedTaxon = taxonSampler.sample();
    }

    public Taxon sample(int start, int end) {
        return taxonSampler.sample(start, end);
    }

    public void sample(SiteModel siteModel, Taxon parent, boolean fixed) {
        if (!(children.isEmpty() || fixed)) {
            fixedTaxon = taxonSampler.conditionalSample(parent, siteModel, branchLength);
        }
    }

    public double evaluateTargetLogLikelihood(Taxon taxon) {
        return taxonSampler.evaluateLogLikelihood(taxon);
    }

    public double evaluateTargetLogLikelihood(Taxon taxon, int start, int end) {
        return taxonSampler.evaluateLogLikelihood(taxon, start, end);
    }

    public double evaluateTargetConditionalLogLikelihood(Taxon conditionalTaxon, Taxon sampledTaxon, SiteModel siteModel, double T) {
        return taxonSampler.evaluateConditionalLogLikelihood(conditionalTaxon, sampledTaxon, siteModel, T);
    }

    public double getBranchLength() {
        return rate * branchLength;
    }

    public int getChildCount() {
        return children.size();
    }

    public TreeNode getChild(int i) {
        return children.get(i);
    }

    public boolean isTip() {
        return children.isEmpty();
    }

    public String constructNewick() {
        StringBuilder sb = new StringBuilder();
        if (logWeight != Double.NEGATIVE_INFINITY) {
            sb.append("[weight=").append(Math.exp(logWeight));
        }
        sb.append("] = ");
        sb.append(this);
        sb.append(";");
        return sb.toString();
    }

    public String getTipName() {
        return tipName;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        if (!children.isEmpty()) {
            sb.append("(");
            for (int i = 0; i < children.size(); i++) {
                if (i > 0) sb.append(",");
                sb.append(children.get(i).toString());
            }
            sb.append(")");
        } else {
            sb.append(tipName);
        }
        if (rate > 0.0) {
            sb.append(":").append(branchLength);
        }
        return sb.toString();
    }

    public String sequences() {
        StringBuilder sb = new StringBuilder();
        sb.append(fixedTaxon.toString()).append("\n");
        if (!children.isEmpty()) {
            for (int i = 0; i < children.size(); i++) {
                sb.append(children.get(i).sequences());
            }
        }
        return sb.toString();
    }

    public String pmfToString() {
        return taxonSampler.toString();
    }
}