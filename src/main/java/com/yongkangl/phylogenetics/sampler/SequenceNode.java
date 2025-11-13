package com.yongkangl.phylogenetics.sampler;

import com.yongkangl.phylogenetics.VDJ.VDJRecombinationModel;
import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.model.SMCParticle;
import com.yongkangl.phylogenetics.substitution.*;
import com.yongkangl.phylogenetics.util.Taxon;
import org.apache.commons.math3.random.RandomDataGenerator;
import jeigen.DenseMatrix;

import java.util.ArrayList;
import java.util.List;

public class SequenceNode extends SMCParticle {
    private final VDJRecombinationModel vdjRecombinationModel;
    private Taxon fixedTaxon;
    private double logVDJprob;
    private final double branchLength;
    private final SequenceNode parent;
    private final List<SequenceNode> children;
    private double logWeight = Double.NEGATIVE_INFINITY;
    private final int mutationSteps;
    private MetropolisHastingsPath path;

    private final SiteModel[] modelSequence;

    public SequenceNode(SequenceNode parent, TreeNode node, int mutationSteps, SiteModel[] modelSequence, VDJRecombinationModel vdjRecombinationModel) {
        super(modelSequence.length - 1);
        this.vdjRecombinationModel = vdjRecombinationModel;
        this.branchLength = node.getBranchLength();
        this.parent = parent;
        this.children = new ArrayList<>();
        this.mutationSteps = mutationSteps;
        this.modelSequence = modelSequence;

        fixedTaxon = new Taxon(node.getFixedTaxon());
        if (parent == null) {
            logVDJprob = vdjRecombinationModel.computeVDJProb(fixedTaxon.toString());
        }
        if (parent != null) {
            this.path = new MetropolisHastingsPath(parent.fixedTaxon, fixedTaxon, branchLength, (IndependentSiteModel) modelSequence[0], modelSequence[1]);
        }
        for (int i = 0; i < node.getChildCount(); i++) {
            TreeNode child = node.getChild(i);
            addChild(new SequenceNode(this, child, mutationSteps, modelSequence, vdjRecombinationModel));
        }
    }

    public void copyFrom(SequenceNode from) {
        fixedTaxon.copySequence(from.fixedTaxon);
        logWeight = from.logWeight;
        if (parent != null) {
            path.copyFrom(from.path);
        }
        for (int i = 0; i < children.size(); i++) {
            children.get(i).copyFrom(from.children.get(i));
        }
    }

    public void mutate() {
        for (int i = 0; i < mutationSteps; i++) {
            mutationStep();
        }
    }

    public void mutationStep() {
        if (!children.isEmpty()) {
            for (SequenceNode child : children) {
                child.path.setTarget(modelSequence[child.getCurrentStep()]);
                if (!child.children.isEmpty()) {
                    child.mutationStep();
                }
            }
        }

        int contextLength = modelSequence[getCurrentStep()].getContextLength();
        int len = fixedTaxon.len();
        for (int site = contextLength; site < len - contextLength; site++) {
            if (parent == null && ((site < (290 - 27)) || (site > len - 21))) {
                continue;
            }
            int originalState = fixedTaxon.get(site);
            double originalLogVDJprob = logVDJprob;
            double logAcceptanceRatio = 0.0;
            if (parent != null) {
                logAcceptanceRatio -= path.hobolthForward(site);
            } else {
                logAcceptanceRatio -= originalLogVDJprob;
            }
            for (SequenceNode child: children) {
                logAcceptanceRatio -= child.path.hobolthBackward(site);
            }
            List<DenseMatrix> childrenBackwards = new ArrayList<>();
            List<Integer> childrenStates = new ArrayList<>();
            for (SequenceNode child : children) {
                childrenBackwards.add(child.path.getFirstBackwardMatrice());
                childrenStates.add(child.path.getTerminalState(site));
            }
            int state;
            if (parent != null) {
                logAcceptanceRatio += path.hobolthEndpointedOnChildren(site, path.getLastForwardMatrice(), childrenBackwards, childrenStates);
                state = path.getTerminalState(site);
            } else {
                state = MetropolisHastingsPath.hobolthEndpointedOnChildren(childrenBackwards, childrenStates);
            }
            fixedTaxon.set(site, state);
            if (parent == null) {
                logVDJprob = vdjRecombinationModel.computeVDJProb(fixedTaxon.toString());
                logAcceptanceRatio += logVDJprob;
            }
            for (SequenceNode child : children) {
                child.path.setInitial(site, state);
                logAcceptanceRatio += child.path.hobolthEndpointed(site, child.path.getBackwardMatrices());
            }
            boolean accept = new RandomDataGenerator().nextUniform(0.0, 1.0) <= Math.exp(logAcceptanceRatio);
            if (Double.isNaN(Math.exp(logAcceptanceRatio))) {
                System.err.println("position: " + site);
                System.err.println("Current VDJ: " + originalLogVDJprob);
                System.err.println("Propose VDJ: " + logVDJprob);
            }
            if (!accept) {
                fixedTaxon.set(site, originalState);
            }
            if (parent != null) {
                path.accept(accept);
                if (!accept) {
                    path.setTerminal(site, originalState);
                }
            }
            for (SequenceNode child : children) {
                child.path.accept(accept);
                if (!accept) {
                    child.path.setInitial(site, originalState);
                    logVDJprob = originalLogVDJprob;
                }
            }
        }
    }

    public double calculateStepLogWeight(int step) {
        logWeight = 0.0;
        if (parent != null) {
            logWeight = path.logImportanceWeight(modelSequence[step], modelSequence[step + 1]);
        }
        for (SequenceNode child : children) {
            logWeight += child.calculateStepLogWeight(step);
        }
        return logWeight;
    }

    @Override
    public void moveForward() {
        super.moveForward();
        for (SequenceNode child : children) {
            child.moveForward();
        }
    }

    @Override
    public double getStepLogWeight(boolean advance) {
        double currentLogWeight = logWeight;
        if (advance) {
            moveForward();
        }
        return currentLogWeight;
    }

    public void setFixedTaxon (Taxon taxon) {
        this.fixedTaxon = taxon;
    }

    public void addChild(SequenceNode child) {
        children.add(child);
    }

    public void setLogWeight(double logWeight) {
        this.logWeight = logWeight;
    }

    public Taxon getFixedTaxon() {
        return fixedTaxon;
    }

    public double getBranchLength() {
        return branchLength;
    }

    public int getChildCount() {
        return children.size();
    }

    public SequenceNode getChild(int i) {
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

    public String mutations() {
        StringBuilder sb = new StringBuilder();
        if (parent != null) {
            sb.append(path.toString()).append("\n");
        }
        if (!children.isEmpty()) {
            for (SequenceNode child : children) {
                sb.append(child.mutations());
            }
        }
        return sb.toString();
    }

    public double logFullPi() {
        double logPi = 0.0;
        if (parent != null) {
            logPi += path.logFullPi();
        }
        for (SequenceNode child : children) {
            logPi += child.logFullPi();
        }
        return logPi;
    }

    public String evolutionaryHistory(int level) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < level; i++) {
            sb.append("    ");
        }
        if (parent != null) {
            sb.append(branchLength);
            sb.append(": ");
            sb.append(path.toString());
        }
        sb.append("\n");
        if (!children.isEmpty()) {
            for (SequenceNode child : children) {
                sb.append(child.evolutionaryHistory(level + 1));
            }
        }
        return sb.toString();
    }
}