package com.yongkangl.phylogenetics.sampler;

import com.yongkangl.phylogenetics.VDJ.VDJRecombinationModel;
import com.yongkangl.phylogenetics.io.TreeNode;
import com.yongkangl.phylogenetics.model.SMCParticle;
import com.yongkangl.phylogenetics.substitution.*;

import java.util.concurrent.ExecutorService;

public class HiddenNodes extends SMCParticle {
    private final VDJRecombinationModel vdjRecombinationModel;
    private final SequenceNode root;

    public HiddenNodes(int nSteps, TreeNode root, int mutationSteps, SiteModel[] modelSequence) {
        super(nSteps);
        this.vdjRecombinationModel = new VDJRecombinationModel("../resources/10X");
        this.root = new SequenceNode(null, root, mutationSteps, modelSequence, vdjRecombinationModel);
    }

    public void copyFrom(HiddenNodes from) {
        super.copyFrom(from);
        root.copyFrom(from.root);
    }

    public void mutate() {
        root.mutate();
    }

    @Override
    public void moveForward() {
        super.moveForward();
        root.moveForward();
    }

    @Override
    public double calculateStepLogWeight(int step) {
        root.calculateStepLogWeight();
        return root.getStepLogWeight(false);
    }

    public String mutations() {
        return root.mutations();
    }

    @Override
    public String toString() {
        return root.logFullPi() + "\n" + root.mutations() + root.sequences();
    }
}