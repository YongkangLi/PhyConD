package com.yongkangl.phylogenetics.substitution;

import com.yongkangl.phylogenetics.util.Mutation;
import com.yongkangl.phylogenetics.util.Taxon;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BlockwiseSubstitutionModel {
    private final SiteModel siteModel;
    private final Taxon initial;
    private final int length;
    private final int contextLength;
    private final int maxBlockSize;
    private final int[] boundaries;
    private final int effectiveLength;
    private final BlockwiseSiteModel[] blockwiseSiteModels;

    private Taxon current;

    public static double[][] getKernel(SiteModel siteModel, Taxon taxon, int start, int end) {
        int contextLength= siteModel.getContextLength();
        int[] leftNeighbors = taxon.subSequence(start - contextLength, start);
        int[] rightNeighbors = taxon.subSequence(end, end + contextLength);
        return BlockwiseSiteModel.getKernel(end - start, leftNeighbors, rightNeighbors, siteModel);
    }

    public BlockwiseSubstitutionModel(SiteModel siteModel, int contextLength, int maxBlockSize, double threshold, Taxon initial) {
        this.siteModel = siteModel;
        this.length = initial.len();
        this.contextLength = contextLength;
        this.initial = new Taxon(length);
        this.initial.copySequence(initial);
        this.maxBlockSize = maxBlockSize;

        SubstitutionModel substitutionModel = new SubstitutionModel(siteModel, initial);
        boolean[] divisibilities = new boolean[initial.len() + 1];
        Arrays.fill(divisibilities, true);
        for (int i = contextLength; i < length - contextLength; i++) {
            if (substitutionModel.getMutability(i) > threshold) {
                for (int j = i - contextLength + 1; j <= i + contextLength; j++) {
                    divisibilities[j] = false;
                }
            }
        }
        divisibilities[contextLength] = true;
        divisibilities[length - contextLength] = true;
        List<BlockwiseSiteModel> blockwiseSiteModels = new ArrayList<>();
        List<Integer> boundaries = new ArrayList<>();
        for (int start = contextLength; start < length - contextLength;) {
            int end = start + 1;
            while (!divisibilities[end] && end <= length - contextLength) {
                end++;
            }
            boundaries.add(end);
            blockwiseSiteModels.add(new BlockwiseSiteModel(getKernel(siteModel, initial, start, end)));
            start = end;
        }
        this.effectiveLength = blockwiseSiteModels.size();
        this.blockwiseSiteModels = new BlockwiseSiteModel[effectiveLength];
        for (int i = 0; i < effectiveLength; i++) {
            this.blockwiseSiteModels[i] = blockwiseSiteModels.get(i);
        }
        this.boundaries = new int[boundaries.size() + 1];
        this.boundaries[0] = contextLength;
        for (int i = 0; i < boundaries.size(); i++) {
            this.boundaries[i + 1] = boundaries.get(i);
        }

        current = new Taxon(this.length);
        current.copySequence(initial);
    }

    public BlockwiseSubstitutionModel (BlockwiseSubstitutionModel blockwiseSubstitutionModel) {
        this.siteModel = blockwiseSubstitutionModel.getSiteModel();
        this.initial = blockwiseSubstitutionModel.getInitial();
        this.length = initial.len();
        this.maxBlockSize = blockwiseSubstitutionModel.getMaxBlockSize();
        this.contextLength = blockwiseSubstitutionModel.getContextLength();
        this.effectiveLength = blockwiseSubstitutionModel.getEffectiveLength();
        this.blockwiseSiteModels = new BlockwiseSiteModel[effectiveLength];
        for (int i = 0; i < effectiveLength; i++) {
            this.blockwiseSiteModels[i] = blockwiseSubstitutionModel.getBlockwiseSiteModel(i);
        }
        this.boundaries = blockwiseSubstitutionModel.getBoundaries();
        this.current = new Taxon(length);
        this.current.copySequence(this.initial);
    }

    public Taxon getInitial() {
        return initial;
    }

    public int getMaxBlockSize() {
        return maxBlockSize;
    }

    public int[] getBoundaries() {
        return boundaries;
    }

    public int getContextLength() {
        return contextLength;
    }

    public int getEffectiveLength() {
        return effectiveLength;
    }

    public BlockwiseSiteModel getBlockwiseSiteModel(int i) {
        return blockwiseSiteModels[i];
    }

    public double logLikelihood(Taxon initial, Taxon terminal, double T) {
        double ll = 0.0;
        for (int i = 0; i < effectiveLength; i++) {
            int start = boundaries[i];
            int end = boundaries[i + 1];
            ll += blockwiseSiteModels[i].logLikelihood(initial.getEncodedBlock(start, end), terminal.getEncodedBlock(start, end), T);
        }
        return ll;
    }

    public SiteModel getSiteModel() {
        return siteModel;
    }

    public double getMutabilitiesSum() {
        double sum = 0.0;
        for (int i = 0; i < effectiveLength; i++) {
            sum += blockwiseSiteModels[i].getMutability(BlockwiseSiteModel.encode(current.subSequence(boundaries[i], boundaries[i+1])));
        }
        return sum;
    }

    public double getRate(int site, int to) {
        int index = 0;
        while (site >= boundaries[index + 1]) {
            index++;
        }
        int start = boundaries[index];
        int end = boundaries[index + 1];
        int original = current.get(site);
        int from = BlockwiseSiteModel.encode(current.subSequence(start, end));
        current.set(site, to);
        int toEncoded = BlockwiseSiteModel.encode(current.subSequence(start, end));
        current.set(site, original);
        return blockwiseSiteModels[index].getRate(from, toEncoded);
    }

    public double mutate(int site, int to) {
        double rate = getRate(site, to);
        current.set(site, to);
        return rate;
    }

    public void reset() {
        current.copySequence(initial);
    }

    public double logLikelihoodFullPath(List<Mutation> mutations, double T) {
        // THIS IS VERY IMPORTANT!!!
        reset();
        double ll = 0;
        double time = 0;
        for (Mutation mutation : mutations) {
            int site = mutation.getSite();
            ll -= (mutation.getTime() - time) * getMutabilitiesSum();
            ll += Math.log(mutate(site, mutation.getTo()));
            time = mutation.getTime();
        }
        ll -= (T - time) * getMutabilitiesSum();
        return ll;
    }

    public double geometricPathMutability(SubstitutionModel substitutionModel, int site, double temperature) {
        double mutability = 0.0;
        int current = substitutionModel.getCurrent().get(site);
        for (int to = 0; to < 4; to++) {
            if (to != current) {
                mutability += Math.pow(getRate(site, to), 1 - temperature) * Math.pow(substitutionModel.getRate(site, to), temperature);
            }
        }
        return mutability;
    }
}
