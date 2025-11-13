package com.yongkangl.phylogenetics.sampler;

import com.yongkangl.phylogenetics.substitution.BlockwiseSiteModel;
import com.yongkangl.phylogenetics.substitution.SiteModel;
import com.yongkangl.phylogenetics.substitution.SubstitutionModel;
import com.yongkangl.phylogenetics.util.Taxon;
import com.yongkangl.phylogenetics.util.Utils;
import jeigen.DenseMatrix;
import org.apache.commons.math3.distribution.EnumeratedDistribution;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SplitTaxonSampler {
    private final int length;
    private List<Integer> boundaries;
    private List<BlockSampler> blockSamplers;
    private int effectiveLength;
    private double[][] pmf;
    private Taxon fixedTaxon;

    public SplitTaxonSampler (Taxon taxon, int contextLength, List<Integer> boundaries) {
        this.fixedTaxon = taxon;
        this.length = taxon.len();
        this.boundaries = boundaries;
        this.effectiveLength = boundaries.size() - 1;
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            int start = boundaries.get(i);
            int end = boundaries.get(i + 1);
            blockSamplers.add(new BlockSampler(fixedTaxon.subSequence(start - contextLength, start), fixedTaxon.subSequence(start, end), fixedTaxon.subSequence(end, end + contextLength)));
        }
        pmf = calculatePmf();
    }

    public SplitTaxonSampler (Taxon taxon, int contextLength) {
        this.fixedTaxon = taxon.len() % contextLength == 0 ? taxon : new Taxon(taxon.subSequence(0, taxon.len() - taxon.len() % contextLength));
        this.length = taxon.len();
        this.boundaries = new ArrayList<>();
        for (int i = 1; i <= length / contextLength - 1; i++) {
            boundaries.add(i * contextLength);
        }
        this.effectiveLength = boundaries.size() - 1;
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            int start = boundaries.get(i);
            int end = boundaries.get(i + 1);
            blockSamplers.add(new BlockSampler(fixedTaxon.subSequence(start - contextLength, start), fixedTaxon.subSequence(start, end), fixedTaxon.subSequence(end, end + contextLength)));
        }
        pmf = calculatePmf();
    }

    public SplitTaxonSampler (SiteModel proposal, Taxon taxon, int contextLength, double threshold) {
        this.fixedTaxon = taxon;
        this.length = taxon.len();

        SubstitutionModel substitutionModel = new SubstitutionModel(proposal, taxon);
        boolean[] divisibilities = new boolean[this.length + 1];
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
        List<BlockSampler> blockSamplers = new ArrayList<>();
        this.boundaries = new ArrayList<>();
        boundaries.add(contextLength);
        for (int start = contextLength; start < length - contextLength;) {
            int end = start + 1;
            while (!divisibilities[end] && end <= length - contextLength) {
                end++;
            }
            boundaries.add(end);
            blockSamplers.add(new BlockSampler(fixedTaxon.subSequence(start - contextLength, start), fixedTaxon.subSequence(start, end), fixedTaxon.subSequence(end, end + contextLength)));
            start = end;
        }
        this.effectiveLength = blockSamplers.size();
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            this.blockSamplers.add(blockSamplers.get(i));
        }
        pmf = calculatePmf();
    }

//    public SplitTaxonSampler(Taxon taxon, int blockSize, int offset, int contextLength) {
//        fixedTaxon = taxon;
//        this.length = fixedTaxon.len();
//        this.blockSize = blockSize;
//        offset = offset == 0 ? blockSize : offset;
//        effectiveLength = 1 + (length - 2 * contextLength - offset + blockSize - 1) / blockSize;
//        this.blockSamplers = new BlockSampler[effectiveLength];
//        int start = contextLength;
//        int end = contextLength + offset;
//        for (int j = 0; j < effectiveLength; j++) {
//            blockSamplers[j] = new BlockSampler(fixedTaxon.subSequence(start - contextLength, start), fixedTaxon.subSequence(start, end), fixedTaxon.subSequence(end, end + contextLength));
//            start = end;
//            end = start + blockSize;
//            if (end > length - contextLength && start < length - contextLength) {
//                end = length - contextLength;
//            }
//        }
//        pmf = calculatePmf();
//    }

    public static List<Integer> determineBoundaries(SplitTaxonSampler splitTaxonSampler1, SplitTaxonSampler splitTaxonSampler2) {
       List<Integer> boundaries1 = splitTaxonSampler1.getBoundaries();
       List<Integer> boundaries2 = splitTaxonSampler2.getBoundaries();
       List<Integer> indeces1 = new ArrayList<>();
       List<Integer> indeces2 = new ArrayList<>();
       int i = 1, j = 1;
       while (i < boundaries1.size() && j < boundaries2.size()) {
           if (boundaries1.get(i) < boundaries2.get(j)) {
               i++;
           } else if (boundaries1.get(i) > boundaries2.get(j)) {
               j++;
           } else {
               indeces1.add(i);
               indeces2.add(j);
               i++;
               j++;
           }
       }
       int contextLength = splitTaxonSampler1.blockSamplers.get(0).getContextLength();
       List<Integer> boundaries = new ArrayList<>();
       boundaries.add(boundaries1.get(0));
       for (int index = 0; index < indeces1.size(); index++) {
           boundaries.add(boundaries1.get(indeces1.get(index)));
       }
//       for (i = 0; i < indeces1.size(); i++) {
//           boolean keep = true;
//           int[] r1 = splitTaxonSampler1.blockSamplers.get(indeces1.get(i) - 1).getRightNeighbors();
//           int[] r2 = splitTaxonSampler2.blockSamplers.get(indeces2.get(i) - 1).getRightNeighbors();
//           for (j = 0; j < contextLength; j++) {
//               if (r1[j] != r2[j]) {
//                   keep = false;
//                   break;
//               }
//           }
//           if (!keep) {
//               continue;
//           }
//           if (indeces1.get(i) < indeces1.size()) {
//               int[] l1 = splitTaxonSampler1.blockSamplers.get(indeces1.get(i)).getLeftNeighbors();
//               int[] l2 = splitTaxonSampler2.blockSamplers.get(indeces2.get(i)).getLeftNeighbors();
//               for (j = 0; j < contextLength; j++) {
//                   if (l1[j] != l2[j]) {
//                       keep = false;
//                       break;
//                   }
//               }
//               if (!keep) {
//                   continue;
//               }
//           }
//           boundaries.add(boundaries1.get(indeces1.get(i)));
//       }
       return boundaries;
    }

    public void reconcileBoundaries(List<Integer> newBoundaries) {
        int start = 0;
        int end = 1;
        List<BlockSampler> newBlockSamplers = new ArrayList<>();
        for (int i = 1; i < newBoundaries.size(); i++) {
            while (boundaries.get(end) < newBoundaries.get(i)) {
                end++;
            }
            newBlockSamplers.add(BlockSampler.concatenateBlockSamplers(blockSamplers, start, end));
            start = end;
            end = start + 1;
        }
        this.effectiveLength = newBlockSamplers.size();
        this.blockSamplers = newBlockSamplers;
        this.boundaries = newBoundaries;
    }

    public void simplify() {
        int contextLength = blockSamplers.get(0).getContextLength();
        for (int i = 0; i < effectiveLength; i++) {
            if (blockSamplers.get(i).getBlockSize() >= 5) {
                int start = boundaries.get(i);
                int end = boundaries.get(i + 1);
                List<Integer> potentials = new ArrayList<>();
                int[] states = new int[end - start + 2 * contextLength];
                for (int site = 0; site < contextLength; site++) {
                    states[site] = blockSamplers.get(i).getLeftNeighbor(site);
                }
                for (int site = start; site < end; site++) {
                    double max;
                    double second;
                    int state;
                    if (pmf[site][0] > pmf.length) {
                        max = pmf[site][0];
                        second = pmf[site][1];
                        state = 0;
                    } else {
                        max = pmf[site][1];
                        second = pmf[site][0];
                        state = 1;
                    }
                    for (int j = 2; j < 4; j++) {
                        if (pmf[site][j] > max) {
                            second = max;
                            max = pmf[site][j];
                            state = j;
                        }
                    }
                    if (max > 2 * second) {
                        potentials.add(site);
                    }
                    states[site - start + contextLength] = state;
                }
                for (int site = 0; site < contextLength; site++) {
                    states[contextLength + end - start + site] = blockSamplers.get(i).getRightNeighbor(site);
                }
                for (int j = contextLength; j < potentials.size() - contextLength; j++) {
                    boolean splittable = true;
                    for (int k = j - contextLength; k < j + contextLength - 1; k++) {
                        if (potentials.get(k) + 1 != potentials.get(k + 1)) {
                            splittable = false;
                            break;
                        }
                    }
                    if (splittable) {
                        int site = potentials.get(j);
                        int length1 = site - start;
                        Pair<BlockSampler, BlockSampler> splitted = blockSamplers.get(i).split(length1, Arrays.copyOfRange(states, j - contextLength, j + contextLength));
                        blockSamplers.remove(i);
                        blockSamplers.add(i, splitted.getFirst());
                        blockSamplers.add(i + 1, splitted.getSecond());
                        effectiveLength++;
                        boundaries.add(i + 1, site);
                        break;
                    }
                }
            }
        }
    }

    public SplitTaxonSampler(SplitTaxonSampler splitTaxonSampler, DenseMatrix transition, boolean backward) {
        this.length = splitTaxonSampler.getLength();
        this.boundaries = new ArrayList<>();
        this.boundaries.addAll(splitTaxonSampler.getBoundaries());
        this.effectiveLength = splitTaxonSampler.getEffectiveLength();
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            this.blockSamplers.add(new BlockSampler(splitTaxonSampler.blockSamplers.get(i), transition, backward));
        }
        pmf = calculatePmf();
    }

    public SplitTaxonSampler(SplitTaxonSampler splitTaxonSampler, SiteModel siteModel, double T, boolean backward) {
        this.length = splitTaxonSampler.getLength();
        this.boundaries = new ArrayList<>();
        this.boundaries.addAll(splitTaxonSampler.getBoundaries());
        this.effectiveLength = splitTaxonSampler.getEffectiveLength();
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            this.blockSamplers.add(new BlockSampler(splitTaxonSampler.blockSamplers.get(i), siteModel, T, backward));
        }
        pmf = calculatePmf();
    }

    public SplitTaxonSampler(SplitTaxonSampler[] splitTaxonSamplers) {
        this.length = splitTaxonSamplers[0].getLength();
        this.boundaries = splitTaxonSamplers[0].getBoundaries();
        this.effectiveLength = splitTaxonSamplers[0].getEffectiveLength();
        this.blockSamplers = new ArrayList<>();
        for (int i = 0; i < effectiveLength; i++) {
            BlockSampler[] messages = new BlockSampler[splitTaxonSamplers.length];
            for (int j = 0; j < splitTaxonSamplers.length; j++) {
                messages[j] = splitTaxonSamplers[j].blockSamplers.get(i);
            }
            this.blockSamplers.add(new BlockSampler(messages));
        }
        pmf = calculatePmf();
//        simplify();
    }

    public double[][] calculatePmf () {
        double[][] pmf = new double[length][4];
        int contextLength = blockSamplers.get(0).getContextLength();
        for (int site = 0; site < contextLength; site++) {
            Arrays.fill(pmf[site], 0.0);
            pmf[site][blockSamplers.get(0).getLeftNeighbor(site)] = 1.0;
        }
        for (int index = 0; index < effectiveLength; index++) {
            double[][] individualPmfs = blockSamplers.get(index).getIndividualPmfs();
            int start = boundaries.get(index);
            int end = boundaries.get(index + 1);
            for (int site = start; site < end; site++) {
                double[] individualPmf = individualPmfs[site - start];
                System.arraycopy(individualPmf, 0, pmf[site], 0, 4);
                double sum = Arrays.stream(pmf[site]).sum();
                for (int state = 0; state < 4; state++) {
                    pmf[site][state] /= sum;
                }
            }
        }
        for (int site = boundaries.get(effectiveLength); site < boundaries.get(effectiveLength) + contextLength; site++) {
            Arrays.fill(pmf[site], 0.0);
            pmf[site][blockSamplers.get(effectiveLength - 1).getRightNeighbor(site - boundaries.get(effectiveLength))] = 1.0;
        }
        return pmf;
    }

    public void iterate (SplitTaxonSampler child, SiteModel siteModel, double T) {
        int contextLength = blockSamplers.get(0).getContextLength();
        int length = getEffectiveLength();
        int size = 1 << (2 * contextLength);

        double[][] iterated = new double[length][size];

        double[][] parentBlockLogPmf = getBlockLogPmf();
        double[][] childBlockLogPmf = child.getBlockLogPmf();

        int leftMost = BlockwiseSiteModel.encode(child.getLeftMost());
        int rightMost = BlockwiseSiteModel.encode(child.getRightMost());

        Arrays.fill(iterated[0], Double.NEGATIVE_INFINITY);
        for (int rightState = 0; rightState < size; rightState++) {
            DenseMatrix kernel = new DenseMatrix(BlockwiseSiteModel.getKernel(contextLength, BlockwiseSiteModel.decode(contextLength, leftMost), BlockwiseSiteModel.decode(contextLength, rightState), siteModel));
            DenseMatrix transition = kernel.mul(T).mexp();
            for (int startState = 0; startState < size; startState++) {
                for (int endState = 0; endState< size; endState++) {
                    iterated[0][startState] = Utils.logAdd(iterated[0][startState],  parentBlockLogPmf[1][rightState] + Math.log(transition.get(startState, endState)) + childBlockLogPmf[0][endState]);
                }
            }
        }
        double lsl = Utils.logSum(iterated[0]);
        for (int state = 0; state < size; state++) {
            iterated[0][state] -= lsl;
        }

        for (int block = 1; block < length - 1; block++) {
            Arrays.fill(iterated[block], Double.NEGATIVE_INFINITY);
            for (int leftState = 0; leftState < size; leftState++) {
                for (int rightState = 0; rightState < size; rightState++) {
                    DenseMatrix kernel = new DenseMatrix(BlockwiseSiteModel.getKernel(contextLength, BlockwiseSiteModel.decode(contextLength, leftState), BlockwiseSiteModel.decode(contextLength, rightState), siteModel));
                    DenseMatrix transition = kernel.mul(T).mexp();
                    for (int startState = 0; startState < size; startState++) {
                        for (int endState = 0; endState< size; endState++) {
                            iterated[block][startState] = Utils.logAdd(iterated[block][startState], parentBlockLogPmf[block - 1][leftState] + parentBlockLogPmf[block + 1][rightState] + Math.log(transition.get(startState, endState)) + childBlockLogPmf[block][endState]);
                        }
                    }
                }
            }
            double logSum = Utils.logSum(iterated[block]);
            for (int state = 0; state < size; state++) {
                iterated[block][state] -= logSum;
            }
        }

        Arrays.fill(iterated[length - 1], Double.NEGATIVE_INFINITY);
        for (int leftState = 0; leftState < size; leftState++) {
            DenseMatrix kernel = new DenseMatrix(BlockwiseSiteModel.getKernel(contextLength, BlockwiseSiteModel.decode(contextLength, leftState), BlockwiseSiteModel.decode(contextLength, rightMost), siteModel));
            DenseMatrix transition = kernel.mul(T).mexp();
            for (int startState = 0; startState < size; startState++) {
                for (int endState = 0; endState < size; endState++) {
                    iterated[length - 1][startState] = Utils.logAdd(iterated[length - 1][startState], parentBlockLogPmf[length - 2][leftState] + Math.log(transition.get(startState, endState)) + childBlockLogPmf[length - 1][endState]);
                }
            }
        }
        double lsr = Utils.logSum(iterated[length - 1]);
        for (int state = 0; state < size; state++) {
            iterated[length - 1][state] -= lsr;
        }

        fillBlockPmfs(iterated);
    }

    public void fillBlockPmfs(double[][] blockLogPmfs) {
        for (int i = 0; i < effectiveLength; i++) {
            blockSamplers.get(i).fillLogPmf(blockLogPmfs[i]);
        }
        pmf = calculatePmf();
    }

    public double getLogLikelihood() {
        double ll = 0.0;
        for (BlockSampler blockSampler : blockSamplers) {
            ll += blockSampler.getLogLikelihood();
        }
        return ll;
    }

    public Taxon sample() {
        int[] states = new int[length];
        int contextLength = blockSamplers.get(0).getContextLength();
        System.arraycopy(blockSamplers.get(0).getLeftNeighbors(), 0, states, 0, contextLength);
        for (int site = contextLength, i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            System.arraycopy(blockSamplers.get(i).sample(), 0, states, site, b);
            site += b;
        }
        System.arraycopy(blockSamplers.get(effectiveLength - 1).getRightNeighbors(), 0, states, length - contextLength, contextLength);
        return new Taxon(states);
    }

    public Taxon sample(int start, int end) {
        int[] states = new int[length];
        int contextLength = blockSamplers.get(0).getContextLength();
        System.arraycopy(blockSamplers.get(0).getLeftNeighbors(), 0, states, 0, contextLength);
        for (int site = blockSamplers.get(0).getContextLength(), i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            if (start - b < site && site < start && site + b <= end) {
                System.arraycopy(blockSamplers.get(i).sample(), start - site, states, start, site + b - start);
            } else if (start <= site && site + b <= end) {
                System.arraycopy(blockSamplers.get(i).sample(), 0, states, site, b);
            } else if (start <= site && site < end && site + b > end) {
                System.arraycopy(blockSamplers.get(i).sample(), 0, states, site, end - site);
            }
            site += b;
        }
        System.arraycopy(blockSamplers.get(effectiveLength - 1).getRightNeighbors(), 0, states, length - contextLength, contextLength);
        return new Taxon(states);
    }

    public double evaluateLogLikelihood(Taxon taxon) {
        int[] states = taxon.getSequence();
        double ll = 0.0;
        for (int site = blockSamplers.get(0).getContextLength(), i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            ll += blockSamplers.get(i).evaluateLogLikelihood(Arrays.copyOfRange(states, site, site + b));
            site += b;
        }
        return ll;
    }

    public double evaluateLogLikelihood(Taxon taxon, int start, int end) {
        double ll = 0.0;
        for (int site = blockSamplers.get(0).getContextLength(), i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            if (start - b < site && site < start && site + b <= end) {
                ll += blockSamplers.get(i).evaluateRightLogLikelihood(Arrays.copyOfRange(taxon.getSequence(), start, site + b));
            } else if (start <= site && site + b <= end) {
                ll += blockSamplers.get(i).evaluateLogLikelihood(Arrays.copyOfRange(taxon.getSequence(), site, site + b));
            } else if (start <= site && site < end && site + b > end) {
                ll += blockSamplers.get(i).evaluateLeftLogLikelihood(Arrays.copyOfRange(taxon.getSequence(), site, end));
            }
            site += b;
        }
        return ll;
    }

    public double evaluateMeanHammingDistance(Taxon taxon) {
        double distance = 0.0;
        for (int site = 0; site < length; site++) {
            distance += 1.0 - (pmf[site][taxon.get(site)] / Arrays.stream(pmf[site]).sum());
        }
        return distance;
    }

    public Taxon conditionalSample(Taxon taxon, SiteModel siteModel, double T) {
        int contextLength = blockSamplers.get(0).getContextLength();
        SplitTaxonSampler fixed = new SplitTaxonSampler(siteModel, taxon, contextLength, 3.0);
        SplitTaxonSampler evolved = new SplitTaxonSampler(fixed, siteModel, T, false);
        List<Integer> newBoundaries = determineBoundaries(this, evolved);
        reconcileBoundaries(newBoundaries);
        evolved.reconcileBoundaries(newBoundaries);
        int[] states = new int[length];
        System.arraycopy(blockSamplers.get(0).getLeftNeighbors(), 0, states, 0, contextLength);
        for (int site = contextLength, i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            List<Pair<Integer, Double>> possibleStates = new ArrayList<>();
            int size = 1 << (2 * b);
            for (int state = 0; state < size; state++) {
                possibleStates.add(new Pair<>(state, Math.exp(blockSamplers.get(i).getLogPmfValue(state) + evolved.blockSamplers.get(i).getLogPmfValue(state))));
            }
            EnumeratedDistribution<Integer> distribution = new EnumeratedDistribution<>(possibleStates);
            System.arraycopy(BlockwiseSiteModel.decode(b, distribution.sample()), 0, states, site, b);
            site += b;
        }
        System.arraycopy(blockSamplers.get(effectiveLength - 1).getRightNeighbors(), 0, states, length - contextLength, contextLength);
        return new Taxon(states);
    }

    public double evaluateConditionalLogLikelihood(Taxon conditionalTaxon, Taxon taxon, SiteModel siteModel, double branchLength) {
        int contextLength = blockSamplers.get(0).getContextLength();
        SplitTaxonSampler fixed = new SplitTaxonSampler(conditionalTaxon, contextLength, boundaries);
        SplitTaxonSampler evolved = new SplitTaxonSampler(fixed, siteModel, branchLength, false);
        int[] states = taxon.getSequence();
        double ll = 0.0;
        for (int site = contextLength, i = 0; i < effectiveLength; i++) {
            int b = blockSamplers.get(i).getBlockSize();
            int size = 1 << (2 * b);
            double[] logProbabilities = new double[size];
            for (int state = 0; state < size; state++) {
                logProbabilities[state] = blockSamplers.get(i).getLogPmfValue(state) + evolved.blockSamplers.get(i).getLogPmfValue(state);
            }
            int state = BlockwiseSiteModel.encode(Arrays.copyOfRange(states, site, site + b));
            ll += logProbabilities[state] / Utils.logSum(logProbabilities);
            site += b;
        }
        return ll;
    }

    public double[][] getPmf() {
        return pmf;
    }

    public double[][] getBlockLogPmf() {
        double[][] blockLogPmf = new double[effectiveLength][];
        for (int i = 0; i < effectiveLength; i++) {
            blockLogPmf[i] = blockSamplers.get(i).getLogPmf();
        }
        return blockLogPmf;
    }

    public int[] getLeftMost() {
        return blockSamplers.get(0).getLeftNeighbors();
    }

    public int[] getRightMost() {
        return blockSamplers.get(effectiveLength - 1).getRightNeighbors();
    }

    public int getLength() {
        return length;
    }

    public List<Integer> getBoundaries() {
        return boundaries;
    }

    public int getEffectiveLength() {
        return effectiveLength;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        for (int i = 0; i < length; i++) {
            sb.append("[");
            for (int j = 0; j < 4; j++) {
                sb.append(pmf[i][j]);
                if (j < 3) {
                    sb.append(",");
                }
            }
            sb.append("]");
            if (i < length - 1) {
                sb.append(",");
            }
        }
        sb.append("]");
        return sb.toString();
    }
}
