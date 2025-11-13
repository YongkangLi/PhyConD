package com.yongkangl.phylogenetics.model;

public abstract class SMCParticle {
    private final int nSteps;
    private final double stepLength;
    private final double[] logWeights;
    private int currentStep;

    public SMCParticle(int nSteps) {
        this.nSteps = nSteps;
        this.stepLength = 1.0 / nSteps;
        logWeights = new double[nSteps];
        currentStep = 0;
    }

    public void copyFrom(SMCParticle from) {
        System.arraycopy(from.logWeights, 0, logWeights, 0, nSteps);
        currentStep = from.currentStep;
    }

    public abstract double calculateStepLogWeight(int step);
    public abstract void mutate();

    public int getCurrentStep() {
        return currentStep;
    }

    public void moveForward() {
        currentStep++;
    }

    public double getCurrentTemperature() {
        return currentStep * stepLength;
    }

    public void calculateStepLogWeight() {
        double currentTemperature = getCurrentTemperature();
        logWeights[currentStep] = calculateStepLogWeight(currentStep);
    }

    public double getStepLogWeight(boolean advance) {
        double currentLogWeight = logWeights[currentStep];
        if (advance) {
            moveForward();
        }
        return currentLogWeight;
    }

    public double[] getLogWeights() {
        return logWeights;
    }

    public double getLogWeight(int i) {
        return logWeights[i];
    }
}
