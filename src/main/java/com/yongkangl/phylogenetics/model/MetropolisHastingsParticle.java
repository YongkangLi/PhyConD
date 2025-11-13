package com.yongkangl.phylogenetics.model;

import org.apache.commons.math3.random.RandomDataGenerator;

public abstract class MetropolisHastingsParticle extends Particle {

    private int accepted = 0;

    public abstract MetropolisHastingsParticle propose();
    public abstract double logPiOverQ();
    public abstract void copyFrom(MetropolisHastingsParticle from);
    public abstract void log();

    @Override
    public void step() {
        MetropolisHastingsParticle proposal = propose();
        double metropolisHastingsRatio = Math.exp(proposal.logPiOverQ() - logPiOverQ());
        double u = new RandomDataGenerator().nextUniform(0.0, 1.0);
        if (u < metropolisHastingsRatio) {
            copyFrom(proposal);
            accepted++;
        }
        log();
    }

    public int getAccepted() {
        return accepted;
    }
}
