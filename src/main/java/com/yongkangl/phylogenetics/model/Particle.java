package com.yongkangl.phylogenetics.model;

public abstract class Particle {
    public abstract void step();

    public void run(int steps) {
        for (int i = 0; i < steps; i++) {
            step();
        }
    }
}
