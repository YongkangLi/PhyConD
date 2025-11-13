package com.yongkangl.phylogenetics.model;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public class MarkovChain<T extends Particle> {
    private ExecutorService executorService;
    private final T[] particles;
    private final int nParticles;

    public MarkovChain(T[] particles) {
        this.nParticles = particles.length;
        this.particles = particles;
    }

    public void run(int steps) {
        ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
        Future<?>[] futures = new Future[nParticles];
        for (int i = 0; i < nParticles; i++) {
            final int index = i;
            futures[i] = executorService.submit(() -> {
                particles[index].run(steps);
            });
        }

        for (Future<?> future : futures) {
            try {
                future.get();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        executorService.shutdown();
    }
}
