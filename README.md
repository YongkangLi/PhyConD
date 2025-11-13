# PhyConD

PhyConD is a program designed to reweight phylogenetic trees in the presence of context-dependent mutations.

## Dependencies

Java 11.0.8

## Example Usage

### Running SMC for a single tree
```
cd example
java -Xss1024m -jar smc.jar -s 32 -p 1024 -m 32 -l 1 1> out/1.out 2> err/1.err
```

| Argument | Meaning                                | Typical Value |
|----------|----------------------------------------|---------------|
| -s       | \#steps / \#intermediate distributions | 32            |
| -p       | \#particles                            | 1024          |
| -m       | \#mutation steps after resampling      | 32            |
| -l       | the tree in `treeLines.txt`            | 1             |

In this example, the posterior samples of root and internal sequences (X, Z | g) are output in `out/1.out` for the tree in line 1 of `treeLines.txt` along with the log-weight of the tree after marginalizing out the internal sequences.

Note that this program is already parallelized over multiple cores by default. You may limit the number of threads used by setting the `JAVA_OPTS` environment variable.

### Running SMC for multiple trees

Basically the same as above, but run the program multiple times for different trees.

In particular, you may use scheduling tools such as slurm to parallelize over multiple trees. For example, using the `--array` option in slurm:
```
cd example
sbatch --array=1-45 run.sh
```
This will run the SMC algorithm for trees in lines 1 to 45 of `treeLines.txt`, each as a separate job. The output and error files will be stored in the `out` and `err` directories respectively.

### Processing Output

#### Ranking Trees under ISM/DSM

`linearham_run.log` contains the log-likelihoods of each tree under the ISM and the output files under the `out/` directory contain the log-weights of each tree.

You may rank the trees based on these log-likelihoods or log-weights to compare their rankings under different models:
```
cd example
python scripts/process_lls.py linearham_run.log out 1 45
```
The output contains the rankings of the trees from high to log-likelihood under ISM and DSM, respectively.

#### Reconstructing the Ancestral \& Internal Sequences

You may use the script `scripts/process_smc.py` to process the posterior samples of root and internal sequences output by the SMC algorithm. It will output the logo plots for the specified segments of particular nodes in the tree. Directly running `python scripts/process_smc.py` will show the usage information.

