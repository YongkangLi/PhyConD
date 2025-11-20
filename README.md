# PhyConD

PhyConD is a program designed to reweight phylogenetic trees in the presence of context-dependent mutations.

## Dependencies

- Java 11.0.8
- Python Libraries:
  - numpy
  - pandas
  - matplotlib
  - logomaker
  - biopython

## Example Usage

### Running SMC for trees sampled from LinearHam

#### Building the artifact

The pre-built artifact `smc.jar` is provided in the `linearham_example/` directory. If you want to modify the source code or build it from source, you may open the project in IntelliJ IDEA and build the artifact `PhyConD:LinearHam` from there.

#### Preparing input files

Use LinearHam to sample trees under the ISM. The sampled trees should be stored in a file named `linearham_run.trees`, which we rename to `treeLines.txt` in this example. Each line in this file corresponds to a single tree in Newick format.

For testing the software, an example input file `treeLines.txt` is provided in the `linearham_example/` directory.

#### Running SMC for a single tree

You may have to rename the output file from LinearHam to `treeLines.txt` as `smc.jar` expects this filename.
```
cd linearham_example
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

#### Running SMC for multiple trees

Basically the same as above, but run the program multiple times for different trees.

In particular, you may use scheduling tools such as slurm to parallelize over multiple trees. For example, using the `--array` option in slurm:
```
cd linearham_example
sbatch --array=1-45 run.sh
```
This will run the SMC algorithm for trees in lines 1 to 45 of `treeLines.txt`, each as a separate job. The output and error files will be stored in the `out` and `err` directories respectively.

#### Processing Output

##### Ranking Trees under ISM/DSM

`linearham_run.log` contains the log-likelihoods of each tree under the ISM and the output files under the `out/` directory contain the log-weights of each tree.

You may rank the trees based on these log-likelihoods or log-weights to compare their rankings under different models:
```
cd example
python scripts/process_lls.py linearham_run.log out 1 45
```
The output contains the rankings of the trees from high to log-likelihood under ISM and DSM, respectively.

##### Reconstructing the Ancestral \& Internal Sequences

You may use the script `scripts/process_smc.py` to process the posterior samples of root and internal sequences output by the SMC algorithm. It will output the logo plots for the specified segments of particular nodes in the tree. Directly running `python scripts/process_smc.py` will show the usage information.

### Running SMC for trees sampled from BEAST

If you want to infer the phylogenetic trees from shorter sequences instead of full length BCRs, you may use BEAST to sample trees under the ISM first, and then reweight the sampled trees under the DSM using the SMC algorithm.

#### Building the artifact

The pre-built artifact `smc.jar` is provided in the `beast_example` directory. You may also build the artifact `PhyConD:Beast` from source using IntelliJ IDEA and rename the artifact jar file to `smc.jar`.

#### Preparing input files

Use BEAST to sample trees under the ISM. The sampled trees should be stored in a file named `clone_name.trees`, where `clone_name` is specified by your BEAST configuration file. We will rename this `clone_name.trees` to `treeLines.txt` in this example. The format of this file is a little bit different from the standard Newick format. 

Run the following command to process this file and get `treeLines.txt`:
```
bash scripts/preprocess.sh clone_name.trees
```
Unlike LinearHam, BEAST does not store sequences in the trees and renamed your sequences. Therefore, you have to provide a separate FASTA file named `seqs.fasta` that contains the observed sequences for the leaves in the trees. Run the following command to process your original FASTA file and get `seqs.fasta`:
```
python scripts/translate.py clone_name.trees your.fasta seqs.fasta
```

In a nutshell, you have to provide two files, `treeLines.txt` and `seqs.fasta`, for running the SMC algorithm.

#### Running SMC for a single tree & multiple trees

Like the LinearHam case, you may run the SMC algorithm for a single tree or multiple trees using the following command:
```
cd beast_example
java -Xss1024m -jar smc.jar -s 32 -p 1024 -m 32 -l 1 1> out/1.out 2> err/1.err
```
Similarly, the arguments have the same meaning as in the LinearHam case.

Use the `--array` option in slurm to parallelize over multiple trees:
```
cd beast_example
sbatch --array=1-45 run.sh
```

## Caution
For 55 sequences with sequence length 50, running the SMC with 1024 particles and 32 mutation steps for a single tree may take around 12 hours on a machine with 45 cores. You should not expect to run SMC algorithm for many large trees within a short time.

## License
This project is licensed under the Apache License 2.0 — see the [LICENSE](LICENSE) file for details.
