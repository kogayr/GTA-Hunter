# GTA_Hunter
  GTA_Hunter is a 'support-vector machine' (SVM) classifier that distinguishes gene transfer agent (GTA) genes from their viral homologs. 
  
  The classifier was developed to detect GTA genes that are homologous to 11 of the genes that encode the GTA in *Rhodobacter capsulatus* (RcGTA) (see Kogay et al, 2019). However, the classifier can be trained to detect other RcGTA genes, other GTAs, or even other virus-like elements, as long as there is training data available and the two classes have distinct amino acid composition. 

## System Requirements
 `Python v3.5.1` 

## Python packages
- NumPy v. 1.13.1
- CVXOPT v. 1.1.8
- matplotlib v. 2.1.0
- pandas v. 0.21.0
- Biopython v. 1.69

## Third party requirements
- BLAST v. 2.2.31+

**GTA_Hunter.py is a command-line program with the following options:**

```
usage: GTA_Hunter_blast.py [-h] [-b] [-g GTA] [-v VIRUS] [-q QUERIES]
                           [-o OUTDIR] [-k [KMER]] [-p [PSEAAC]] [-y] [-m]
                           [-O] [-f FOLDER] [-W] [-w WEIGHT WEIGHT]
                           [-t CLUSTER_TYPE] [-d DIST] [-c C] [-x [XVAL]]
                           [-e KERNEL KERNEL] [-s]

Gene Classification Using SVM.

optional arguments:
  -h, --help            show this help message and exit
  -b, --blast           Run BLASTP search to first identify GTA homologs for
                        the classification.
  -g GTA, --GTA GTA     The FASTA-formatted (.faa) "true GTA" sequences
                        used for training.
  -v VIRUS, --virus VIRUS
                        The FASTA-formatted (.faa) "true viruses"
                        sequences used for training.
  -q QUERIES, --queries QUERIES
                        The FASTA-formatted (.faa) sequences to be
                        classified.
  -o OUTDIR, --outdir OUTDIR
                        The folder path in which to store output.
  -k [KMER], --kmer [KMER]
                        The size of k-mer (default=4).
  -p [PSEAAC], --pseaac [PSEAAC]
                        Expand feature set to include pseudo amino acid
                        composition (default=3). As the parameter, specify
                        value of lambda. Weight = 0.05 (after Chou 2001).
  -y, --physico         Expand feature set to include physicochemical
                        properties of amino acids.
  -m, --min             Print bare minimum results.
  -O, --optimal         Use the optimal parameters for the RcGTA gene homolog
                        classification as listed in Table 2 in Kogay et al
                        2019.
  -f FOLDER, --folder FOLDER
                        Provide a folder with one or multiple proteomes (*.faa
                        files).
  -W                    Weight training datasets. Pre-calculated distances 
  			will be used automatically.
  -w WEIGHT WEIGHT, --weight WEIGHT WEIGHT
                        Weight training datasets. Specify the two pairwise
                        distance files needed for weighting (first file for
                        GTAs, second file for viruses).
  -z CLUSTER_TYPE, --cluster_type CLUSTER_TYPE
                        Specify 'farthest' or 'nearest' neighbors clustering
                        (default='farthest').
  -t DIST, --dist DIST  Specify the cutoff distance for clustering in the
                        weighting scheme (default=0.01).
  -c C, --soft_margin C
                        The soft margin for the SVM (default=1.0).
  -x [XVAL], --xval [XVAL]
                        Performs cross validation of training set. Specify
                        folds over 10 repetitions (default=5).
  -e KERNEL KERNEL, --kernel KERNEL KERNEL
                        (for advanced users) Specify kernel to be used and sigma for the gaussian kernel.
                        (options: 'linear', 'gaussian') (default='linear', 0).
  -s, --svs             Show support vectors.
  
```

## How to use it - Example 1: Classify homologs of RcGTA-like gene using the provided training data
In this example, we will classify five homologs of the RcGTA portal protein (g3) using k-mer of size 3 as a feature and setting raw softness of the SVM margin to 100. The sequences to classify are in FASTA format in one file, and are located in the "example_run" folder; the training data is in "data/training" folder. 

```
python GTA_Hunter.py -g data/training/gta/3_gta.faa -v data/training/viral/3_viral.faa -c 10000 -w data/training/gta/3_gta.dist data/training/viral/3_viral.dist -t 0.02 -k 3 -q example_run/g3_example.faa
```

In the output, the first two sequences are classified as a "GTA" and the remaining three are classified as a "virus":

```
Gene                                                                                           Score          Classification    GTA Gene
>WP_100283140.1 phage portal protein [Sphingomonas sp. Cra20]                                   -0.824313      GTA               3
>WP_100367008.1 phage portal protein [Yoonia maricola]                                          -0.947548      GTA               3
>WP_102855116.1 phage portal protein [Phaeobacter inhibens]                                     0.873095       virus             3
>WP_121219731.1 phage portal protein [Oceanibaculum indicum]                                    0.965013       virus             3
>WP_128262196.1 phage portal protein [Bradyrhizobium yuanmingense]                              0.702121       virus             3
time to run: 10.004424810409546
```

The obtained scores in the "Score" column are derived as described at Figure 2 in Kogay et al. 2019.

The output is in the file "example_output/example_1.txt".

## How to use it - Example 2: Cross-validate classifier's performance for a gene

In this example, we will use five-fold cross-validation to examine the accuracy of the classifier to separate viral and GTA homologs of gene g4 with a specific set of parameters.
```
python GTA_Hunter.py -g data/training/gta/4_gta.faa -v data/training/viral/4_viral.faa -c 10000 -w data/training/gta/4_gta.dist data/training/viral/4_viral.dist -t 0.02 -k 3 -p -x 5
```
In the output, we found that on average all sequences in the training dataset are correctly classified, indicating that classifer has 100% of accuracy for the gene g4.

```
We correctly classified (on average) 62.00/62 GTA and 465.00/465 Viral genes.
time to run: 61.76529502868652
```

The output is in the file "example_output/example_2.txt".

## How to use it - Example 3: Find homologs of RcGTA-like gene in a genome and classify them using the provided training data and the most optimal SVM parameters
In this example, we will first identify homologs of the RcGTA genes among the encoded proteins in a *Rhodobacter sphaeroides* strain HJ genome using a BLASTP search (Altschul et al. 1997), and then classify them using training data provided in data/taining folder and the optimal parameters that were detected via cross validation (see Table 2 in Kogay et al. 2019). The input file is FASTA-formatted file with amino acid sequences of encoded proteins (*.faa format provided by GenBank@NCBI) and is located in the "example_blast" folder. 

Note: please create 'example_outdir' folder before running the program.

```
python GTA_Hunter.py -b -f example_blast/ -o example_outdir/ -O
```
The results will be written to the file 'result_(file_name).out' in the output directory.
This exampe search identified 27 RcGTA homologs and 17 of were classified as RcGTA-like genes.

```
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457497.1 hypothetical protein [Rhodobacter sphaeroides]                                  -0.797741      GTA              15
>WP_137457888.1 host specificity protein [Rhodobacter sphaeroides]                              -1.700612      GTA              15
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457108.1 phage portal protein [Rhodobacter sphaeroides]                                  0.766007       virus             3
>WP_137457895.1 phage portal protein [Rhodobacter sphaeroides]                                  -0.980008      GTA               3
>WP_137459338.1 phage portal protein [Rhodobacter sphaeroides]                                  0.740713       virus             3
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457892.1 DUF3168 domain-containing protein [Rhodobacter sphaeroides]                     -1.814797      GTA               8
>WP_137458772.1 DUF3168 domain-containing protein [Rhodobacter sphaeroides]                     0.071541       virus             8
>WP_137459210.1 DUF3168 domain-containing protein [Rhodobacter sphaeroides]                     0.243832       virus             8
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457496.1 peptidase [Rhodobacter sphaeroides]                                             -0.412156      GTA              14
>WP_137457889.1 peptidase [Rhodobacter sphaeroides]                                             -1.302067      GTA              14
Gene                                                                                           Score          Classification    GTA Gene
>WP_009564481.1 MULTISPECIES: phage major tail protein, TP901-1 family [Rhodobacter]            -1.020971      GTA               9
Gene                                                                                           Score          Classification    GTA Gene
>WP_101328016.1 MULTISPECIES: phage major capsid protein [Rhodobacter]                          -0.948470      GTA               5
>WP_137457105.1 phage major capsid protein [Rhodobacter sphaeroides]                            0.511603       virus             5
>WP_137459300.1 phage major capsid protein [Rhodobacter sphaeroides]                            0.706092       virus             5
Gene                                                                                           Score          Classification    GTA Gene
>WP_115474016.1 HK97 family phage prohead protease [Rhodobacter sphaeroides]                    -1.040761      GTA               4
>WP_137458228.1 HK97 family phage prohead protease [Rhodobacter sphaeroides]                    -0.244565      GTA               4
>WP_137458765.1 HK97 family phage prohead protease [Rhodobacter sphaeroides]                    -0.759410      GTA               4
>WP_137459215.1 HK97 family phage prohead protease [Rhodobacter sphaeroides]                    0.017206       virus             4
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457103.1 hypothetical protein [Rhodobacter sphaeroides]                                  0.037486       virus             6
>WP_137457893.1 hypothetical protein [Rhodobacter sphaeroides]                                  -0.952592      GTA               6
>WP_137458770.1 phage gp6-like head-tail connector protein [Rhodobacter sphaeroides]            0.344344       virus             6
>WP_137459212.1 phage gp6-like head-tail connector protein [Rhodobacter sphaeroides]            0.085210       virus             6
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457495.1 DUF2163 domain-containing protein [Rhodobacter sphaeroides]                     -1.427688      GTA              13
>WP_137457890.1 DUF2163 domain-containing protein [Rhodobacter sphaeroides]                     -1.316051      GTA              13
Gene                                                                                           Score          Classification    GTA Gene
>WP_115474024.1 TIGR02217 family protein [Rhodobacter sphaeroides]                              -1.328590      GTA              12
>WP_137457494.1 TIGR02217 family protein [Rhodobacter sphaeroides]                              -0.455782      GTA              12
Gene                                                                                           Score          Classification    GTA Gene
>WP_137457896.1 ATP-binding protein [Rhodobacter sphaeroides]                                   -1.470279      GTA               2
```

The output is in the file "example_output/example_3.txt".

## How to use it - Example 4: Are the detected GTA-like genes located in the same neighborhood in a genome?

`GTA_Hunter.py` only evaluates one gene at a time, even if queries come from the same genome or even when the whole genome was scanned. To see if the detected "GTA" genes are found in the same neighborhood in a genome, one can run an additional script (`Clustering_genes.py`; located in the "Clustering genes" folder). The script uses the DBSCAN algorithm (Ester et al. 1996).

**Clustering_genes.py is a command-line program with the following options:**  

```  
usage: Clustering_genes.py [-h] --data DATA --feature FEATURE [-s SIZE]
                           [-e DIST]

Clustering of genes via DBSCAN to get RcGTA-like clusters

optional arguments:
  -h, --help            Show this help message and exit
  --data DATA           File that containes one line with space or tab-separated entries. First entry is the genome or contig ID, 
			the rest are RefSeq IDs of the encoded proteins to cluster
  
  --feature FEATURE     Genome feature table in the NCBI feature table file format.
  -s SIZE, --Cluster SIZE
                        The minimum cluster size to consider as RcGTA-like
                        cluster (default=6).
  -e DIST, --Max DIST   The maximum spatial distance to cluster genes
                        (default=8000)

```

  For example, let's examine 17 RcGTA-like homologs detected in the *Rhodobacter sphaeroides* strain HJ genome from the example 3 above. To call a region an "RcGTA-like cluster", we will require that at least 6 of the inputted genes have no more than 8,000 base pairs between the adjacent genes. 
  
```
python Clustering_genes.py --data example_cl/data.txt --feature example_cl/Rhodobacter_sphaeroides_HJ.txt -s 6 -e 8000
```
  The program found that 11 out of the 17 inputted RcGTA-like genes are in the same neighborhood in the genome, and therefore can be designated as an "RcGTA-like cluster":

```
NZ_CP036419.1 has RcGTA-like cluster; The cluster size is 11 genes
WP_137457888.1  WP_137457889.1  WP_137457890.1  WP_115474024.1  WP_009564481.1  WP_137457892.1  WP_137457893.1  WP_101328016.1  WP_115474016.1  WP_137457895.1  WP_137457896.1
```
The RefSeq IDs of the 11 genes in the cluster are listed. The output is in the file "example_output/example_4.txt".

## References

- Altschul SF, et al. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25: 3389-3402
- Chou KC. 2001. Prediction of protein cellular attributes using pseudo‐amino acid composition. Proteins 43: 246-255.
- Ester M, Kriegel H-P, Sander J, Xu X. 1996. A density-based algorithm for discovering clusters a density-based algorithm for discovering clusters in large spatial databases with noise. In Simoudis E, Han J, Fayyad U editors. Proceedings of the Second International Conference on Knowledge Discovery and Data Mining: AAAI Press. p. 226-231.
- Kogay, R., Neely, T. B., Birnbaum, D. P., Hankel, C. R., Shakya, M., & Zhaxybayeva, O. (2019). Machine-learning classification suggests that many alphaproteobacterial prophages may instead be gene transfer agents. Genome biology and evolution, 11(10), 2941-2953.
