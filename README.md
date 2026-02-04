# Selection Pressure Analysis
____________________________________________________________________________________________________________________________________

To identify the selection pressures affecting the H5N1 2.3.4.4b clade's genome.

- The *Dataset* provided is after eliminating all the poor quality sequences and is meant to be used as the input for _Cawlign_, the first step in the pipeline.

Creating the Tree for the Analysis -
-   Running **Cawlign** on all the sequences to perform codon-level amino-acid level alignments directly
-   **tn93 cluster** _(-l 500 -t 0.005 -f)_ to select clusters all closely related sequences (all at 0.5% or less nucleotide distance from every other member of the cluster) and represent them with a single sequence.
-   Create a tree using IQTree2 _(-T 8 -m GTR+I+G)_ and annotate branches belonging to the focal clade (2.3.3.4b) with 2.3.3.4b using the [LabelTree.bf script](https://github.com/veg/hyphy/blob/master/res/TemplateBatchFiles/lib/label-tree.bf) in HyPhy v2.5.62

The Analysis -
-   **Branch-Site Unrestricted Statistical Test for Episodic Diversification BUSTED[S]** method is to seek gene-wide evidence of episodic diversifying selection (EDS) in the focal clade.
-   **Fixed Effects Likelihood (FEL)** applied to internal branches labeled 2.3.3.4b to identify individual sites evolving non-neutrally
-   **Mixed Effects Model of Evolution (MEME) method** applied to internal branches labeled 2.3.3.4b to identify individual sites evolving subject to episodic diversifying positive selection
-   **Contrast-FEL** applied to internal branches labeled 2.3.3.4b or others to identify individual sites that have significantly different ω (i.e., evolve subject to different selective pressures) between the focal (2.3.3.4b) and reference clades.
-   **RELAX** method to compare the intensity of selective forces acting on the 2.3.4.4b clade and reference clades.
____________________________________________________________________________________________________________________________________
### Authors:
Aniket Naik, Yoshita Jadalanki

----------------------------------
```diff
+ Aniket Naik | aniketsnaik7@gmail.com | ann279@pitt.edu
+ 
````
