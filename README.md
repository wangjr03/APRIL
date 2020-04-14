# APRIL
Genetic association analysis based on long-range enhancer-promoter regulatory networks
## Summary
This package discovers novel and known disease genes by integrating connected graphs of promoter, enhancers, TFs and GWAS SNPs with machine learning methods. It can be applied to any chromatin contact data, e.g. ChIA-PET, Capture-C, Hi-C and computational predicted enhancer-gene interactions.

## Introduction
The program takes the chromatin contact maps and GWAS SNPs as input and discovers novel disease genes. It can be applied to :
1. Discover and visualize high-order chromatin interactions from pair-wise chromatin contact datas.
2. Integrate GWAS SNPs with chromatin interactions data.
3. Predict new disease genes from the chromatin contact data and GWAS SNPs.

## Dependencies
The program is developed under `R/3.5.1` and depends on six R packages: `igraph`, `Matrix`, `randomForest`, `ROCR`, `caret`, `dplyr` and `stringr`.

## Input file formats
1. Chromatin contact maps: The chromatin contact maps generated by ChIA-PET, Capture C, Hi-C and computational predictions should be in tab separated file with five columns:


	| col | abbrv. | type | description |
	| --- | --- | --- | --- |
	| 1 | chr | string | Name of the chromosome |
	| 2 | frag1.start | int | Fragment 1 start |
	| 3 | frag1.enh | int | Fragment 1 enh |
	| 4 | frag2.start | int | Fragment 2 start |
	| 5 | frag2.end | int |Fragment 2 end |

2. SNPs effect size: SNPs effect size file provides the effect score on the exact disease, which should be in tab separated file with three columns:


	col | abbrv. | type | description
	--- | --- | --- | ---
	1 | chr | string | Name of the chromosome
	2 | pos | int | Genetic variant position  
	3 | effect | int | Effect size of the SNP in exact disease

3. Significant eQTLs: significant eQTLs file provides the p-value of significant eQTLs on the exact tissue, which should be in tab separated file with three columns:
		
		
	col | abbrv. | type | description
	--- | --- | --- | ---
	1 | chr | string | Name of the chromosome
	2 | pos | int | Genetic variant position  
	3 | gene | string | Target gene of eQTL
	4 | p-value | float | P-value of eQTL

4. Disease gene annotation: Disease-related genes set is used in supervised learning to generate labelled data. The disease gene file just has column:


	col | abbrv. | type | description
	--- | --- | --- | ---
	1 | gene | string | Name of disease gene


## Pre-calculated data
1. Gene annotation: The comprehensive gene annotation with GENCODE V19 has been integrated with the program. The gene annotation file should have at least seven columns:


	col | abbrv. | type | description
	--- | --- | --- | ---
	1 | chr | string | Name of the chromosome
	2 | start | int | Gene body start
	3 | end | int | Gene body end
	4 | stand | string | Orientation of the gene (+/-)
	5 | gene_type | string | The type of the gene, e.g. protein coding, lincRNA, pseudogene and so on.
	6 | gene_symbol | string | HGNC gene symbol
	7 | ENSEMBL_ID | string | ENSEMBL gene ID

	Note: Using the alternative gene annotation file is not supported with the current version of APRIL, since it will also change multiple other pre-calculated input data. 
With the gene annotation file, the promoters of protein-coding genes are calculated (“../promoter.txt”). The promoters are defined as the +/1 kb flanking region of the TSS.

2. Enhancer annotation: A consensus enhancer annotation of hg19 from the Roadmap (https://personal.broadinstitute.org/meuleman/reg2map/HoneyBadger2-intersect_release/) has been integrated into the software for users’ convenience (“./consensus_enhancer_annotation.txt”). Briefly, it includes ~470,000 enhancers across 127 cell lines/tissues, which are identified with DNase-seq peaks in 53 cell-lines and chromatin states learnt by ChromHMM in 127 cell lines/tissues. 
The enhancer annotation file has 3 columns:


	col | abbrv. | type | description
	--- | --- | --- | ---
	1 | chr | string | Name of the chromosome
	2 | start | int | Enhancer start
	3 | enh | int | Enhancer end


3. TF motif hits:
	The genome-wide motif annotation data from the Roadmap (http://compbio.mit.edu/encode-motifs/matches-with-controls-0.3.txt.gz) has been integrated into the software. Controls (indicated by “_C”) were removed. This results in ~13.6M motif hits for ~500 TFs. The TF motif hits are pre-overlapped with enhancers and promoters. The overlapping result for each chromosome can be found at the “promoter_tf.Rdata” and “enh_tf.Rdata”. Each .Rdata file contains two objects: enh_tf and use_pos. The “enh_tf” object is a matrix, where each row is an enhancer in the corresponding chromosome and each column is a TF. The “use_pos” object is a vector, which contains the index of rows of “enh_tf” object in all enhancers of the chromosome. For example, the first three elements of the “use_pos” in chr1 is 16,18,23, which means the first three rows of “enh_tf” correspond to the 16th, 18th and 23rd enhancers of all enhancers within chr1 in the consensus enhancer annotation file. 

4. The enhancer activity profile (DNase-seq):
	The DNase-seq signal track was downloaded in .bigwig format and then converted to .bedGraph format with bedtools (cite). The averaged DNase-seq signal for each enhancer across 127 cell lines/tissues has been pre-calculated as described in the method section. The result (“/data/enh_act_mat.Rdata”) is organized into a matrix where each row is an enhancer and each column is one of 127 cell lines. The order of enhancers is strictly the same as the enhancer annotation file.

5. The gene expression profile (log2(RPKM)):
	The log2(RPKM) signal for each gene across 127 cell lines/tissues has been pre-calculated. The result (“/data/gene_act_updated.Rdata”) is organized into a matrix, where each row is a gene and each column is one of 127 cell lines. The order of genes is strictly the same as the promoter file.

## Fast run
```APRIL.py```
* `-f1`: required, path to the chromatin contact data
* `-f2`: required, path to the SNP effect size data
* `-f3`: required, path to the disease gene data
* `-f4`: required, path to the significat eQTL data
* `-i`: required, cell line index
* `-t`: not required, TF expression threshold. DEFAULT: 1
* `-n`: not required, number of clusters. DEFAULT: number of modules/20
* `-m`: not required, methods using in label propagation. options: GENE, OTHER, BOTH. DEFAULT: BOTH
* `-tree`: not required, number of trees using in random forest. DEFAULT: 50


## Description of scripts
The software consists of 6 sequential scripts. A wrapper is provided to run the whole pipeline. A detailed description of each piece is provided in case the user only wants to run a part of the pipeline. The order of the scripts is indicated in the name.
1. 1_network_construction.R: This step takes chromatin contact maps as input and constructs a list of connected 3D modules.
	* Input data: chromatin contact maps
	* Output: 
		* A list of adjacency matrices, which is stored as “./output/network_adjacency_matrix.Rdata”. Each adjacency matrix corresponds to one 3D module.
		* A list of strings indicated the coordinates of the DNA fragments (nodes) in each module. The order is the same as that in i. The result is stored as “../output/network_nodes_name.Rdata”.
		* A union set of the genomic fragments (nodes), sorted  by chromosome and genomic locations. The result is stored as “../output/island_loc.txt”.
	* Command line usage: 
	
		`Rscript 1_network_construction.R <path to the chromatin contact data>`

2. 2_annotate_fragment.R: This script annotates each node as enhancers or promoters. The active gene, enhancer and TFs are selected here.
	* Input data: 
		* Coordinates of genomic fragments (“../output/island_loc.txt”)
		* Enhancer coordinates ("../data/consensus_enhancer_coord.txt")
		* Promoter coordinates ("../data/promoter.txt")
	* Output:
		Two objects are integrated in "../output/island_annotation.Rdata": “enh_id” and “gene_id”. Both objects are a list object, which store the index of enhancers/promoters that can overlap with the fragments. The order of the fragments is the same as the input DNA fragments. The index of enhancers/promoters is relative to the consensus enhancer annotation file/promoter file. If the index is 0, it means the fragment cannot overlap with any enhancers/promoters. 
	* Command line usage:
	
		`Rscript 2_annotate_fragment.R`
		
		Note:
		
		the path of input data has been specified in the script.

3. 3_extract_frag_TF_matrix.R: This script identifies TF motif hits within fragments that are annotated as enhancers or promoters.
	* Input data:
		* Consensus enhancer annotation
		* enhancer -TF overlapping results and promoter-TF overlapping results
		* Promoter coordinates
		* Enhancer activity profile
		* Gene expression profile
		* Gene annotation
	* Output:
		* A fragment-TF matrix is stored in "../output/frag_TF_mat.Rdata", where each row is a fragment and each column is an expressed TF in corresponding cell line.
		* A vector contains functional annotation of fragments, which can be “E” (enhancer), “P”(promoter),”T”(TF) and “F”(other fragments)
	* Command line usage:
	
		`Rscript 3_extract_frag_TF_matrix.R <cell_index>`
		
		Notes: 
		
		The first argument corresponds to the Roadmap index of the input cell line. A full list of index of 127 cell lines and descriptions are provided as “ENCODE_cell_type.csv”

4. 4_merge_network.R: This script will merge multiple 3D modules with similar TF profiles together and construct merged sub-networks.
	* Input data:
		* Adjacency matrix for each module.
		* Fragment-TF matrix generated in step 3.
	* Output:
		* Final APRIL subgraphs are stored as igraph objects in “../output/merged_skeleton_igraph.Rdata”. Only TFs linking two different modules are kept, other TFs are masked from the subgraph.
	* Command line usage: 
	
		`Rscript 4_merge_network.R <TF_expression_threshold> <number_of_clusters>`
		
		Notes: 
		
		1. The first argument specifies the threshold of the expression for TFs. The default value will be 1, which is used to generate all results in this paper.
		
		2. The second argument specifies the number of module clusters. We used hierarchical clustering with the “complete” method based on  correlations.  The default value of clusters is “number of modules/20”. The result shown in this paper is generated with default value.

* 5_annotate_effect.R: This script annotates effect size for each genomic fragment based on GWAS SNPs and maps the effect size of nodes in each network.
	* Input data:
		* SNPs effect size data.
		* Coordinates of genomic fragments in step 1 (“../output/island_loc.txt”)
		* APRIL sub-networks igraph object in step 4 (“../output/merged_skeleton_igraph.Rdata”).
	* Output:
		* A union set of the genomic fragments (nodes) , and the added effect size of each fragment is in the last column. The result is stored in “../output/frag_effect.Rdata”
		* A list stores the effect size of nodes in each network. The result is stored in "../output/graph_effect.Rdata".
	* Command line usage: 
		
		`Rscript 5_annotate_effect.R <path to the SNP effect size data>`

* 6_label_propagation.R: This script uses HotNet algorithm to propagate disease genes.
	* Input data:
		* APRIL subgraphs igraph object in step 4 (“../output/merged_skeleton_igraph.Rdata”).
		* The effect size of nodes in each network in step 5 ("../output/graph_effect.Rdata").
		* Promoter coordinates ("../data/promoter.txt").
		* Gene annotation ("../data/gene_annotation_V19.txt").
	* Output:
		* A boxplot shows the performance of label propagation
		* A list shows the rank of prediction score of disease-related genes in each network, which is stored in "../output/propagate_rank.Rdata".
	* Command line usage: 
	
		`Rscript 6_label_propagation.R <method_flag>`
		
		Notes: 
		
		The argument determines which version we use to propagate disease genes. It can be “GENE”, “OTHER” or “BOTH”. “GENE” means model just uses the effect size of genes; “OTHER” means model uses the effect size of other fragments except genes; “BOTH” means model uses the effect size of all fragments. The default value will be “BOTH”.

* 7_random_forest.R: This script uses random forest algorithm to predict disease genes.
	* Input data:
		* APRIL subgraphs igraph object in step 4 (“../output/merged_skeleton_igraph.Rdata”).
		* The fragment-TF matrix in step 3 ("../output/frag_TF_mat.Rdata").
		* The index of enhancers/promoters that can overlap with the fragments in step 2 ("../output/island_annotation.Rdata").
		* Coordinates of DNA fragments in step 1 (“../output/island_loc.txt”).
		* enhancer -TF overlapping results and promoter-TF overlapping results
		* Promoter coordinates ("../data/promoter.txt").
		* Enhancer activity profile ("../data/enh_act_mat.Rdata").
		* Gene expression profile ("../data/gene_act_mat_updated.Rdata").
		* Gene annotation ("../data/gene_annotation_V19.txt").
		* The effect size of nodes in each network in step 5 ("../output/graph_effect.Rdata").
	* Output:
		* An ROC shows the 10-fold cross validation performance
		* The probability of each gene within the network to be the disease-related gene, which is stored in "../output/prob_disease_gene.txt".
	* Command line usage: 
		
		`Rscript 7_random_forest.R <cell line index> <number_trees> <path to the disease gene data> <path to the eQTL data>`
		
		Notes: 
		
		1. The first argument corresponds to the Roadmap index of the input cell line.
			
		2. The second argument corresponds to the number of decision trees used in the random forest algorithm.
