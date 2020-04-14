#!/usr/bin/env python

import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='APRIL script. Annotate Potential disease Related genes by Integrating 3D chromatin Loops')
    parser.add_argument("-f1", type = str, required = True, help = "path to the chromatin contact data" )
    parser.add_argument("-i", type = str, required = True, help = "cell line index")
    parser.add_argument("-t", default = "1", type = str, required = False, help = "TF expression threshold. DEFAULT: 1")
    parser.add_argument("-n", default = "default", type = str, required = False, help = "number of clusters. DEFAULT: number of modules/20")
    parser.add_argument("-f2", type = str, required = True, help = "path to the SNP effect size data")
    parser.add_argument("-m", default = "BOTH", type = str, required = False, help = "methods using in label propagation. options: GENE, OTHER, BOTH. DEFAULT: BOTH")
    parser.add_argument("-tree", default="50", type=str, required=False, help="number of trees using in random forest. DEFAULT: 50")
    parser.add_argument("-f3", type=str, required=True, help="path to the disease gene data")
    args = parser.parse_args()
    return args

def build_graphs(args):
    os.system("Rscript 1_network_construction.R {0}".format(args.f1))
    os.system("Rscript 2_annotate_fragment.R")
    os.system("Rscript 3_extract_frag_TF_matrix.R {0} {1}".format(args.i, args.t))
    os.system("Rscript 4_merge_network.R {0}".format(args.n))

def predict_gene(args):
    os.system("Rscript 5_annotate_effect.R {0}".format(args.f2))
    os.system("Rscript 6_label_propagation.R {0}".format(args.m))
    os.system("Rscript 7_random_forest.R {0} {1} {2}".format(args.i, args.tree, args.f3))

def main():
    args = parse_args()
    os.system("mkdir -p ../output")
    build_graphs(args)
    predict_gene(args)

main()