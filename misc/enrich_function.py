#!/usr/bin/python
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import os
import sys

import sys
var1 = sys.argv[1]
var2 = sys.argv[2]
var3 = sys.argv[3]
var4 = sys.argv[4]
var5 = sys.argv[5]

def enrich_function(gene_list, universe, gene_sets, directory, filename):
        ##this function takes in a gene_list that is a df with one column with list of gene symbols. universe is the genes in the background one column list of gene symbols. gene_sets is a gmt file, directory is a path and filename is the name the file that will be saved
        os.chdir(directory)
        glist = gene_list.squeeze().str.strip().to_list() #format the genelist
        glist_universe = universe.squeeze().str.strip().to_list()
        enr2 = gp.enrich(gene_list=glist, # or gene_list=glist
                        gene_sets=gene_sets, # kegg is a dict object
                        background=glist_universe, # or "hsapiens_gene_ensembl", or int, or text file, or a list of genes
                        outdir=None,
                        verbose=True)
        df = pd.DataFrame(enr2.results)
        df.to_csv(filename, index=False, sep = ",")

gene_list = pd.read_csv(var1)
universe = pd.read_csv(var2)
gene_sets = var3
directory = var4
filename = var5
enrich_function(gene_list, universe, gene_sets, directory, filename)
