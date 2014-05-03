#!/usr/bin/python

import csv
import os
import sys
import glob
import re
from collections import defaultdict

def open_csv_file(file):
	csvfile = open(file, 'rb')
	csvdat = csv.reader(csvfile, delimiter='\t', quotechar='"')
	csvdat.next()
	return csvdat

def open_csv_file_comma(file):
	csvfile = open(file, 'rb')
	csvdat = csv.reader(csvfile, delimiter=',', quotechar='"')
	csvdat.next()
	return csvdat

def cmap_drug_list(cmap_file,num):
	cmap_drugs=[]
	for lines in cmap_file:
		drugs = lines[num].upper()
		cmap_drugs.append(drugs)
	cmap_drugs = list(set(cmap_drugs))
	return cmap_drugs

def cmap_upanddown_list(cmap_file):
	up_genes=[]
	down_genes=[]

	for lines in cmap_file:
		if float(lines[1]) < 0 :
			genes = lines[0]	
			down_genes.append(genes)
		elif float(lines[1]) > 0:
			genes = lines[0]	
			up_genes.append(genes)				
	up_genes = list(set(up_genes))
	down_genes = list(set(down_genes))
	return up_genes,down_genes

def drug_targ_dic(drug_targ):
	ddict = defaultdict(list)
	geneids_drugs = {}
	for druglines in drug_targ:
		hgnc_id = druglines[2]
		druglist = druglines[12].split(';')
		#print hgnc_id
		try:
			for drugs in druglist:
				hgnc_id.upper()
				ddict[hgnc_id].append(drugs)
		except:
			pass
		geneids_drugs = dict((k, tuple(v)) for k, v in ddict.iteritems())
	return geneids_drugs


def main():
	expres_file = open_csv_file("exrpression_results_final1000.txt")
	#gene_sys_cmap = open_csv_file("gene_symbols.txt")

	#express_list = cmap_drug_list(expres_file, 7)
	#cmap_sys_list = cmap_drug_list(gene_sys_cmap,1)

	#intersect_drugs =  list(set(express_list).intersection(cmap_sys_list))
	#drug_targ = open_csv_file_comma("drugbank/all_target_ids_with_known_action.csv")
	#drug_targ_dict = drug_targ_dic(drug_targ)
	'''
	for genes in intersect_drugs:
		try:
			print drug_targ_dict[genes]
		except:
			pass
	'''
	up_genes,down_genes = cmap_upanddown_list(expres_file)
	up_regulated_genes = open("up_regulated_genes.grp", 'wb')
	down_regulated_genes = open("down_regulated_genes.grp", 'wb')
	up_regulated_genes.write("\n".join(up_genes))
	down_regulated_genes.write("\n".join(down_genes))

	'''
	up_regulated_genes = open("up_regulated_genes.grp", 'wb')
	down_regulated_genes = open("down_regulated_genes.grp", 'wb')

	up_regulated_genes_grp = open_csv_file("up_regulated_probesets.grp")
	down_regulated_genes_grp = open_csv_file("down_regulated_probesets.grp")

	up_genes =[]
	down_genes=[]
	for lines in up_regulated_genes_grp:
		if lines[1]:
			up_genes.append(lines[1])
	for lines in down_regulated_genes_grp:
		if lines[1]:
			down_genes.append(lines[1])

	up_genes = list(set(up_genes))
	down_genes=list(set(down_genes))

	inter_genes = list(set(up_genes).intersection(set(down_genes)))
	up_genes = list(set(up_genes) - set(inter_genes))
	down_genes = list(set(down_genes) - set(inter_genes))

	up_regulated_genes.write("\n".join(up_genes))
	down_regulated_genes.write("\n".join(down_genes))
	'''

if __name__ == "__main__":
        main()
