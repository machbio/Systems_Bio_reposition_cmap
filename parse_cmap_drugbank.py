#!/usr/bin/python

import csv
import os
import sys
import glob
import re
from collections import defaultdict

def open_csv_file(file):
	csvfile = open(file, 'rb')
	csvdat = csv.reader(csvfile, delimiter=',', quotechar='"')
	csvdat.next()
	return csvdat

def cmap_drug_list(cmap_file):
	cmap_drugs=[]
	for lines in cmap_file:
		drugs = lines[2].upper()
		cmap_drugs.append(drugs)
	cmap_drugs = list(set(cmap_drugs))
	return cmap_drugs

def drugbank_list(drug_id_file):
	bank_drugs = {}
	drugbank_names = []
	for lines in drug_id_file:
		bank_drugs[lines[1].upper()] = lines[0].upper()
		drugbank_names.append(lines[1].upper())
	#bank_drugs = list(set(bank_drugs))
	return bank_drugs, drugbank_names

def common_drugs(cmap_drugs,drugbank_names, bank_drugs):
	common_drugs_ls = []
	intersect_drugs =  list(set(cmap_drugs).intersection(drugbank_names))
	for drugs in intersect_drugs:
		common_drugs_ls.append(bank_drugs[drugs])
	return common_drugs_ls

def drug_targ_dic(drug_targ):
	ddict = defaultdict(list)
	geneids_drugs = {}
	for druglines in drug_targ:
		hgnc_id = druglines[10].split(':')
		druglist = druglines[12].split(';')
		try:
			for drugs in druglist:
				ddict[hgnc_id[1]].append(drugs)
		except:
			pass
		geneids_drugs = dict((k, tuple(v)) for k, v in ddict.iteritems())
	return geneids_drugs

def prot_target_comm(geneids_drugs, common_drugs_ls):
	geneid_target = []
	for drugs in common_drugs_ls:
		if (drugs in [x for v in geneids_drugs.values() for x in v]):
			geneid_target.extend(geneids_drugs.keys())
	return geneid_target

def main():
	#cmap_file = open_csv_file("cmap/cmap_instances.csv")
	#drug_id_file = open_csv_file("drugbank/drug_links.csv")
	drug_targ = open_csv_file("drugbank/all_target_ids_with_known_action.csv")

	'''
	#Cmap list of drugs
	cmap_drugs = cmap_drug_list(cmap_file)
	#Drugbank list of drugs
	bank_drugs, drugbank_names = drugbank_list(drug_id_file)
	#Common Drug List between CMAP and DRUGBANK
	common_drugs_ls = common_drugs(cmap_drugs,drugbank_names,bank_drugs)
	geneids_drugs = drug_targ_dic(drug_targ)
	# Common Protein Targets
	geneid_target = prot_target_comm(geneids_drugs, common_drugs_ls)
	geneid_target_uniq = list(set(geneid_target))

	#for targets in geneid_target_uniq:
	#	print targets
	'''
	


if __name__ == "__main__":
        main()
