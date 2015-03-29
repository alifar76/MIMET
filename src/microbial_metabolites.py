"""
MIMET v0.0.1: MIcrobial METabolites
Copyright 2014 Ali A. Faruqi
This pipeline predicts metabolites from output produced through PICRUSt software (http://picrust.github.io/picrust/)
Based on KEGG database flatfiles, it produces a standard metabolites by sample table 
that can be used for further analyses.
"""


import argparse
import re
import numpy as np
import os
from collections import Counter
from operator import mul
from datetime import datetime


def comp_map(filename,kegg_compound):
	""" Create a file matching compound IDs (C) with their D and G synonyms """
	infile = open(kegg_compound, 'rU')
	outfile = open(filename, "w")
	for line in infile:
		if line.startswith("ENTRY"):
			g_id = []
			d_id = []
			while 1:
				next_line = infile.next().strip()
				if next_line.startswith("REMARK"):
					comp_search2 = re.compile(r"\bG\d\d\d\d\d\b").findall(next_line)
					comp_search3 = re.compile(r"\bD\d\d\d\d\d\b").findall(next_line)
					if comp_search2:
						g_id = comp_search2
					if comp_search3:
						d_id = comp_search3
				if next_line.startswith('///'):
					break
			comp_search = re.compile(r"\b(C\d\d\d\d\d)\b").search(line)
			if len(g_id) != 0:
				outfile.write(comp_search.group(1)+"\t"+'\t'.join(g_id)+"\n")
			if len(d_id) != 0:
				outfile.write(comp_search.group(1)+"\t"+'\t'.join(d_id)+"\n")
	outfile.close()
	return outfile


def replace_all(text, dic):
	for i, j in dic.iteritems():
		text = text.replace(i, j)
	return text


def equation_creator(mapfile,reactionfile,outputfile):
	id_list = []
	infile = open(mapfile, 'rU')				# Open comp_id_map.txt file (created with compMap function)
	for line in infile:
		spline = line.strip().split("\t")
		comp_search2 = re.compile(r"\bG\d\d\d\d\d\b").search(line)
		if comp_search2:
			for x in spline[1:]:
				comp_name = {}
				comp_name[x] = spline[0]
				id_list.append(comp_name)
	# id_list is a list of dictionaries. The dictionaries have G-ids as keys and C-ids as values derived from "comp_id_map.txt" file		
	g_ids = []
	for x in id_list:
		for y in x.keys():
			g_ids.append(y)						# All G-ids from id_list are appended to g_ids list
	gids = []
	for x in g_ids:
		if g_ids.count(x) > 1:
			gids.append(x)						# Append those G-ids to gids list that have count greater than 1
	multi_g_id = list(set(gids))				# Create a list of unique G-ids having count greater than 1 i.e., G-ids map to more than 1 C-ids
	conflict_list = []
	for m in multi_g_id:
		for x in id_list:
			try:
				if x[m]:
					conflict_list.append(x)		# Add dictionaries from id_list that have same G-ids corresponding to different C-ids
			except KeyError:
				pass
	final_list = []
	for x in id_list:
		if x not in conflict_list:
			final_list.append(x)				
	# final_list G-ids/C-ids dictionaries from id_list that do not have same G-ids mapping to different C-ids
	super_dict = {}
	for d in final_list:
		for k, v in d.iteritems():
			super_dict[k] = v					# Store G-ids as keys and C-ids as values. G-ids that have multiple C-ids have already been filtered.
	new_dict = {}
	for m in conflict_list:					
		new_dict[m[m.keys()[0]]] = m.keys()[0]	# Store C-ids as keys and G-ids as values from conflict_list in new_dict
	# Length of new_dict and conflict_list is same. It's 15.
	infile = open(reactionfile, 'rU')
	equations_store = []
	outfile = open(outputfile, "w")
	for line in infile:
		if line.startswith("EQUATION"):
			comp_search = re.compile(r"\bG\d\d\d\d\d\b").findall(line)	# Find all G-ids
			if comp_search:
				line_replace = replace_all(line,super_dict)			# Replace all G-ids with their corresponding C-ids
				comp_search2 = re.compile(r"\bC\d\d\d\d\d\b").findall(line_replace)
				if comp_search2:
					equation_conf = replace_all(line_replace,new_dict)		# Replace all already replaced equation's C-ids with G-ids that had multiple C-ids mapped
					equations_store.append(equation_conf.strip().split("EQUATION")[1].strip())
				else:
					equations_store.append(line_replace.strip().split("EQUATION")[1].strip())
			else:
				line_replace2 = replace_all(line,new_dict)
				if line_replace2 != line: 
					equations_store.append(line_replace2.split("EQUATION")[1].strip())
				else:
					equations_store.append(line.strip().split("EQUATION")[1].strip())
			ko_list = []	
			while 1:
				next_line = infile.next().strip()
				ko_search = re.compile(r"\bK\d\d\d\d\d\b").findall(next_line)		# Regex selects only KO
				if ko_search:
					ko_list.append(ko_search[0])
				if next_line.startswith('///'):
					break
			if len(ko_list) != 0:
				outfile.write(equations_store[-1]+"\t"+'\t'.join(ko_list)+"\n")
			else:
				outfile.write(equations_store[-1]+"\t"+"None"+"\n")
	outfile.close()
	return


def kos_and_comps(picrustfile,reactmapfile):
	""" Creates list of total KOs in PICRUSt and corresponding total metabolites involved with those KOs """
	kos = []
	infile = open(picrustfile,'rU')
	for line in infile:
		if not line.startswith("#"):
			kos.append(line.strip().split("\t")[0])
	comp_list = []
	reaction_ko_list = []
	for x in kos:
		reactions = []
		infile = open(reactmapfile, 'rU')
		for line in infile:
			if x in line.strip():
				spline = line.strip().split("\t")
				reactions.append(spline[0])
		if len(reactions) != 0:
			reaction_ko_list.append(x)
			for eqs in reactions:
				for react_prod in eqs.strip().split("<=>"):
					for compounds in react_prod.strip().split("+"):
						moles = re.compile("C\d\d\d\d\d|G\d\d\d\d\d").split(compounds.strip())
						if moles[0] != '':
							comp_id = re.compile(r"(C\d\d\d\d\d)|(G\d\d\d\d\d)").search(compounds.strip())
							if comp_id:
								comp_list.append(comp_id.group(1))
						if moles[0] == '':
							comp_list.append(compounds.strip())
	compounds_total = list(set(comp_list))
	return compounds_total, reaction_ko_list




def ko_compound_matrix(total_list,filename):
	outfile = open(filename,"w")
	outfile.write("Metabolite"+"\t"+"\t".join(total_list[0])+"\n")
	""" Very import piece of code for calculations of metabolites vs. KOs """
	for x in total_list[1]:
		reactions = []
		final_list = []
		compound_entries = ['0']*len(total_list[0])
		infile = open("KO_Reactions_Map.txt", 'rU')
		for line in infile:
			if x in line.strip():
				spline = line.strip().split("\t")
				reactions.append(spline[0])
		for reacts in reactions:
			for react_prod in reacts.strip().split("<=>"):
				for compounds in react_prod.strip().split("+"):
					moles = re.compile("C\d\d\d\d\d|G\d\d\d\d\d").split(compounds.strip())
					if moles[0] != '':
						comp_id = re.compile(r"(C\d\d\d\d\d)|(G\d\d\d\d\d)").search(compounds.strip())
						if comp_id:
							try:
								final_list.append([(comp_id.group(1),int(moles[0]))])
							except ValueError:
								#print moles[0]+"\t"+x+"\t"+reacts
								final_list.append([(comp_id.group(1),int('1'))])		# Set value 1 for all n
					if moles[0] == '':
						final_list.append([(compounds.strip(),int('1'))])
		main_dict = dict(sum((Counter(dict(x)) for x in final_list),Counter()))
		for comps in main_dict.keys():
			compound_entries[total_list[0].index(comps)] = main_dict[comps]
		outfile.write(x+"\t"+"\t".join(str(e) for e in compound_entries)+"\n")
	outfile.close()
	return outfile



def matrix_multiplication(compmatrix,komatrix,outputname):
	""" Multiply compound and KO matrices """
	metab_frame = []
	infile = open(compmatrix, 'rU')
	for line in infile:
		spline = line.strip().split("\t")
		metab_frame.append(spline)
	picrust_frame = []
	infile = open(komatrix, 'rU')
	for line in infile:
		if not line.startswith("# Constructed"):
			spline = line.strip().split("\t")
			picrust_frame.append(spline)
	picrust_subset = []						# Store those KOs from PICRUSt output that map to KOs of compound matrix
	picrust_subset.append(picrust_frame[0][:-1])
	for kos in metab_frame[1:]:
 		for kos2 in picrust_frame[1:]:
 			if kos2[0] == kos[0]:
 				picrust_subset.append(kos2[:-1])		# Ensuring that arrangement of KOs of PICRUSt data is same as KOs of compounds data
	picrust_values_data = []
	for m in picrust_subset[1:]:
		picrust_values_data.append([float(j) for j in m[1:]]) 
	picrust_array = np.array(picrust_values_data)
	outfile = open(outputname, "w")
	outfile.write("Metabolite"+"\t"+'\t'.join(picrust_subset[0][1:])+"\n")
	for comp_id in metab_frame[0][1:]:
		index_comp = metab_frame[0].index(comp_id)
		compound_col = [float(m) for m in zip(*metab_frame)[index_comp][1:]]
		for_multiply = np.repeat(np.array(compound_col),np.shape(picrust_array)[1])
		final_mat = np.reshape(for_multiply,np.shape(picrust_array))
		outfile.write(comp_id+"\t"+'\t'.join(str(e) for e in list(sum(final_mat*picrust_array)))+"\n")
	outfile.close()
	return outfile



def add_name_compound(routput,keggcomp,outputname):
	""" Add names to compound-ids """
	compounds = {}
	infile = open(routput,'rU')
	for line in infile:
		if not line.startswith("Metabolite"):
			spline = line.strip().split("\t")
			nm = re.compile(r"n|m").search(spline[0])
			if nm:
				comp_id = re.compile(r"(C\d\d\d\d\d)|(G\d\d\d\d\d)").search(spline[0])
				if comp_id:
					compounds[spline[0]] = comp_id.group()
			else:
				compounds[spline[0]] = spline[0]
	name_dict = {}
	infile = open(keggcomp, 'rU')
	for line in infile:
		if line.startswith("ENTRY"):
			allines = []
			compline = line.strip().split("ENTRY")[1].strip()
			if compline.strip().split()[0] in compounds.values():	#['C18233']
				while 1:
					next_line = infile.next().strip()
					allines.append(next_line)
					if (next_line.startswith("FORMULA")) or (next_line.startswith("REACTION")) or (next_line.startswith("SEQUENCE")) or (next_line.startswith("COMMENT")) or (next_line.startswith("REMARK")):
						break
				names = []
				for k in allines[:-1]:
					if len(k.split("NAME")) == 2:
						names.append(k.split("NAME")[1].strip())
					else:
						names.append(k.split("NAME")[0])
				name_dict[compline.strip().split()[0]] = str(names)
	for ids in compounds.values():
		if ids not in name_dict.keys():
			name_dict[ids] = "['None']"
	outfile = open(outputname, "w")
	infile = open(routput,'rU')
	for line in infile:
		if line.startswith("Metabolite"):
			spline = line.strip().split("\t")
			outfile.write("#OTU ID"+"\t"+'\t'.join(spline[1:])+"\n")
		else:
			spline = line.strip().split("\t")
			outfile.write(spline[0]+" "+name_dict[compounds[spline[0]]]+"\t"+'\t'.join(spline[1:])+"\n")		
	outfile.close()
	return outfile


def main_metabolite(picrustout,outname,react,comps):
	""" Function to call all other functions """
	comp_map("comp_id_map.txt",comps)
	equation_creator("comp_id_map.txt",react,"KO_Reactions_Map.txt")
	total_list = kos_and_comps(picrustout,"KO_Reactions_Map.txt")		# Replace high_lgg_metabolite_pipeline.txt with PICRUSt result
	ko_compound_matrix(total_list,"ko_compound_matrix.txt")
	matrix_multiplication("ko_compound_matrix.txt",picrustout,"Metabolite_PICRUSt.txt")		# Replace high_lgg_metabolite_pipeline.txt with PICRUSt result
	add_name_compound("Metabolite_PICRUSt.txt",comps,outname)				# Replace high_lgg_metabolite_pipeline_predicted.spf with name of final output file
	os.system('rm comp_id_map.txt KO_Reactions_Map.txt ko_compound_matrix.txt Metabolite_PICRUSt.txt')
	return


if __name__ == '__main__':
	startTime = datetime.now()
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
	description='Microbial Metabolomics',
	epilog='''
An example to run the pipeline:
python microbial_metabolites.py -i all_samples_metabolite_pipeline.tab -o all_samples_metabolite_pipeline_predicted.txt -r reaction -c compound
	''')
	parser.add_argument('-i', metavar='Input file', nargs=1, help='Name of the PICRUSt output file that\
						serves as input the for metabolomics pipeline',required=True)
	parser.add_argument('-o', metavar='Output file', nargs=1, help='Name of the output file produced\
						by the metabolomics pipeline',required=True)
	parser.add_argument('-r', metavar='Reaction',nargs=1, help='Full path to the KEGG reaction\
						database file',required=True)
	parser.add_argument('-c', metavar='Compound', nargs=1, help='Full path to the KEGG compound\
						database file',required=True)
	args = parser.parse_args()
	main_metabolite(args.i[0],args.o[0],args.r[0],args.c[0])
	print "\n"+"Task Completed! Completion time: "+ str(datetime.now()-startTime)