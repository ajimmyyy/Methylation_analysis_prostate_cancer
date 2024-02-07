from goatools.semantic import TermCounts, common_parent_go_ids, min_branch_length, semantic_distance, semantic_similarity
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools.base import get_godag
from collections import defaultdict
from goatools.semantic import deepest_common_ancestor
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score
from scipy.cluster.hierarchy import dendrogram 
import scipy.cluster.hierarchy as shc
import numpy as np
import pandas as pd
from collections import namedtuple
from sklearn import metrics
from scipy.spatial.distance import squareform
import sys
import random
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
from goatools.semantic import TermCounts, get_info_content
from goatools.associations import read_associations
from pathlib import Path
import os

file_dir = Path(__file__).resolve().parent / "_resources"
resource_dir = file_dir / "go-basic.obo"
if not os.path.exists(str(file_dir)):
    os.makedirs(str(file_dir))
fin_dag = download_go_basic_obo(str(resource_dir))
go = GODag(fin_dag, optional_attrs={'relationship', 'replaced_by'}, load_obsolete=True)

#Relationship Semantic Contribution 
is_a = 0.8
part_of = 0.6
regulates = 0.7
negatively_regulates = 0.7
positively_regulates = 0.7
reg = 'reg0.7'
Gene = namedtuple('Gene', 'GeneName GOAnnotations')
c = 0.67

### New Functions 
similarities_down = {}

def all_common_parent_go_ids(goids, godag):
	'''
		This function finds the common ancestors in the GO
		tree of the list of goids in the input.
	'''
	# Find candidates from first
	rec = godag[goids[0]]
	candidates = rec.get_all_upper()
	candidates.update({goids[0]})

	# Find intersection with second to nth goid
	for goid in goids[1:]:
		rec = godag[goid]
		parents = rec.get_all_upper()
		parents.update({goid})

		# Find the intersection with the candidates, and update.
		candidates.intersection_update(parents)

	return candidates
def lowest_common_ancestor(goterms, godag):
	'''
		This function gets the nearest common ancestor
		using the above function.
		Only returns single most specific - assumes unique exists.
	'''
	# Take the element at maximum depth.
	return max(all_common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)
def min_branch_length(go_id1, go_id2, godag, branch_dist):
	'''
		Finds the minimum branch length between two terms in the GO DAG.
	'''
	# First get the deepest common ancestor
	goterm1 = godag[go_id1]
	goterm2 = godag[go_id2]
	if goterm1.namespace == goterm2.namespace:
		dca = lowest_common_ancestor([go_id1, go_id2], godag)
		# Then get the distance from the DCA to each term
		dca_depth = godag[dca].depth
		depth1 = goterm1.depth - dca_depth
		depth2 = goterm2.depth - dca_depth

		# Return the total distance - i.e., to the deepest common ancestor and back.
		return depth1 + depth2

	elif branch_dist is not None:
		return goterm1.depth + goterm2.depth + branch_dist
def semantic_distance(go_id1, go_id2, godag, branch_dist=None):
	'''
		Finds the semantic distance (minimum number of connecting branches)
		between two GO terms.
	'''
	return min_branch_length(go_id1, go_id2, godag, branch_dist)
def semantic_similarity(go_id1, go_id2, godag, branch_dist=None):
	'''
		Finds the semantic similarity (inverse of the semantic distance)
		between two GO terms.
	'''
	dist = semantic_distance(go_id1, go_id2, godag, branch_dist)
	if dist is not None:
		return 1.0 / float(dist)	
		
def semantic_value_resnik_lin_IEA_MF():
	associations = read_associations("association_IEA_MF.txt")	
	termcounts = TermCounts(go, associations)
	#infocontent = get_info_content(go_id, termcounts)
	#print('Information content ({}) = {}'.format(go_id, infocontent))
	return termcounts 
	
def semantic_value_resnik_lin_nonIEA_MF():
	associations = read_associations("association_nonIEA_MF.txt")	
	termcounts = TermCounts(go, associations)
	return termcounts 
	

def all_paths_to_top(term, godag):
	# inputs: term_id and Go dag with 'relationship' as optional attributes
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper()
			for parent in parents:
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths

		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
def relationships_in_all_paths_to_top(go_id):
	all_all_paths = all_paths_to_top(go_id, go)
	for index, path in enumerate(all_all_paths):
		print("index = " + str(index))
		path = list(reversed(path))
		#path[::-1] #shallow copying a list in reverse order
		for idx, term in enumerate(path):
			if idx == 0:	#which is on idx = 0 
				print ('Term Itself')
			if term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150':
				print ('Root')
			print (term.item_id)
			if idx < len(path)-1:
				if term.relationship != {}:
					#print(term.relationship.keys())
					if 'part_of' in term.relationship:
						if path[idx+1] in term.relationship['part_of']:
							print ('part_of')
						else: 
							print ('is_a')		
					elif 'regulates' in term.relationship:
						if path[idx+1] in term.relationship['regulates']:
							print ('regulates')
						else: 
							print ('is_a')
					elif 'negatively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['negatively_regulates']:
							print ('negatively_regulates')	
						else: 
							print ('is_a')	
					elif 'positively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['positively_regulates']:
							print ('positively_regulates')
						else: 
							print ('is_a')
				else: 
					print ('is_a')
def all_paths_to_top_wang(term, godag, optional_relationships):
	# inputs: term_id and Go dag with 'relationship' as optional attributes
		""" Returns all possible paths to the root node"""
		if term not in godag:
			sys.stderr.write("Term %s not found!\n" % term)
			return
		def _all_paths_to_top_recursive(rec):
			if rec.level == 0:
				return [[rec]]
			paths = []
			parents = rec.get_goterms_upper_rels(optional_relationships)
			#parents = rec.get_goterms_upper()
			for parent in parents:
				#parent = go[parent1]
				top_paths = _all_paths_to_top_recursive(parent)
				for top_path in top_paths:
					top_path.append(rec)
					paths.append(top_path)
			return paths
		go_term = godag[term]
		return _all_paths_to_top_recursive(go_term)
def relationships_in_all_paths_to_top_wang(go_id, go, optional_relationships):
	all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
	for index, path in enumerate(all_all_paths):
		print("index = " + str(index))
		path = list(reversed(path))
		#path[::-1] #shallow copying a list in reverse order
		for idx, term in enumerate(path):
			if idx == 0:	#which is on idx = 0 
				print ('Term Itself')
			if term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150':
				print ('Root')
			print (term.item_id)
			if idx < len(path)-1:
				if term.relationship != {}:
					#print(term.relationship.keys())
					if 'part_of' in term.relationship:
						if path[idx+1] in term.relationship['part_of']:
							print ('part_of')
						else: 
							print ('is_a')		
					elif 'regulates' in term.relationship:
						if path[idx+1] in term.relationship['regulates']:
							print ('regulates')
						else: 
							print ('is_a')
					elif 'negatively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['negatively_regulates']:
							print ('negatively_regulates')	
						else: 
							print ('is_a')	
					elif 'positively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['positively_regulates']:
							print ('positively_regulates')
						else: 
							print ('is_a')
				else: 
					print ('is_a')				
def Semantic_Value(go_id, method):
	''' input: goterm_id 
		returns all of the weighted (all relatinships in go-basic) paths to root 
		#relationship types are global variables with appropriate weights
	'''
	if method == 'wang':
		# calculates all paths to top (with all relationships)
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
#								print ('part_of')
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else: 
#								print ('is_a')
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						else: 
#							print ('is_a')
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
					else: 
#						print ('is_a')
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
		return final_values(S_values,'max')
	elif method == 'GOGO':
		optional_relationships = {'part_of'}
		all_all_paths = all_paths_to_top_wang(go_id, go, optional_relationships)
		S_values = list()
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
				if idx < len(path)-1:
					if term.relationship != {}:
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
#								print ('part_of')
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + part_of
						# 		print(path[idx+1].item_id)
# 								print(len(go[path[idx+1].item_id].children))
# 								print(weight)
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight))
							else: 
#								print ('is_a')
								weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
							# 	print(path[idx+1].item_id)
# 								print(len(go[path[idx+1].item_id].children))
# 								print(weight)
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
						else: 
#							print ('is_a')
							weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
				# 			print(path[idx+1].item_id)
# 							print(len(go[path[idx+1].item_id].children))
# 							print(weight)
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
					else: 
#						print ('is_a')
						weight = (1 / (c + len(go[path[idx+1].item_id].children))) + is_a
				# 		print(path[idx+1].item_id)
# 						print(len(go[path[idx+1].item_id].children))
# 						print(weight)
						S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * weight ))
		return final_values(S_values,'max')
	#Modified Wang's Measure - Almost Same as Semantic Value; only Difference is the weight of root as ancestor = 0 and different realtionships
	else:
		all_all_paths = all_paths_to_top(go_id, go)
		S_values = list()
		#print(go_id)
		for index, path in enumerate(all_all_paths):
			S_values.append([])
			#print("index = " + str(index))
			path.reverse()
			for idx, term in enumerate(path):
				#print (term.item_id)
				# Semantic Value of the term itself is always 1 
				if idx == 0:	#which is on idx = 0 
					S_values[index].append((go_id, 1))
					#print (term.item_id)
				if idx < len(path)-1 and path[idx+1].item_id != 'GO:0003674' and path[idx+1].item_id != 'GO:0005575' and path[idx+1].item_id != 'GO:0008150':
					if term.relationship != {}: 
						#print(term.relationship.keys())
						if 'part_of' in term.relationship:
							if path[idx+1] in term.relationship['part_of']:
								# print ('part_of')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
							else: 
								# print ('is_a')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'regulates' in term.relationship:
							if path[idx+1] in term.relationship['regulates']:
								# print ('regulates')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
							else: 
								# print ('is_a')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'negatively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['negatively_regulates']:
								# print ('negatively_regulates')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
							else: 
								# print ('is_a')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
						elif 'positively_regulates' in term.relationship:
							if path[idx+1] in term.relationship['positively_regulates']:
								# print ('positively_regulates')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
							else: 
								# print ('is_a')
#								print (path[idx+1].item_id)
								S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
						else: 
							# print ('is_a')
#							print (path[idx+1].item_id)
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					else: 
						# print ('is_a')
#						print (path[idx+1].item_id)
						S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
				if (term.item_id == 'GO:0003674' or term.item_id == 'GO:0005575' or term.item_id == 'GO:0008150'):
						# print('is-a')
#						print(term.item_id)
						S_values[index].append((term.item_id, 0))
		if method == 'Baseline':
			return final_values(S_values, 'max')
		if method == 'Baseline_LCA_max' or method == 'Baseline_LCA_avg' :
			svaluesfinal = final_values(S_values, 'max')
			# Replace the max or min s-values in all paths to assign only one value to each node, consistent in all paths
			S_values_Modified = list()
			for index, path in enumerate(S_values):
			#	print("index: ", index)
			#	print('path: ', path)
				S_values_Modified.append([])
				for idx, term1 in enumerate(path):
			#		print(list(value for (term, value) in svaluesfinal if term == term1[0]))
					ind = [value for (term, value) in svaluesfinal if term == term1[0]]
					S_values_Modified[index].append((path[idx][0],ind[0]))
			SumOfNodesOnEachPath = list()
			for index, path in enumerate(S_values_Modified):
			#	print("Modified index and path: ", index, path)
				SumOfNodesOnEachPath.append(sum(x[1] for x in path))
			#print(SumOfNodesOnEachPath)
			#if isMax == 'max':
			maxPath = SumOfNodesOnEachPath.index(max(SumOfNodesOnEachPath))
			return S_values_Modified[maxPath], SumOfNodesOnEachPath[maxPath]

def final_values(S_values, isMax):

	''' helper function to assign the max of the weights assigned to each term'''
	#S_values = sorted(S_values, key=lambda x: x[0])
	unique_terms_s_values = []
	for path in S_values:
		for term in path:
			unique_terms_s_values.append(term)
			#print(term)
	unique_terms_s_values = sorted(unique_terms_s_values, key=lambda x: x[0])
	_s_values = {}
	for y, x in unique_terms_s_values: 
		if y in _s_values: 
			_s_values[y].append((y,x)) 
		else: 
			_s_values[y] = [(y, x)]
	final_s_values = []
	if(isMax == 'max'):
		for node in _s_values:
			final_s_values.append(max(_s_values[node]))
	elif (isMax == 'min'):
		for node in _s_values:
			final_s_values.append(min(_s_values[node]))
	return final_s_values
def intersection(lst1, lst2): 
	#[[x for x in sublist if x in list1] for sublist in list2]
	'''Helper Function to find intersecting terms from the two input lists of (term, s_Value)'''
	da = {v:k for v, k in lst1}
	db = {v:k for v, k in lst2} 
	# check to see if intersecting GO_ids are correct
	#[(k, da[k],db[k]) for k in da.keys() & db.keys()]
	return [(da[k],db[k]) for k in da.keys() & db.keys()]
## Downward Graph 				  
def common_children_go_ids(goids, godag):
	'''
		This function finds the common children in the GO
		tree of the list of goids in the input.
	'''
	# Find candidates from first
	rec = godag[goids[0]]
	#candidates = rec.get_all_children()
	candidates = rec.get_all_lower()
	candidates.update({goids[0]})
	# Find intersection with second to nth goid
	for goid in goids[1:]:
		rec = godag[goid]
		#children = rec.get_all_children()
		children = rec.get_all_lower()
		children.update({goid})
		# Find the intersection with the candidates, and update.
		candidates.intersection_update(children)
	return candidates
def highest_common_descendant(goterms, godag):
	'''
		This function gets the nearest common descendant
		using the above function.
		Only returns single most specific - assumes unique exists.
	'''
	# Take the element at minimum depth.
	common_children = common_children_go_ids(goterms, godag)
	if len(common_children) != 0:
		return min(common_children, key=lambda t: godag[t].depth)
	else:
		return 0
def all_paths_to_bottom(term, godag, x):
		# inputs: term_id and Go dag with 'relationship' as optional attributes
	""" Returns all possible paths to the root node"""
	if term not in godag:
		sys.stderr.write("Term %s not found!\n" % term)
		return

	def _all_paths_to_bottom_recursive(rec):

		if rec.depth == godag[term].depth+x:
			return [[rec]]
		else:
			paths = []
			children = rec.get_goterms_lower()
			for child in children:
				bottom_paths = _all_paths_to_bottom_recursive(child)
				for bottom_path in bottom_paths:
					bottom_path.append(rec)
					paths.append(bottom_path)
			return paths
	go_term = godag[term]
	return _all_paths_to_bottom_recursive(go_term)
def Downward_Semantic_Value(go_id, go, x):
	''' input: goterm_id 
		returns all of the weighted nodes in path to the Goterm from the children at the xth level below
		#relationship types are global variables with appropriate weights
	'''
	all_all_paths = all_paths_to_bottom(go_id, go, x)
	#print(len(all_all_paths))
	S_values = list()
	for index, path in enumerate(all_all_paths):
		S_values.append([])
		#print("index = " + str(index))
		path.reverse()
		for idx, term in enumerate(path):
			#print (term.item_id)
			if idx == 0:	#which is on idx = 0 
				S_values[index].append((go_id, 1))
			if idx < len(path)-1:
				if term.relationship != {}: 
					#print(term.relationship.keys())
					if 'part_of' in term.relationship:
						if path[idx+1] in term.relationship['part_of']:
							#print ('part_of')
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * part_of))
						else: 
							#print ('is_a')
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'regulates' in term.relationship:
						if path[idx+1] in term.relationship['regulates']:
							#print ('regulates')
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * regulates))
						else: 
							#print ('is_a')
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'negatively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['negatively_regulates']:
							#print ('negatively_regulates')
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * negatively_regulates))
						else: 
							#print ('is_a')
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
					elif 'positively_regulates' in term.relationship:
						if path[idx+1] in term.relationship['positively_regulates']:
							#print ('positively_regulates')
							S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * positively_regulates))
						else: 
							#print ('is_a')
							S_values[index].append((path[idx+1].item_id, S_values[index][idx][1] * is_a))
				else: 
					#print ('is_a')
					S_values[index].append((path[idx+1].item_id,S_values[index][idx][1] * is_a))
	#print(S_values)
	return final_values(S_values, 'max')
### Calculating Similarity 
# def Similarity_of_Two_GOTerms(go_id1, go_id2, go, method):
# #	if go_id1 == go_id2:
# #			return 1.0000
# 	if method == 'Baseline_LCA_max' or method == 'Baseline_LCA_min':
# 		lca = lowest_common_ancestor((go_id1, go_id2), go)
# 		sim = Semantic_Value(lca, go, method)
# 		return sim[1]
# 	if method == 'Baseline_LCA_avg':
# 		lca = lowest_common_ancestor((go_id1, go_id2), go)
# 		#print(lca)
# 		sim_lca = Semantic_Value(lca, go, method)
# 		#print(sim_lca[1])
# 		sim1 = Semantic_Value(go_id1, go, method)
# 		sim2 = Semantic_Value(go_id2, go, method)
# 		
# 		#print(sim1, sim2)
# 		# sim1[1] and sim2[1] are the sums of all the nodes on the path with the max s-values. 
# 		avg_div = (sim1[1] + sim2[1])/2
# 		return (sim_lca[1]/avg_div)
# 	
# 	elif method == 'GOntoSim':
# 		hcd = highest_common_descendant((go_id1, go_id2), go)
# 		if hcd != 0:
# 			hcd_depth = go[hcd].depth
# 			go1_depth = go[go_id1].depth
# 			go2_depth = go[go_id2].depth
# 			x = hcd_depth - go1_depth
# 			y = hcd_depth - go2_depth
# 			sv_a = Downward_Semantic_Value(go_id1, go, x)
# 			sv_b = Downward_Semantic_Value(go_id2, go, y)
# 			intersecting_terms = intersection(sv_a,sv_b)
# 			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
# 			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
# 			sim_down = (numerator/denominator)
# 		else:
# 			sim_down = 0
# 		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline_LCA_avg')
# 		sim = (sim_down*0.5) + sim_lca_upper
# 		return sim
# 		
# 	elif method == 'wang' or method == 'Baseline' or method == 'GOGO':
# 		sv_a = Semantic_Value(go_id1, go, method)
# 		sv_b = Semantic_Value(go_id2, go,  method)
# 		#print(sv_a[1], sv_b[1])
# 		intersecting_terms = intersection(sv_a,sv_b)
# 		#print([x for t in intersecting_terms for x in t])
# 		numerator = sum([x for t in intersecting_terms for x in t])
# 		#print(numerator)
# 		denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
# 		#print (denominator)
# 		Similarity = (numerator/denominator)
# 		return Similarity
# 	elif method == 'Baseline_Desc':
# 		hcd = highest_common_descendant((go_id1, go_id2), go)
# 		if hcd != 0:
# 			hcd_depth = go[hcd].depth
# 			go1_depth = go[go_id1].depth
# 			go2_depth = go[go_id2].depth
# 			x = hcd_depth - go1_depth
# 			y = hcd_depth - go2_depth
# 			sv_a = Downward_Semantic_Value(go_id1, go, x)
# 			sv_b = Downward_Semantic_Value(go_id2, go, y)
# 			intersecting_terms = intersection(sv_a,sv_b)
# 			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
# 			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
# 			sim_down = (numerator/denominator)
# 		else:
# 			sim_down = 0
# 		sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline')
# 		sim = (sim_down*0.5) + sim_upper
# 		return sim
# 	elif method == 'Baseline_Desc_only':
# 		hcd = highest_common_descendant((go_id1, go_id2), go)
# 		if hcd != 0:
# 			hcd_depth = go[hcd].depth
# 			go1_depth = go[go_id1].depth
# 			go2_depth = go[go_id2].depth
# 			x = hcd_depth - go1_depth
# 			y = hcd_depth - go2_depth
# 			sv_a = Downward_Semantic_Value(go_id1, go, x)
# 			sv_b = Downward_Semantic_Value(go_id2, go, y)
# 			intersecting_terms = intersection(sv_a,sv_b)
# 			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
# 			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
# 			sim_down = (numerator/denominator)
# 		else:
# 			sim_down = 0
# 		#sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline')
# 		#sim = (sim_down*0.5) + sim_upper
# 		return sim_down
# 			
def Similarity_of_Two_GOTerms(go_id1, go_id2, go, method, S_values):

	if method == 'Baseline_LCA_max' or method == 'Baseline_LCA_min':
		lca = lowest_common_ancestor((go_id1, go_id2), go)
		sim = Semantic_Value(lca, method)
		return sim[1]
	if method == 'Baseline_LCA_avg':
		lca = lowest_common_ancestor((go_id1, go_id2), go)
		#print(lca)
#		sim_lca = Semantic_Value(lca, go, method)
		#print(sim_lca[1])
#		sim1 = Semantic_Value(go_id1, go, method)
#		sim2 = Semantic_Value(go_id2, go, method)

		sim_lca = S_values.get(lca, 'none')
		#print(sim_lca[1])
		if sim_lca == 'none':
			sim_lca = Semantic_Value(lca, method)
			S_values[lca] = sim_lca
			print(lca)
			
		sim1 = S_values.get(go_id1)
		sim2 = S_values.get(go_id2)

		#print(sim1[1], sim2[1])
		# sim1[1] and sim2[1] are the sums of all the nodes on the path with the max s-values. 
		avg_div = (sim1[1] + sim2[1])/2
		return (sim_lca[1]/avg_div)	
	elif method == 'GOntoSim':
		sim = similarities_down.get((go_id1, go_id2), 'none')
		if sim != 'none':
			return sim
		else:	
			if go_id1 != 'GO:0003674' and go_id1 !='GO:0005575' and go_id1 != 'GO:0008150' and  go_id2 != 'GO:0003674' and go_id2 !='GO:0005575' and go_id2 != 'GO:0008150':
				hcd = highest_common_descendant((go_id1, go_id2), go)
				if hcd != 0:
					hcd_depth = go[hcd].depth
					go1_depth = go[go_id1].depth
					go2_depth = go[go_id2].depth
					x = hcd_depth - go1_depth
					y = hcd_depth - go2_depth
					sv_a = Downward_Semantic_Value(go_id1, go, x)
					sv_b = Downward_Semantic_Value(go_id2, go, y)
					intersecting_terms = intersection(sv_a,sv_b)
					numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
					denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
					# print(go_id1)
# 					print(go_id2)
# 					print(numerator)
# 					print(denominator)
					if go_id1 != go_id2:
						denominator = denominator - numerator
# 					print(denominator)
					sim_down = (numerator/denominator)
				else:
					sim_down = 0
			else: 
				sim_down = 0
			sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline_LCA_avg', S_values)
			sim = (sim_down*0.5) + (sim_upper*0.5)
			#sim = (sim_down*0.3) + (sim_upper*0.7)            
			similarities_down[(go_id1,go_id2)] = sim
			return sim
	elif method == 'wang' or method == 'GOGO':
		sv_a = S_values.get(go_id1)
		sv_b = S_values.get(go_id2)
		intersecting_terms = intersection(sv_a,sv_b)
		#print([x for t in intersecting_terms for x in t])
		numerator = sum([x for t in intersecting_terms for x in t])
		#print(numerator)
		denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
		#print (denominator)
		Similarity = (numerator/denominator)
		return Similarity
	elif method == 'Baseline':
		sv_a = S_values.get(go_id1)
		sv_b = S_values.get(go_id2)
		intersecting_terms = intersection(sv_a,sv_b)
		#print([x for t in intersecting_terms for x in t])
		numerator = sum([x for t in intersecting_terms for x in t])
		#print(numerator)
		denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
	# 	print(go_id1)
# 		print(go_id2)
# 		print(numerator)
# 		print(denominator)
#		denominator = denominator - numerator
		if go_id1 != go_id2:
			denominator = denominator - numerator
		Similarity = (numerator/denominator)
		return Similarity
	elif method == 'Baseline_Desc':
		sim = similarities_down.get((go_id1, go_id2), 'none')
		if sim != 'none':
			return sim
		else:	
			if go_id1 != 'GO:0003674' and go_id1 !='GO:0005575' and go_id1 != 'GO:0008150' and  go_id2 != 'GO:0003674' and go_id2 !='GO:0005575' and go_id2 != 'GO:0008150':
				hcd = highest_common_descendant((go_id1, go_id2), go)
				if hcd != 0:
					hcd_depth = go[hcd].depth
					go1_depth = go[go_id1].depth
					go2_depth = go[go_id2].depth
					x = hcd_depth - go1_depth
					y = hcd_depth - go2_depth
					sv_a = Downward_Semantic_Value(go_id1, go, x)
					sv_b = Downward_Semantic_Value(go_id2, go, y)
					intersecting_terms = intersection(sv_a,sv_b)
					numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
					denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
					if go_id1 != go_id2:
						denominator = denominator - numerator
 					
					sim_down = (numerator/denominator)
				else:
					sim_down = 0
			else: 
				sim_down = 0
			sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline', S_values)
			sim = (sim_down*0.5) + (sim_upper * 0.5)
			#sim = (sim_down*0.3) + (sim_upper * 0.7)
			
			# similarity with itself got to 1.5 which should not be > 1
			#sim = (sim_down*0.5) + (sim_upper)

			similarities_down[(go_id1,go_id2)] = sim
			return sim
	elif method == 'Baseline_Desc_only':
		sim = similarities_down.get((go_id1, go_id2), 'none')

		hcd = highest_common_descendant((go_id1, go_id2), go)
		if hcd != 0:
			hcd_depth = go[hcd].depth
			go1_depth = go[go_id1].depth
			go2_depth = go[go_id2].depth
			x = hcd_depth - go1_depth
			y = hcd_depth - go2_depth
			sv_a = Downward_Semantic_Value(go_id1, go, x)
			sv_b = Downward_Semantic_Value(go_id2, go, y)
			intersecting_terms = intersection(sv_a,sv_b)
			numerator = sum([x for t in intersecting_terms for x in t]) # Sum of common terms in all paths wrt to each of the 2 terms
			denominator = sum(x for y, x in sv_a) + sum(x for y, x in sv_b) #(where sv_a has 2 values for each term, second being the SV)
			sim_down = (numerator/denominator)
		else:
			sim_down = 0
		#sim_upper = Similarity_of_Two_GOTerms(go_id1,go_id2, go, 'Baseline')
		#sim = (sim_down*0.5) + sim_upper
		return sim_down	
	

		
def Similarity_of_Set_of_GOTerms(set1, set2, method, S_values):
	if set1 == [] or set2 == []:
		return 0

	Sim1 = []
	Sim2 = []
	#print(len(set1))
	#print(len(set2))
	for idx, goterm in enumerate(set1):
		# print ("=========", goterm)
		Sim1.append([])
		for goid in set2:
			#print (goterm , goid)
			Sim1[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method, S_values))))
	#print(Sim1)
	#print([y[0][0] for y in Sim1])
	#print([[y[1] for y in	x] for x in Sim1])
	#print([y[1] for y in Sim1[0]])
	
	#print(pd.DataFrame(data = [[y[2] for y in	x] for x in Sim1], index = [y[0][0] for y in Sim1], columns = [y[1] for y in Sim1[0]]))

	for idx, goterm in enumerate(set2):
		Sim2.append([])
		for goid in set1:
			#print (goterm , goid)
			Sim2[idx].append((goterm, goid,(Similarity_of_Two_GOTerms(goterm, goid, go, method, S_values))))
			
	#print(pd.DataFrame(data = [[y[2] for y in	x] for x in Sim2], index = [y[0][0] for y in Sim2], columns = [y[1] for y in Sim2[0]]))

	sem1 = []
	sem2 = []

	for index, goterm in enumerate(Sim1):
		sem1.append((max(Sim1[index], key=lambda x: x[2])))
	#print(sem1)
	for index, goterm in enumerate(Sim2):	
		sem2.append((max(Sim2[index], key=lambda x: x[2])))
	#print(sem2)
	#
	similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(set1) + len(set2))

	return round(similarity, 3)
	


def Similarity_of_Two_GOTerms_IC(go_id1, go_id2, go, method, termcounts):
	if method == 'Resnik':
		Similarity = resnik_sim(go_id1, go_id2, go, termcounts)
		return Similarity
	elif method == 'Lin':
		Similarity = lin_sim(go_id1, go_id2, go, termcounts)
		return  Similarity
def Similarity_of_Set_of_GOTerms_IC(set1, set2, method, termcounts):
	
	Sim1 = []
	Sim2 = []

	for idx, goterm in enumerate(set1):
		Sim1.append([])
		for goid in set2:
			Sim1[idx].append((goterm, goid,(Similarity_of_Two_GOTerms_IC(goterm, goid, go, method, termcounts))))
	
	for idx, goterm in enumerate(set2):
		Sim2.append([])
		for goid in set1:
			Sim2[idx].append((goterm, goid,(Similarity_of_Two_GOTerms_IC(goterm, goid, go, method, termcounts))))
			
	sem1 = []
	sem2 = []

	for index, goterm in enumerate(Sim1):	
		sem1.append((max(Sim1[index], key=lambda x: x[2])))
	#print(sem1)
	for index, goterm in enumerate(Sim2):	
		sem2.append((max(Sim2[index], key=lambda x: x[2])))
	#print(sem2)
	#
	similarity = (sum(x[2] for x in sem1)+ sum(x[2] for x in sem2) )/(len(set1) + len(set2))
	
	return round(similarity, 3)
	
		
## Evaluation: Clustering, AMI/ARI Scores

def Similarity_Matrix(genes, method, S_values):
	sim_matrix = []
	for idx,gene in enumerate(genes):
		print(gene)
		sim_matrix.append([(lambda x: Similarity_of_Set_of_GOTerms(x[1],gene[1], method, S_values))(x) for x in genes])
	return sim_matrix

def Similarity_Matrix_IC(genes, method, termcounts):
	sim_matrix = []
	for idx,gene in enumerate(genes):
		print(gene)
		sim_matrix.append([(lambda x: Similarity_of_Set_of_GOTerms_IC(x[1],gene[1], method, termcounts))(x) for x in genes])
	return sim_matrix


def Agglomerative_Clustering(pathway, Genes, n_clusters, method, S_values):
	# Similarity Matrix
	data = Similarity_Matrix(Genes, method, S_values)
	length = len(data)
	data1 = pd.DataFrame(data = data, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	#print('similarity matrisx: ')
	#print(data1)
	#writeSim = pathway + "_SimilarityMatrix.csv"
	#data1.to_csv(writeSim)
	data_matrix = []
	for row in data:
		data_matrix.append([(lambda x: 1-x)(x) for x in row])
	my_data = pd.DataFrame(data = data_matrix, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
		# print('distance matrix: ')
	# print( my_data)
	# writeDist = pathway + "_DistanceMatrix.csv"
	# my_data.to_csv(writeDist)	
	return Agglomerative(my_data,Genes, pathway, n_clusters)
	

def Agglomerative(data, Genes, pathway, n_clusters):
	model = AgglomerativeClustering(n_clusters, affinity='precomputed', linkage='complete').fit_predict(data)
	#GeneNames = [gene.GeneName for gene in Genes]
	return model.tolist()

def AgglomerativeClusteringIC(pathway, Genes, n_clusters, method, termcounts):
	# Similarity Matrix
	data = Similarity_Matrix_IC(Genes, method, termcounts)
	length = len(data)
	data1 = pd.DataFrame(data = data, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	#print('similarity matrisx: ')
	#print(data1)
	#writeSim = pathway + "_SimilarityMatrix.csv"
	#data1.to_csv(writeSim)
	data_matrix = []
	for row in data:
		data_matrix.append([(lambda x: 1-x)(x) for x in row])
	my_data = pd.DataFrame(data = data_matrix, index=[x[0] for x in Genes],columns=[x[0]for x in Genes])
	# print('distance matrix: ')
	# print( my_data)
	# writeDist = pathway + "_DistanceMatrix.csv"
	# my_data.to_csv(writeDist)		
	return Agglomerative(my_data,Genes, pathway, n_clusters)



def plot_dendrogram(model, dist,**kwargs):
	distance1 = dist.to_numpy() 
	distance = squareform(distance1)
	#print('Plot Dendogram')
	#print('distance ', distance)
	# Create linkage matrix and then plot the dendrogram
	linkage_matrix = shc.linkage(distance,'complete' )
	#print(linkage_matrix)
	# Plot the corresponding dendrogram
	dendrogram(linkage_matrix, leaf_rotation=90, **kwargs)


def ARI_Score(label_pred, labels_true):
	return adjusted_rand_score(label_pred, labels_true)

def AMI_Score(label_pred, labels_true):
	return adjusted_mutual_info_score(label_pred, labels_true, average_method='arithmetic')

def cont_matrix (y_true, y_pred):
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	# return purity
	return(contingency_matrix)

def purity_score(y_true, y_pred):
	# compute contingency matrix (also called confusion matrix)
	contingency_matrix = metrics.cluster.contingency_matrix(y_true, y_pred)
	# return purity
	#print(contingency_matrix)
	return np.sum(np.amax(contingency_matrix, axis=0)) / np.sum(contingency_matrix) 

def goterm_replaced_by(term):
	if go[term].replaced_by[:3] == 'GO:':
		return go[term].replaced_by
	else:
		return term

def main():

# 	if sys.argv[1] == 'GOterms':		
	if sys.argv[1] == 'wang':
		method = 'wang'
	elif sys.argv[1] == 'baseline':
		method = 'Baseline'
	elif sys.argv[1] == 'lca':
		method = 'Baseline_LCA_avg'
	elif sys.argv[1] == 'baselineDesc':
		method = 'Baseline_Desc'
	elif sys.argv[1] == 'gontosim':
		method = 'GOntoSim'
	elif sys.argv[1] == 'GOGO' or sys.argv[1] == 'gogo':
		method = 'GOGO'
	elif sys.argv[1] == 'resnik':
		method = 'Resnik' 
	elif sys.argv[1] == 'lin':
		method = 'Lin'
	else:
		method = 'GOntoSim'