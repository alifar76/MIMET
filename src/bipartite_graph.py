import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite


eqs = []
kos = []
edge_store = []
compounds = []
graph_labels = {}

infile = open("metagenome_predction.tab", 'rU')
for line in infile:
	if not line.startswith("#"):
		spline = line.strip().split("\t")
		kos.append(spline[0])
for x in kos:
	infile = open("KO_Reactions_Map.txt",'rU')
	for line in infile:
		spline = line.strip().split("\t")
		if x in spline[1:]:
			#print spline[0]+"\t"+x
			eqs.append(spline[0])
reactions = list(set(eqs))
for x in reactions:
	comp_search = re.compile(r"\bC\d\d\d\d\d\b").findall(x)
	comp_search2 = re.compile(r"\bG\d\d\d\d\d\b").findall(x)
	if comp_search:
		for m in comp_search:
			container = []
			container.append(m)
			container.append(x)
			edge_store.append(tuple(container))
			compounds.append(m)
			print x+"\t"+m
	if comp_search2:
		for m in comp_search2:
			container = []
			container.append(m)
			container.append(x)
			edge_store.append(tuple(container))
			compounds.append(m)
			print x+"\t"+m



B = nx.Graph()
B.add_nodes_from(compounds, bipartite=0) # Add the node attribute "bipartite"
B.add_nodes_from(reactions, bipartite=1)
B.add_edges_from(edge_store)

compound_nodes = set(n for n,d in B.nodes(data=True) if d['bipartite']==0)
reaction_nodes = set(B) - compound_nodes

for x in compound_nodes:
	graph_labels[x] = x

for x in reaction_nodes:
	graph_labels[x] = x

pos=nx.graphviz_layout(B)

nx.draw_networkx_nodes(B,pos,node_size=25,nodelist=compound_nodes,node_shape='o')
nx.draw_networkx_nodes(B,pos,node_size=25,nodelist=reaction_nodes,node_color='b',node_shape='s')
nx.draw_networkx_edges(B,pos,alpha=0.5,width=1)	
#nx.draw_networkx_labels(B,pos,graph_labels,font_size=16)
plt.axis('off')
plt.show()