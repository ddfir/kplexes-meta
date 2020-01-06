import networkx as nx
import os

for filename in os.listdir('.'):
   print(filename)
   if not (filename.endswith(".adj") or filename.endswith(".csv") or filename.endswith(".txt") or filename.endswith(".txtw") or filename.endswith(".tsv")):
       continue

   print(filename, sum(1 for line in open(filename)))
   aaname = "%s.aa"% os.path.splitext(filename)[0]

   if filename.endswith(".adj"):
       G = nx.read_adjlist(filename)

   if filename.endswith(".csv"):
       G = nx.read_edgelist(filename, delimiter=',', data=(('weight1', int),('weight2',float)))

   if filename.endswith(".txt"):
       G = nx.read_edgelist(filename)

   if filename.endswith(".txtw"):
       G = nx.read_edgelist(filename, data=(('weight1', long),))

   if filename.endswith(".tsv"):
       G = nx.read_edgelist(filename, delimiter='\t', data=(('weight1', float),))

   print("letto", G.number_of_nodes(), G.number_of_edges())
   G.remove_edges_from(G.selfloop_edges())
   print("scrivo", G.number_of_nodes(), G.number_of_edges())
   nx.write_edgelist(G,aaname,data=False)       

   G = nx.read_edgelist(aaname)
   K = list(nx.find_cliques(G))
   print("grafo %s, cricche: %d, cricca max: %d" % (aaname,len(K),max(len(k) for k in K)))
