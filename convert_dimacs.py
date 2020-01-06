import networkx as nx
import os

for filename in os.listdir('.'):
   if filename.endswith(".adj"):
       print(filename)
       G = nx.read_adjlist(filename)
       print("letto", G.number_of_nodes(), G.number_of_edges())
       nx.write_edgelist(G,"%s.nde"%filename,data=False)

   if filename.endswith(".csv"):
       print(filename, sum(1 for line in open(filename)))
       G = nx.read_edgelist(filename, delimiter=',', data=(('weight1', int),('weight2',float)))
       print("letto", G.number_of_nodes(), G.number_of_edges())
       G.remove_edges_from(G.selfloop_edges())
       print("scrivo", G.number_of_nodes(), G.number_of_edges())
       nx.write_edgelist(G,"%s.aa"% os.path.splitext(filename)[0],data=False)

   if filename.endswith(".txt"):
       print(filename, sum(1 for line in open(filename)))
       G = nx.read_edgelist(filename)
       print("letto", G.number_of_nodes(), G.number_of_edges())
       G.remove_edges_from(G.selfloop_edges())
       print("scrivo", G.number_of_nodes(), G.number_of_edges())
       nx.write_edgelist(G,"%s.aa"% os.path.splitext(filename)[0],data=False)
