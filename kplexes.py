# coding: utf-8

# Libs, Settings & Utils

import multiprocessing
import signal
import networkx as nx
import math
#import matplotlib.pyplot as plt
from collections import defaultdict
import os
import time
from subprocess import call, DEVNULL
import uuid
from datetime import datetime
from random import shuffle
import sys
import itertools

indir     = "input_data"
workdir   = "working_dir"
extension = ".aa"
berin     = "%s/%s" % (workdir, "block.edgelist")
berout    = "%s/%s" % (workdir, "out")
berexe    = "berlowitz/kplex.py"
timeout   = 3600*12

class Stats(object):
    def __init__(self):
        self.sum_edges  = multiprocessing.Value('i', 0)
        self.blocks     = multiprocessing.Value('i', 0) #blocks done
        self.blockstodo = multiprocessing.Value('i', 0) #blocks todo (same as done if no timeout)
        self.max_nodes  = multiprocessing.Value('i', 0)      
        self.redundant  = multiprocessing.Value('i', 0)  
        self.max_time   = multiprocessing.Value('d', 0)  
        self.size_H     = multiprocessing.Value('i', 0)          
        self.lock       = multiprocessing.Lock()

    def update(self, sum_edges, new_nodes):
        with self.lock:
            self.sum_edges.value += sum_edges
            self.blocks.value += 1
            self.max_nodes.value = max(self.max_nodes.value, new_nodes)   

    def set_blockstodo(self, blockstodo):
        with self.lock:
            self.blockstodo.value = blockstodo     

    def set_size_H(self, size_H):
        with self.lock:
            self.size_H.value = size_H            

    def add_redundant(self):
        with self.lock:
            self.redundant.value += 1

    def update_max_time(self, new_time):
        with self.lock:
            self.max_time.value = max(self.max_time.value, new_time)               

def handler(signum, frame):
    raise Exception("timeout")

def test(G,H):
    color_map = []
    for u in G:
        if u  in H:
            color_map.append('green')
        else:
            color_map.append('gray') 
    nx.draw_spectral(G, node_color = color_map, with_labels = True)
    plt.show()

def neighborhood(G, nodes):
    S = set(nodes)
    for u in nodes:             #TODO: Overl Cliques
        S.update(G.neighbors(u))
    return S

def clique_number(G):
    Cl = defaultdict(int)
    for K in nx.find_cliques(G):
        for u in K:
            Cl[u] = max(Cl[u], len(K))
    return Cl

def all_plexes(G, k, stat): #state-of-the-art exhaustive enum
    stat.update(G.number_of_edges(), G.number_of_nodes())

    #berlowitz
    unique_filename = uuid.uuid4().hex
    filein  = "%s.%s" % (berin, unique_filename)
    fileout = "%s.%s" % (berout, unique_filename)

    nx.write_edgelist(G, filein)
    argstring = "--file=%s --k=%d --type=connected --num_of_kplex=%d --output=%s --size=0"           % (filein, k, 999999999999, fileout)
    call("python2.7 %s %s" % (berexe, argstring), shell=True, stdout=DEVNULL)
    
    with open("%s_connected" % fileout, "r") as f:
        for line in f:
            P = line.rstrip('\n').split(',')
            yield P
            
    os.remove(filein)
    os.remove("%s_connected" % fileout)

    #stub
    #yield from nx.find_cliques(G)

# Algs

def complete(G, clique, candidates):
    S = set(clique)
    for u in sorted(candidates):
        if u in S:
            continue
        compatible = True
        for v in S:
            if not G.has_edge(u,v):
                compatible = False
                break
        if compatible:
            #print(u, "OK")
            S.add(u)
        #else:
        #    print(u, "KO")
    return S

def is_parent_clique(G, target_nodes, candidate_parent):
    parent_aux = complete(G, [min(target_nodes)], target_nodes)
    #print("----")
    parent = complete(G, parent_aux, G.nodes())
    #print("*", candidate_parent, "->", target_nodes, "*", parent_aux, parent)
    return parent == set(candidate_parent)

def graph_filter(G, k, m):
    Co = nx.core_number(G)
    Cl = clique_number(G)
    surv = [u for u in G.nodes() if Co[u] >= m-k and Cl[u] >= math.floor(m/k)]
    H = G.subgraph(surv)
    return H

def graph_filter_recursive(G, k, m):
    Co = nx.core_number(G)
    #Cl = clique_number(G)
    
    Cl = defaultdict(int)
    J  = nx.Graph() #clique graph
    for K in nx.find_cliques(G):
        for u in K:
            Cl[u] = max(Cl[u], len(K))
        if len(K) >= math.floor(m/k):
            J.add_edges_from((u,v) for (u,v) in itertools.combinations(K,2))

    surv = [u for u in G.nodes() if Co[u] >= m-k and Cl[u] >= math.floor(m/k)]
    H = G.subgraph(surv).copy()
    cond = [(u,v) for (u,v) in H.edges() if not J.has_edge(u,v)]
    print("%s remove %d edges" % (datetime.today().strftime('%Y-%m-%d %H:%M:%S'), len(cond)))
    H.remove_edges_from(cond)

    core_only = set()
    cliq_only = set()
    surv = set()
    for K in nx.find_cliques(H):
        B = neighborhood(H, K)
        HB = H.subgraph(B)  

        Co = nx.core_number(HB)
        Cl = clique_number(HB)

        core_only.update([u for u in HB.nodes() if Co[u] >= m-k])
        cliq_only.update([u for u in HB.nodes() if Cl[u] >= math.floor(m/k)])
        surv.update([u for u in HB.nodes() if Co[u] >= m-k and Cl[u] >= math.floor(m/k)])                                

    H = G.subgraph(surv)

    #Co = nx.core_number(G)
    #Cl = clique_number(G)
    #surv = [u for u in G.nodes() if Co[u] >= m-k and Cl[u] >= math.floor(m/k)]
    #H = G.subgraph(surv)
    return H

def large_plexes_sing(H, k, m, stat):
    stat.set_size_H(H.number_of_nodes())
    stat.set_blockstodo(1)

    start = time.time()
    for P in all_plexes(H, k, stat):
        if len(P) >= m:
            yield P
        else:
            stat.add_redundant()
    end = time.time()
    stat.update_max_time(end-start)

def process_clique_batch(H, k, m, cliques, l, r, batch_id, resultdic, stat):
    results = []

    for K in cliques[l:r]:
        B = neighborhood(H, K)
        HB = H.subgraph(B)

        start = time.time()        
        for P in all_plexes(HB, k, stat):
            if len(P) >= m and is_parent_clique(H, P, K):
                results.append(P)
            else:
                stat.add_redundant()   

        end = time.time()
        stat.update_max_time(end-start)                             

    resultdic[batch_id] = results
    #print("sono %d e vado da %d a %d, leggo %d cricche e trovo %d kplessi" % (batch_id, l, r, len(cliques[l:r]), len(resultdic[batch_id])))

def large_plexes_thre(H, k, m, procnum, threaddict, stat):
    stat.set_size_H(H.number_of_nodes())    

    Q = list(nx.find_cliques(H))
    shuffle(Q)

    stat.set_blockstodo(len(Q))

    #---------------------------------------#
    batch_size = math.ceil(len(Q)/procnum)  #
    #batch_size = 1
    #---------------------------------------#
    manager    = multiprocessing.Manager()
    resultdic  = manager.dict()

    if (batch_size > 0):
        for (i,l) in enumerate(range(0, len(Q), batch_size)):
            threaddict[i] = multiprocessing.Process(target=process_clique_batch, args=(H, k, m, Q, l, l + batch_size, i, resultdic, stat), daemon=True)
            threaddict[i].start()
        print("%s started %d threads (batch_size=%d, procnum=%d)" % (datetime.today().strftime('%Y-%m-%d %H:%M:%S'), i+1, batch_size, procnum), file = sys.stderr)        

    for (i,p) in threaddict.items():
        p.join()        
        print("%s ended thread %d" % (datetime.today().strftime('%Y-%m-%d %H:%M:%S'), i), file = sys.stderr)        
        yield from resultdic[i]

def max_plexes_greedy(G, k):
    o = max(len(K) for K in nx.find_cliques(G))
    m = o
    print("[greedy] try", m)
    Plist = list(large_plexes(G,k,m)) #TODO: reuse clique computation if needed
    maxsz = max(len(P) for P in Plist)
    for P in Plist:
        if len(P) == maxsz:
            yield P
                
def max_plexes_binary(G, k):
    o = max(len(K) for K in nx.find_cliques(G))
    m = k * o
    print("[binary] try", m)
    Plist = list(large_plexes(G,k,m)) 
    #--
    while not Plist:
        m = m - math.ceil((m-o) / 2)
        print("[binary] try", m)
        Plist = list(large_plexes(G,k,m)) 
    #--
    maxsz = max(len(P) for P in Plist)
    for P in Plist:
        if len(P) == maxsz:
            yield P

# Run
def experiments():
    signal.signal(signal.SIGALRM, handler)
    t          = datetime.today().strftime('%Y%m%d_%H%M%S')
    log_file   = open("%s.log"%t,"w")
    #log_file   = sys.stdout

    print("%20s %20s %3s %6s %20s %5s %5s %7s %10s %7s %7s %7s %10s %8s %8s %5s %5s %10s %10s %10s %10s %10s" %\
         ("timestamp", "gf", "k", "const", "alg", "o", "m", "n", "edges", "H", "H_core", "H_cliq", "sumedges", "blocks", "maxblock", "found", "smax", "tottime", "redundant", "maxtime", "coretime", "cliquetime"), \
          file = log_file)
    log_file.flush()

    all_files = (os.path.join(basedir, filename) for basedir, dirs, files in os.walk(indir) for filename in files)
    sorted_files = sorted(all_files, key = os.path.getsize)
    print(sorted_files)

    for filename in sorted_files:
        if filename.endswith(extension):
            G  = nx.read_edgelist(filename)
            #G.remove_edges_from(G.selfloop_edges())
            gf = os.path.basename(filename).split('.')[0]                
        
            for k in [2, 3]:
                signal.alarm(timeout)
                try:   
                    o  = max(len(K) for K in nx.find_cliques(G))
                except: 
                    t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                    print("%20s %20s %3s %6s %20s %5s %5s %7d %10d %7s %7s %7s %10s %8s %8s %33s %10s" %\
                         (t, gf, k, '----', "INIT", "", "", G.number_of_nodes(), G.number_of_edges(), '', '', '', '', '', '', "*****"+str(timeout)+"******", ""), \
                          file = log_file)
                    log_file.flush()
                signal.alarm(0) 

                ## filtering stats ##
                for (fatt,label) in []:#[(100,'100'), (50,'50'), (10,'10')]:
                    m  = max(k**2, math.ceil(fatt))

                    signal.alarm(timeout)
                    try:   
                        core_time = time.time()
                        Co = nx.core_number(G)
                        core_time = time.time() - core_time

                        clique_time = time.time()                        
                        Cl = clique_number(G)
                        clique_time = time.time() - clique_time    

                        core_only = len([u for u in G.nodes() if Co[u] >= m-k])
                        cliq_only = len([u for u in G.nodes() if Cl[u] >= math.floor(m/k)])

                        surv = [u for u in G.nodes() if Co[u] >= m-k and Cl[u] >= math.floor(m/k)]
                        H = G.subgraph(surv)
                        size_H_local = H.number_of_nodes()                         

                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3s %6s %20s %5d %5d %7d %10d %7d %7d %7d %10s %8s %8s %33s %10s %10.4f %10.4f" %\
                             (t, gf, k, label, "F-STATS", o, m, G.number_of_nodes(), G.number_of_edges(), size_H_local, core_only, cliq_only, '', '', '', "", "", core_time, clique_time), \
                              file = log_file)
                        log_file.flush()

                        sum_edges_local = 0

                        dist = defaultdict(int)
                        blocks_local, max_nodes_local, max_block, max_block_clique = 0, 0, 0, 0
                        for K in nx.find_cliques(H):
                            B = neighborhood(H, K)
                            HB = H.subgraph(B)
                            sum_edges_local += HB.number_of_edges()
                            blocks_local += 1
                            if len(B) > max_nodes_local:
                                max_nodes_local = len(B)
                                max_block = B
                                max_block_clique = K

                        recursive_time = time.time()
                        H = graph_filter_recursive(G, k, m)
                        recursive_time = time.time() - recursive_time

                        size_H_local = H.number_of_nodes()                                        

                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7d %7d %10d %8d %8d %11s %10.2f %10s %10s %10s" %\
                              (t, gf, k, label, "2SF-STATS", o, m, G.number_of_nodes(), G.number_of_edges(), size_H_local, core_only, cliq_only, sum_edges_local, blocks_local, max_nodes_local, "", recursive_time, "", "", ""), \
                              file = log_file)                        
                        log_file.flush()                                                           
                    except Exception as e:
                        print(e) 
                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7s %7s %10s %8d %42s" %\
                              (t, gf, k, label, "2SF-STATS", o, m, G.number_of_nodes(), G.number_of_edges(), size_H_local, "", "", "", blocks_local, "*****"+str(timeout)+"******"), \
                              file = log_file)                        
                        log_file.flush()
                        continue
                    signal.alarm(0) 

                    ## maxblock stats ##
                    #signal.alarm(timeout)
                    #try:
                    #    start = time.time()
                    #    if (size_H_local > 0):
                    #        nx.write_edgelist(H.subgraph(max_block), berin)

                    #        argstring = "--file=%s --k=%d --type=connected --num_of_kplex=%d --output=%s --size=0"           % (berin, k, 999999999999, berout)
                    #        call("python2.7 %s %s" % (berexe, argstring), shell=True, stdout=DEVNULL)

                    #        with open("%s_connected" % berout, "r") as f:
                    #            for line in f:
                    #                P = line.rstrip('\n').split(',')
                    #                check = len(P) >= m and is_parent_clique(H, P, K)
                    #    end = time.time()
                    #    max_time_local = end-start
                    #except: 
                    #    continue
                    #signal.alarm(0)
                    ###################

                ## enumeration stats ##
                for (fatt,label) in [(1.0*o,'100%'), (0.75*o,'75%'),(0.5*o,'50%')]: #(o+1,'max'), 
                    m  = max(k**2, math.ceil(fatt))                    
                       
                    ##Filter only
                    size_H, redundant, max_time = 0, 0, 0
                    stat = Stats()
                    signal.alarm(timeout)
                    try:
                        start = time.time()
                        Pnum, Smax = 0, 0 
                        for P in large_plexes_sing(graph_filter_recursive(G, k, m), k, m, stat):
                            Pnum +=1
                            print(P, file=sys.stderr)
                            if (len(P) > Smax):
                                Smax = len(P)
                        end = time.time()
                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7d %7d %10d %8d %8d %5d %5d %10.2f %10d %10.2f" %\
                             (t, gf, k, label, "VANIL(H)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, 0, 0, stat.sum_edges.value, stat.blocks.value, stat.max_nodes.value, Pnum, Smax, end-start, stat.redundant.value, stat.max_time.value), \
                              file = log_file)
                        log_file.flush()                    
                    except Exception as e: 
                        print(e, file = sys.stderr)   
                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7s %7s %10s %8d %42s" % (t, gf, k, label, "VANIL(H)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, "", "", "", stat.blocks.value, "*****"+str(timeout)+"******"), \
                              file = log_file)
                        log_file.flush()                    
                    signal.alarm(0)
                    ###################

                    ##Filter + Blocks
                    size_H = 0
                    stat = Stats()                
                    threaddict = {}
                    signal.alarm(timeout)
                    try:
                        start = time.time()
                        Pnum, Smax = 0, 0 
                        for P in large_plexes_thre(graph_filter_recursive(G, k, m), k, m, 31, threaddict, stat):
                            Pnum +=1
                            if (len(P) > Smax):
                                Smax = len(P)
                        end = time.time()
                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7d %7d %10d %8d %8d %5d %5d %10.2f %10d %10.2f" %\
                             (t, gf, k, label, "BLOCKS(H)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, 0, 0, stat.sum_edges.value, stat.blocks.value, stat.max_nodes.value, Pnum, Smax, end-start, stat.redundant.value, stat.max_time.value), \
                              file = log_file)
                        log_file.flush()                    
                    except:
                        for (i,p) in threaddict.items():
                            if (p.is_alive()):
                                p.terminate()
                                p.join()        
                        t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                        print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7s %7s %1s %8d/%8d %42s" % (t, gf, k, label, "BLOCKS(H)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, "", "", "", stat.blocks.value, stat.blockstodo.value, "*****"+str(timeout)+"******"), \
                              file = log_file)
                        log_file.flush()
                        #########
                        #break    
                        #########                
                    signal.alarm(0)  
                    ###################           

                ## full enumeration stats ##
                ##Blocks Only
                '''
                fatt  = k**2
                label = 'all'
                m     = fatt
                size_H = 0
                stat = Stats()                
                threaddict = {}                
                signal.alarm(timeout)
                try:
                    start = time.time()
                    Pnum, Smax = 0, 0 
                    for P in large_plexes_thre(G, k, m, 31, threaddict, stat):
                        Pnum +=1
                        if (len(P) > Smax):
                            Smax = len(P)
                    end = time.time()
                    t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                    print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7d %7d %10s %8s %8s %5d %5d %10.2f %10d %10.2f" %\
                         (t, gf, k, label, "BLOCKS(G)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, 0, 0, stat.sum_edges.value, stat.blocks.value, stat.max_nodes.value, Pnum, Smax, end-start, stat.redundant.value, stat.max_time.value), \
                          file = log_file)
                    log_file.flush()                    
                except:
                    for (i,p) in threaddict.items():
                        if (p.is_alive()):
                            p.terminate()
                            p.join()                        
                    t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                    print("%20s %20s %3d %6.2f %20s %117s" % (t, gf, k, fatt, "BLOCKS(G)", "*****"+str(timeout)+"******"), \
                          file = log_file)
                    log_file.flush()                    
                signal.alarm(0)
                '''
                ###################  

                ##Plain Berlowitz
                '''
                size_H, redundant, max_time = 0, 0, 0
                stat = Stats()
                signal.alarm(timeout)
                try:
                    start = time.time()
                    Pnum, Smax = 0, 0 
                    for P in large_plexes_sing(G, k, m, stat):
                        Pnum +=1
                        if (len(P) > Smax):
                            Smax = len(P)
                    end = time.time()
                    t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                    print("%20s %20s %3d %6s %20s %5d %5d %7d %10d %7d %7d %7d %10d %8d %8d %5d %5d %10.2f %10d %10.2f" %\
                         (t, gf, k, label, "VANIL(G)", o, m, G.number_of_nodes(), G.number_of_edges(), stat.size_H.value, 0, 0, stat.sum_edges.value, stat.blocks.value, stat.max_nodes.value, Pnum, Smax, end-start, stat.redundant.value, stat.max_time.value), \
                          file = log_file)
                    log_file.flush()                    
                except: 
                    t = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                    print("%20s %20s %3d %6.2f %20s %117s" % (t, gf, k, fatt, "VANIL(G)", "*****"+str(timeout)+"******"), \
                          file = log_file)
                    log_file.flush()                    
                signal.alarm(0) 
                '''
                ###################   

            print("",file = log_file)
            log_file.flush()

        print("",file = log_file)
        log_file.flush()

    log_file.close()

experiments()
