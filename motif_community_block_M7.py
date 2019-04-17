import networkx as nx
import random as rd
import matplotlib.pyplot as plt
from collections import deque
from scipy.cluster.hierarchy import dendrogram

##----N is the number of nodes in G; SIZE is the number of nodes in one group;
##----Pin is the probability to generate a motif instance in the same community; Pout is the probability to generate a motif instance between different communities
DG = nx.DiGraph()
G = nx.Graph()
N=60
SIZE=30
Pin=0.2
Pout=0.002
Peout=0.0
MP=3.0
INF=1000
Stack=[]


###---sugraph that can become a motif instance
DG0=nx.DiGraph()
DG1=nx.DiGraph()
DG2=nx.DiGraph()
DG3=nx.DiGraph()
DG4=nx.DiGraph()
for u in range(3):
    DG0.add_node(u)
    DG1.add_node(u)
    DG2.add_node(u)
    DG3.add_node(u)
    DG4.add_node(u)
DG1.add_edge(0,2)
DG2.add_edge(0,2)
DG2.add_edge(1,2)
DG3.add_edge(0,1)
DG3.add_edge(1,0)
DG4.add_edge(0,1)
DG4.add_edge(1,0)
DG4.add_edge(0,2)



##---generate a motif instance
def add_motif(i,j,k):
    nodes=[i,j,k]
    DH=DG.subgraph(nodes)
    elist=list(DH.edges())
    if nx.is_isomorphic(DH,DG0):
        s=rd.choice(nodes)
        nodes.remove(s)
        DG.add_edge(nodes[0],nodes[1])
        DG.add_edge(nodes[1],nodes[0])
        DG.add_edge(nodes[0],s)
        DG.add_edge(nodes[1],s)
    elif nx.is_isomorphic(DH,DG1):
        e=elist[0]
        nodes.remove(e[0])
        nodes.remove(e[1])
        t=nodes[0]
        DG.add_edge(e[0],t)
        DG.add_edge(t,e[0])
        DG.add_edge(t,e[1])
    elif nx.is_isomorphic(DH,DG2):
        e1=elist[0]
        e2=elist[1]
        DG.add_edge(e1[0],e2[0])
        DG.add_edge(e2[0],e1[0])
    elif nx.is_isomorphic(DH,DG3):
        e=elist[0]
        nodes.remove(e[0])
        nodes.remove(e[1])
        t=nodes[0]
        DG.add_edge(e[0],t)
        DG.add_edge(e[1],t)
    elif nx.is_isomorphic(DH,DG4):
        e1=elist[0]
        e2=elist[1]
        e3=elist[2]
        if e1[1]==e2[0] and e2[1]==e1[0]:
            nodes.remove(e3[0])
            nodes.remove(e3[1])
            t=nodes[0]
            DG.add_edge(t,e3[1])
        elif e1[1]==e3[0] and e1[0]==e3[1]:
            nodes.remove(e2[0])
            nodes.remove(e2[1])
            t=nodes[0]
            DG.add_edge(t,e2[1])
        elif e2[1]==e3[0] and e2[0]==e3[1]:
            nodes.remove(e1[0])
            nodes.remove(e1[1])
            t=nodes[0]
            DG.add_edge(t,e1[1])
        
        

##----draw a graph from the stochastic block model
def initial_network():
    for i in range(SIZE):
        DG.add_node(i,label=i+1,g=1)
        G.add_node(i,val=0.0)
    for i in range(SIZE,N):
        DG.add_node(i,label=i+1,g=2)
        G.add_node(i,val=0.0)
    rd.seed()
    for i in range(N):
        for j in range(i+1,N):
            if DG.node[i]['g']!=DG.node[j]['g']:
                re=rd.random()
                if re<=Peout:
                    DG.add_edge(i,j)
            for k in range(j+1,N):
                tmp=rd.random()
                if DG.node[i]['g']==DG.node[j]['g'] and DG.node[j]['g']==DG.node[k]['g']:
                    if tmp<=Pin:
                        add_motif(i,j,k)
                else:
                    if tmp<=Pout:
                        add_motif(i,j,k)
    for e in list(DG.edges()):
        G.add_edge(e[0],e[1],col='r',w=0.0)



#test motif whether a subgraph is a motif instance of M_7
def is_motif(i,j,k):
    if i in list(DG.successors(j)) and j in list(DG.successors(i)) and k in list(DG.successors(i)) and i not in list(DG.successors(k)) and k in list(DG.successors(j)) and j not in list(DG.successors(k)):
        return True
    elif i in list(DG.successors(k)) and k in list(DG.successors(i)) and j in list(DG.successors(i)) and i not in list(DG.successors(j)) and j in list(DG.successors(k)) and k not in list(DG.successors(j)):
        return True
    elif k in list(DG.successors(j)) and j in list(DG.successors(k)) and i in list(DG.successors(j)) and j not in list(DG.successors(i)) and i in list(DG.successors(k)) and k not in list(DG.successors(i)):
        return True
    else:
        return False


#assign weights to each node
def node_motif_weight():
    for i in range(N):
        for j in range(i+1,N):
            for k in range(j+1,N):
                if is_motif(i,j,k):
                    G.node[i]['val']+=1/MP
                    G.node[j]['val']+=1/MP
                    G.node[k]['val']+=1/MP
                    G.edges[i,j]['col']='b'
                    G.edges[i,k]['col']='b'
                    G.edges[j,k]['col']='b'
                    
                    

#remove edges that are not in M_7 instances
def non_motif_edge_remove():
    temps=[]
    for e in G.edges():
        if G.edges[e[0],e[1]]['col']=='r':
            temps.append(e)
    for e in temps:
        G.remove_edge(e[0],e[1])
        
    

#using BFS to compute motif betweeness of edges
def BFS_betweeness(u):
    colors=[]
    pres=[]
    ws=[]
    dis=[]
    pnums=[]
    for i in range(N):
        colors.append(0)
        dis.append(INF)
        pres.append([])
        ws.append(0.0)
        pnums.append(0)
    colors[u]=1
    dis[u]=0
    pnums[u]=1
    queue = deque()
    stack=[]
    queue.append(u)
    stack.append(u)
    while len(queue)>0:
        i=queue.popleft()
        for j in list(G[i]):
            if colors[j]==0:
                colors[j]=1
                dis[j]=dis[i]+1
                pres[j].append(i)
                queue.append(j)
                stack.append(j)
            elif colors[j]==1 and dis[j]==dis[i]+1:
                pres[j].append(i)
        if len(pres[i])>0:
            for j in pres[i]:
                pnums[i]+=pnums[j]
        colors[i]=2
    while len(stack)>0:
        v=stack.pop()
        weight=(G.node[u]['val']*G.node[v]['val']+ws[v])*1.0/(pnums[v]*1.0)
        for i in pres[v]:
                tmp=weight*pnums[i]*1.0
                G.edges[i,v]['w']+=tmp
                ws[i]+=tmp



    

#computing motif betweeness of edges
def compute_betweeness():
    for u in G.nodes():
        if G.degree[u]>0 and G.node[u]['val']>0:
            BFS_betweeness(u)



#delete the edge with largest motif betweenness
def del_max_edge():
    max_b=-1
    edge=[-1,-1]
    for e in list(G.edges()):
        if G.edges[e[0],e[1]]['w']>max_b:
            max_b=G.edges[e[0],e[1]]['w']
            edge[0]=e[0]
            edge[1]=e[1]
        G.edges[e[0],e[1]]['w']=0.0
    G.remove_edge(edge[0],edge[1])
    Stack.append(edge)




            
#showing clusters
def clusters_present():
    cnum=nx.number_connected_components(G)
    if cnum>=6 and cnum<=10:
        for comp in sorted(nx.connected_components(G), key=len, reverse=True):
            print(*comp)
            print("%%%%%%%%%%%%%%%%%%%%%%%%")
        print("------------------------------------------------------------------")
   


#decpit dendrogram
def draw_dengrogram():
    global N,Stack
    classes=list(range(N))
    cidx=N
    clusters=dict()
    dists=dict()
    Z=[]
    ls=list(range(N))
    for i in range(N):
        ls[i]=i+1
        classes[i]=i
        tmp=[i,]
        clusters[i]=tmp
        dists[i]=0.0
    while len(Stack)>0:
        e=Stack.pop()
        c1=classes[e[0]]
        c2=classes[e[1]]
        tmp1=[]
        if(c1!=c2):
            for u in clusters[c1]:
                tmp1.append(u)
                classes[u]=cidx
            for u in clusters[c2]:
                tmp1.append(u)
                classes[u]=cidx
            clusters[cidx]=tmp1
            dists[cidx]=max(dists[c1],dists[c2])+1.0
            ztmp=[]
            ztmp.append(c1)
            ztmp.append(c2)
            ztmp.append(dists[cidx])
            ztmp.append(len(clusters[cidx]))
            Z.append(ztmp)
            cidx+=1
    fig=plt.gca()
    fig.axes.xaxis.set_ticklabels([])
    fig.axes.yaxis.set_ticklabels([])
    dn=dendrogram(Z,labels=ls,leaf_font_size=10)
    plt.show()




##---the algorithm to detect motif-based community structure
def detect_motif_clusters():
    node_motif_weight()
    non_motif_edge_remove()
    while G.number_of_edges()!=0:
        compute_betweeness()
        del_max_edge()
        #clusters_present()
    
        


def main():
    initial_network()
    #nx.write_gml(DG, "blockmodel2.gml")
    #nx.draw_networkx(DG,with_labels=True,width=0.38,node_color=colors,labels=labels,node_size=100,font_color='r',font_weight='bold',pos=positions)
    #nx.draw_spring(DG,with_labels=True,width=0.38,node_color=colors,labels=labels,node_size=100,font_color='b')
    #limits = plt.axis('off')
    #plt.show()
    detect_motif_clusters()
    draw_dengrogram()
    
     



main()
                        
                        
    
