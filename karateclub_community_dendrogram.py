import networkx as nx
import random as rd
from collections import deque
from scipy.cluster.hierarchy import dendrogram
from matplotlib import pyplot as plt

G=nx.Graph()
N=34
MP=2.0
INF=1000
#########dendrogram#############
Stack=[]

def initial_network(f):
    for i in range(N):
        G.add_node(i,val=0.0)
    for line in f:
        strli=line.split()
        a=int(strli[0])
        b=int(strli[1])
        G.add_edge(a-1,b-1,w=0.0)




#assign weights to each node
def node_motif_weight():
    for u in list(G.nodes()):
        G.node[u]['val']=G.degree[u]/MP
       
                        
                    
                    

        
#BFS to compute edge betweeness
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
        weight=(G.node[u]['val']*G.node[v]['val']+ws[v])/pnums[v]
        for i in pres[v]:
                tmp=weight*pnums[i]
                G.edges[i,v]['w']+=tmp
                ws[i]+=tmp
    
            
    

#compute betweeness of edges
def compute_betweeness():
    for u in G.nodes():
        if G.degree[u]>0:
            BFS_betweeness(u)



#delete largest 
def del_max_edge():
    global Stack
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




            
#show partitioned clusters
def clusters_present():
    cnum=nx.number_connected_components(G)
    if cnum>=2:
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
    dn=dendrogram(Z,labels=ls)
    plt.show()
        
                
                
    



def detect_motif_clusters():
    node_motif_weight()
    while G.number_of_edges()!=0:
        compute_betweeness()
        del_max_edge()
        #clusters_present()
        


def main():
    f=open('karate-club.txt', 'r')
    initial_network(f)
    detect_motif_clusters()
    draw_dengrogram()
    f.close()
    
    


main()
