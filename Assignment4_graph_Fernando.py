#!/usr/bin/env python

import random

"""
Author: Fernando Bueno Gutierrez
Student nr:890605143090
Script to: Determine whether a graph is Eulerian, and if so, find eulerian 
           cycles and eulerian paths. The script also creates graphs from a
           spectrum
"""


def list_of_tupples(graph):
    """
    Transforms graph diccionary into list of tupples.

    Arguments
    ---------
    graph: dict, graph, where the key is the node and the values are the nodes 
           that it points to.

    Returns
    ---------
    tupples: list of tupples, equivalent to "graph".
    """

    tupples = []
    for tuple in graph.items():
        for i in range(0, len(tuple[1])):
            tupples.append((tuple[0],tuple[1][i]))


    return tupples



def is_eulerian(graph):
    """
    Determines whether a graph is eulerian

    Arguments
    ---------
    graph: dict, graph, where the key is the node and the values are the nodes 
           that it points to.

    Returns
    ---------
    eulerian: logical, is "True" if the graph is eulerian.

    balances: list of numbers, where each number refers to the "balance*" 
              of a node.

    balance*: The difference between the number of times that the node is 
              pointed (by other nodes in the graph) and the number of nodes 
              that the node points to.
    """

    eulerian= False
    list_called = []	# #times a node is called
    for key, value in graph.iteritems():
	    list_called.extend(value)

    balances=[]
    for node in graph:
	    calls = len(graph[node])	# calls: #times a node calls
	    balance = list_called.count(node) - calls	# balance: degree
	    balances.append(balance)

    if (all(balance == 0 for balance in balances)):
	    eulerian=True

    return eulerian, balances



def has_eulerian_path(graph):
    """
    Determines whether a graph has any eulerian path.

    Arguments
    ---------
    graph: dict, graph, where the key is the node and the values are the nodes 
           that it points to.

    Returns
    ---------
    returns: logical,is "True" if the graph has any eulerian path.
    """
    nodes = graph.keys()
    non_balanced = 0
    semi_balanced = 0

    balances=is_eulerian(graph)[1]

    for idx, node in enumerate(graph):
        if abs(balances[idx]) == 1:
	        semi_balanced =+ 1
        if abs(balances[idx]) > 1:
	        non_balanced =+ 1
    if non_balanced==0 & semi_balanced<3:
        return True

    return False



def graph_spec(spect):
    """
    Generates graph in the form of a dictionary from a list of k-mers

    Arguments
    ---------
    spect: list of a k-mers, also called "spectrum"

    Returns
    ---------
    graph_spect: dict, graph generated from the spectrum.
    """


    keys=[]
    values=[]
    graph_spect=dict()
    for idx, lmer in enumerate(spect):
        key=lmer[0:2]
        value=lmer[1:3]
        if idx==0:
            graph_spect[key] = [value]			
        else:
            if key in graph_spect:
     		    graph_spect[key].append(value)
            else:	
                graph_spect[key] = [value]
    for idx, lmer in enumerate(spect):
        key=lmer[1:3]
        if key not in graph_spect.keys():
            value='NaN'
            if idx==0:
                 graph_spect[key] = [value]			
            else:
                if key in graph_spect:
         			graph_spect[key].append(value)
                else:	
                 	graph_spect[key] = [value]				


    return(graph_spect)	




def possible_nodes(current,graph):
    """
    Auxiliary function that returns the nodes to which another node 
    is pointing (so, the "possible next nodes").

    Arguments
    ---------
    current: tupple, node of interest

    graph:   dict, graph, where the key is the node and the values are the 
             nodes that it points to.

    Returns
    ---------
    possible_nodes: List of tupples, the nodes in "graph" the node "current" 
                    points to. 
    """

    possible_nodes=[]
    for node in graph:
	    if node[0] == current[1]:
		    possible_nodes.append(node)
    return possible_nodes






def find_path(graph):
    """
    Auxiliary function that returns one of the many possible paths in a 
    graph. Note that the path may not be eulerian.

    Arguments
    ---------
    graph: dict, graph, where the key is the node and the values are the nodes 
           that it points to.

    Returns
    ---------
    path: List of char, each char is a nodes in "graph". The order in which the 
          leters appear in the list corresponds to the order of the path
          (i.e. "['A','B']" means that direction is from node 'A' to node 'B').
    """

    result = None
    while result is None:
	    try:
		    result = None
		    graph_tuple=(list_of_tupples(graph)) 
            #graph_tuple: graph as list of tupples
		    current=random.choice(graph_tuple) #current: first edge
		    graph_tuple.remove(current)	
            #Remove, so that the edge is not used again
		    path=[current[0]]
		    current_list=[current]

		    for i in range(len(graph_tuple)):
			    neighS=possible_nodes(current,graph_tuple)
                #neighS:neighbours, possible nodes to follow current
			    next=random.choice(neighS) 
                #One of the possible nodes, will be the next
			    neighS.remove(next)
			    graph_tuple.remove(next)
                #Remove, so that the node is not used again
			    path.append(next[0])
			    current_list=[current_list,next]
			    current=next
			    result=path
	    except:
		    pass
    return path


def return_eulerian(graph):
    """
    Returns one of the possible eulerian paths in a graph. 

    Arguments
    ---------
    graph: dict, graph, where the key is the node and the values are the nodes 
           that it points to.

    Returns
    ---------
    eulerian: List of char, each char is a nodes in "graph". The order in which 
              the chars appear in the list corresponds to the order of the 
              eulerian path. (i.e. "['A','B']" means that direction is from 
              node 'A' to node 'B').
    """	
    graph_tuple=(list_of_tupples(graph))
    stop=len(graph_tuple)-1 #stop: integer, number of different nodes
    finish=False #whether all the different nodes have been added 
    while not finish:
        pa=find_path(graph)
        finish=False
        if len(pa)>stop:
            finish=True
            eulerian=pa
    if finish:
        #close the cycle by adding the first node, so that start and end match
        if isinstance(eulerian[0], basestring):
            first = str(eulerian[0])
            first = ["".join(first)]
        else:   
            first = [eulerian[0]]
        eulerian.extend(first)
        return eulerian


if __name__ == "__main__":

    # GRAPH FROM FIG 8.22
    graph_822 = {'A':['B'],'B':['C'],'I':['H'],'H':['F'],'F':['G','E'],\
        'C':['I','J'],'G':['A'],'E':['J'],'J':['F','D'],'D':['C']}
    # A SLIGHTLY BIGGER GRAPH, NEEDED FOR Q8
    bigger_graph = {5:[6],6:[7],10:[11],11:[4],4:[5,3],\
        7:[10,9],3:[9,1],9:[4,8],8:[7],1:[2], 2:[3]}
    # SPECTRUM FROM FIG 8.20
    s = ['ATG','TGG','TGC','GTG','GGC','GCA','GCG','CGT']

    # Put function calls, print statements etc. to answer the questions here
    # When we run your script we should see the answers on screen (or file)


# Question 1
print "\nQUESTION 1 \n"
if is_eulerian(graph_822)[0]:
    print "The graph is Eulerian"
else:
    print "The graph is not Eulerian"


# Question 2
print "\nQUESTION 2 \n"
if has_eulerian_path(graph_822):
    print "The graph contains an Eulerian path"
else:
    print "The graph does not contain an Eulerian path"


# Question 3
print "\nQUESTION 3 \n"
print(return_eulerian(graph_822))



# Question 4
print "\nQUESTION 4 \n"
print "run1 ",return_eulerian(graph_822)
print "run2 ",return_eulerian(graph_822)
print "run3 ",return_eulerian(graph_822)




# Question 5
print "\nQUESTION 5 \n"
grp=graph_spec(s)
for k, v in grp.items():
 print k, v


#Question 6
print "\nQUESTION 6 \n"
if is_eulerian(grp)[0]:
    print "The graph is Eulerian"
else:
    print "The graph is not Eulerian" 

if has_eulerian_path(grp):
    print "The graph contains an Eulerian path"
else:
    print "The graph does not contain an Eulerian path"



#Question 7
print "\nQUESTION 7 \n"
print(return_eulerian(grp))
L=return_eulerian(grp)


L2=[]
for ele in L:
    L2.append(str(ele)[0])
L2.append(str(L[-1])[1])

res = "".join(L2)

print(res)




#Question 8
print "\nQUESTION 8 \n"
print(return_eulerian(bigger_graph))




 
