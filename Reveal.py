#!/usr/bin/env python


"""
Author: Fernando Bueno Gutierrez
Script to: Implement the REVEAL (REVerse Engineering ALgorithm) algorithm for 
           gene network inference
"""

from __future__ import division
from numpy import log2
from itertools import combinations
from copy import deepcopy
import math
import itertools
import numpy as np
import re
import time
import csv
import random
import pandas as pd 
import os

start_time = time.time()




eps = 1e-100

def read_tsv_file(tsv):
    """Return dicts of inputs (t = 0,...,T-1) and outputs (t = 1,...,T)

    Arguments
    ---------
    tsv: file, data stream in tsv format with each row a time point and each
        column a gene

    Returns
    ---------
    inputs:  dict of {gene: [float, ->]}, with each float the binary
             state of the gene at a certain time point. The list is ordered
             by time point, and the last time point is not included

    outputs: dict of {gene: [float, ->]}, with each float the binary
             state of the gene at a certain time point. The list is ordered
             by time point, and the first time point is not included
    """

    header = tsv.readline().rstrip() #parse header
    nodes  = header.split("\t")
    
    inputs  = dict([(n,[]) for n in nodes[1:]])
    for line in tsv.readlines():
        vals = line.rstrip().split("\t")
        for i in range(1, len(nodes)):
            inputs[nodes[i]].append(float(vals[i]))
    outputs = deepcopy(inputs)
    
    for i in range(1,len(nodes)):
      del inputs[nodes[i]][-1] #delete last time point
      del outputs[nodes[i]][0] #delete first time point
    
    return inputs, outputs

def check_table(ins, outs):
    """Print entropy and mutual information from Fig 5 of Liang et al, 1998 [1]

    Arguments
    ---------
    ins:    dict of {gene: [float, ->]}, with each float the binary
            state of the gene at a certain time point. The list is ordered
            by time point, and the last time point is not included

    outs:   dict of {gene: [float, ->]}, with each float the binary
            state of the gene at a certain time point. The list is ordered
            by time point, and the first time point is not included

    Returns
    ---------
    Print statements, displaying entropy and mutual information from Fig 5
        of Liang et al, 1998 [1]
    """

    print "H(A), H(B), H(C) =",       \
          entropy(ins['A']), entropy(ins['B']), entropy(ins['C'])
    print "H(A,B), H(B,C), H(A,C) =", \
          entropy(ins['A'],ins['B']), entropy(ins['B'],ins['C']), \
          entropy(ins['A'],ins['C'])
    print "H(A,B,C) =",               \
          entropy(ins['A'],ins['B'],ins['C'])

    print
    print "H(A') =", entropy(outs['A'])
    print "H(A',A), H(A',B), H(A',C) =", \
        entropy(outs['A'],ins['A']), \
        entropy(outs['A'],ins['B']), \
        entropy(outs['A'],ins['C'])
    print "M(A',A), M(A',B), M(A',C) =", \
        mutual_information(outs['A'],ins['A']), \
        mutual_information(outs['A'],ins['B']), \
        mutual_information(outs['A'],ins['C'])

    print
    print "H(B') =", entropy(outs['B'])
    print "H(B',A), H(B',B), H(B',C) =", \
        entropy(outs['B'],ins['A']), \
        entropy(outs['B'],ins['B']), \
        entropy(outs['B'],ins['C'])
    print "H(B',[A,B]), H(B',[B,C]), H(B',[A,C]) =", \
        entropy(outs['B'],ins['A'],ins['B']), \
        entropy(outs['B'],ins['B'],ins['C']), \
        entropy(outs['B'],ins['A'],ins['C'])
    print "M(B',A), M(B',B), M(B',C) =", \
        mutual_information(outs['B'],ins['A']), \
        mutual_information(outs['B'],ins['B']), \
        mutual_information(outs['B'],ins['C'])
                              
    print
    print "H(C') =", entropy(outs['C'])
    print "H(C',A), H(C',B), H(C',C) =", \
        entropy(outs['C'],ins['A']), \
        entropy(outs['C'],ins['B']), \
        entropy(outs['C'],ins['C'])
    print "H(C',[A,B]), H(C',[B,C]), H(C',[A,C]) =", \
        entropy(outs['C'],ins['A'],ins['B']), \
        entropy(outs['C'],ins['B'],ins['C']), \
        entropy(outs['C'],ins['A'],ins['C'])
    print "H(C',[A,B,C]) =", entropy(outs['C'],ins['A'],ins['B'],ins['C'])
    print "M(C',A), M(C',B), M(C',C) =", \
        mutual_information(outs['C'],ins['A']), \
        mutual_information(outs['C'],ins['B']), \
        mutual_information(outs['C'],ins['C'])
    print "M(C',[A,B]), M(C',[B,C]), M(C',[A,C]) =", \
        mutual_information(outs['C'],ins['A'],ins['B']), \
        mutual_information(outs['C'],ins['B'],ins['C']), \
        mutual_information(outs['C'],ins['A'],ins['C'])
    print "M(C',[A,B,C]) =",\
        mutual_information(outs['C'],ins['A'],ins['B'],ins['C'])
    


def all_vectors(*arg):
    """
    Auxiliary function for the function entropy. Returns the proportion of
    time points that show a particular combination of expression of genes (i.e,
    "011" means that 2 out of 3 genes are active at a given point).

    Arguments
    ---------
    arg*: un-fixed number of lists. Each list one integer (0/1) that represents 
          the expression values of a gene in certain time points.
  
    Returns
    ---------
    table_combinations: dict, keys are strings that represent combinations of 
                        expression values (0/1) given *genes ,values are the 
                        proprtion of time points that show that combintaion of
                        expression values (i.e. if two genes are give, all
                        possible combinations are 00,01,10,11)
    """

    table_combinations = {}
    for i in range(0,len(arg[0])):
        key_each_time=[]
        for j in range(0,len(arg)):
            key_each_time.append(arg[j][i])
        key_each_time=''.join(map(str, key_each_time))

        if key_each_time not in table_combinations.keys():
                    table_combinations[key_each_time] = 1
        else:
            table_combinations[key_each_time] += 1

    

    for key, value in table_combinations.items():
        table_combinations[key] = value / len(arg[0])

    return table_combinations

 
   
def entropy(*arg):
    """
    Return the value of entropy corresponding to the a combination of genes 

    Arguments
    ---------
    *arg: list of lists of non-fixed length. Each list is the expression values 
          of a gene. Each integer is an expression value (0 or 1) in a given  
          time point. The outer list represents a combination of genes.

    Returns
    ---------
    MI: Integer, mutual information of the combination of genes

    Note
    ---------
    According to [1] a maximum of 0 or 1 genes should be in output-step, 
    the rest should be in input-state.

    """
    
    table_combinations=all_vectors(*arg)


    to_be_summed=[]
    for value in table_combinations.values():
        combination=-value*math.log(value,2) #base 2 as in paper
        to_be_summed.append(combination)
    
    entropy=round(sum(to_be_summed),4)

    return entropy



def mutual_information(*arg):
    """
    Return the value of mutual information corresponding to the input-state of  
    a combination of genes and the output-state of one gene of interest.

    Arguments
    ---------
    *arg: list of lists of non-fixed length. Each list is the expression values 
          of a gene. Each integer is an expression value (0 or 1) in a given  
          time point. The outer list represents a combination of genes. The 
          first inner list corresponds to the output-state of the gene of 
          interest.

    Returns
    ---------
    MI: Integer, mutual information of the combination of genes

    """

    In=[] #expression values of genes in input-state
    for i in range(1,len(arg)):
        In.append(arg[i])


    In_Out=[] #expression values of genes in input-state and output-state
    for i in range(0,len(arg)):
        In_Out.append(arg[i])

    MI = round(entropy(arg[0])+entropy(*In)-entropy(*In_Out),2)
    # MI: mutual information
    # The '*' is essential, otherwise it will call entropy with a liso of lists
    # instead of a sequence of lists, and result will be different

    return MI




def combinations(inp_remaining,k):
    """
    Auxiliary function that returns a list with all possible k-size unique
    combinations of a list of genes.
    
    Example
    ---------    
    For k=2 and inp=['A','B','C','D'], combination=['ACB','ACD','ABD','CBD']
    (each key has been combined with k other keys)
 
    Arguments
    ---------   
    inp_remaining: names of the genes whouse output-state we aim to be able 
                   to explain

    k:             integer, specifies the size of the combinations. k keys are 
                   combined with each key to create a combination. Note that k 
                   is the string size -1.

    Returns
    ---------
    combinations: list of strings, each string is a unqiue combination of genes 
                  from "inp_remaining". The order of the letters within a 
                  string is irrelevant. 
    """

    combinations=[]   
    for subset in itertools.combinations(inp_remaining, k):
        combinations.append("|".join(subset))
    return combinations





def dics(inp,outp,k,solved):
    """
    Auxiliary function that returns two dictionaries that are used in the 
    function "reveal".

    Arguments
    ---------    
    inp:    dict, each key is a string corresponding to a gene name and each 
            values is a list of integers corresponding to expression values in 
            different time points. The values correspond to input state.

    outp:   dict, each key is a string corresponding to a gene name and each 
            values is a list of integers corresponding to expression values in 
            different time points. The values correspond to outp state
   
    k:      integer, specifies the size of the combinations. k keys are  
            combined with each key to create a combination. Note that k is the 
            string size -1.

    solved: list of strings, each string is a gene. Genes within "solved"
            can be explained by other genes. The size of "solved" grows 
            gradually with each iteration of the function "reveal". The 
            algorithm converges when all genes are in "solved". "solved" allows 
            to ieterate in an efficient way, considering that it is not 
            required to carry the analysis for those genes that can aready be 
            explianed by previous iterations.

    Returns
    ---------
    d_inp:      dict, each key is a string corresponding to a gene name and  
                each values is the entropy value for that gene given its 
                expression  values in different time points.

    d_inp_outp: dict, each key is a string corresponding to a combination of 
                gene names and each values is the entropy value for such  
                combination given the expression values of the genes in 
                different time points.

    Note
    ---------
    Gene names within "inp" should be in capital letter, whereas gene 
    names within "outp" should be lower case. If the function "dics" is called 
    from "reveal", reveal will automatically do this transformation.

    """

    inp_remaining=[x for x in inp.keys() if x not in solved]
    comb_inp=combinations(inp_remaining,k)
    # inp_remaining: names of the genes in inp whouse output-state has not been
    # explained in previous iterations
    # comb_inp: all possible k+1-wise combinations of gene neames in inp

    d_inp={} #entropies of combinations of inp. keys are comb_inp, 
                #values are entropies 
    for combination in comb_inp:
        
        lista=[]
        for i in range(0,len(combination)):
            combi=combination.split("|")
            for gene in combi:
                lista.append(inp.values()[inp.keys().index(gene)]) 

        d_inp[combination] = entropy(*lista)#call with all elements of lista



    # Once d_inp is defined, we define d_inp_outp, which accounts for both 
    # inp and outp. This way we can access every entropy of class H(g,[GGG])
    # where lower letter sare outp-states and upper letters are inp-states
    
    single=outp.keys() 
    single=list(set(single)-set(solved)) #Do not carry for genes in "solved"
    inp_and_oups=[single,comb_inp]
    # inp_and_oups: List of lists. The 1st list is genes in outp, the 2nd is 
    # comb_inp
    comb_inp_and_oups=list(itertools.product(*inp_and_oups))
    comb_inp_and_oups=["|".join(x) for x in comb_inp_and_oups]
    # comb_inp_and_oups is the equivalent of "comb_inp", but considering
    # also outp (note that only one outp gene per combination is allowed). 

    d_inp_outp={}
    for combination in comb_inp_and_oups:
        lista=[]
        for i in range(0,len(combination)):#get the value
            combi=combination.split("|")
            for gene in combi:
                if gene in inp.keys():
                    lista.append(inp.values()[inp.keys().index(gene)]) 
                else:
                    lista.append(outp.values()[outp.keys().index(gene)])      

        d_inp_outp[combination] = entropy(*lista)#entropy of all lista elements

    return d_inp, d_inp_outp






def reveal(inp,outp,K=99999):
    """
    Given the expression values of a group of genes in different time points, 
    tries to find whether the expression of a gne can be explained by the 
    expression of another gene. The norm is: if the entropy of a combination
    of two genes is equal to the entropy of gene_2 and equal to "1", 
    then gene_1 (G1) can be explained by gene_2 (G2), as follows:
    if H(g1,G2) == H(G2), then G2->G1,
    where "g1" is the output state of "G1". *(end of script, for details)

    Arguments
    --------- 
    inp:  dict, each key is a string corresponding to a gene name and each 
          values is a list of integers corresponding to expression values in 
          different time points. The values correspond to input state

    outp: dict, each key is a string corresponding to a gene name and each 
          values is a list of integers corresponding to expression values in 
          different time points. The values correspond to outp state

    Returns
    --------- 
    revealed: dict, each key is a string corresponding to a gene name from 
              "inp" and each values is a list of strings where each string  
              represents a rule to explain the expression of the gene in key.  
              Thus, the relationship between two strings within a list can be 
              regarded as the logical "and". (i.e. 'c': ['A|C|B'] means that 
              the expression in output-state of gene C can be explained by the 
              rule: "A or C or B", where A, B and C are genes whose expression 
              in the inp-state is known).  

    Note
    --------- 
    Here H(a,[B,C]) (from [1]) is represented as H(a|B|C), 
    where gene names of genes in output-state are in lowercase. 
  


    """

    # For simplicity, gene names in "inp" are converted to upper case and gene 
    # names within "outp" are coinverted to lower case

    for key, value in inp.items():
        inp[key.upper()] = inp.pop(key) 
    for key, value in outp.items():
        outp[key.lower()] = outp.pop(key)


    # If no K is given, K will be equal to the number of genes+1
    # so maximum combination size will be equal to the number of genes in "inp" 
    # plus one gene from outp (size fig 5. of [1]). 
    if K==99999: #99999 is a symbolic number

        K=len(outp.keys())
        K+=1

    k=1
    solved=[] #genes whouse output could be explianed in previous iterations
    results=[] #genes that explain the otput of other genes

    while k<K+1: #k: current k; K: maximum k 
        dics_object=dics(inp,outp,k,solved)
        d_inp=dics_object[0]
        d_inp_outp=dics_object[1]

        iterate_over=outp.keys()
        iterate_over=list(set(iterate_over)-set(solved)) 
        # iterate_over: genes whose outp-states are still to be explained by 
        # the input of other genes

        regulations={}
        # regulations: dict, keys are genes, values are combinations of genes 
        # that regulate the gene in key.

        for gene_outp in iterate_over:
    
            # gene_outp: name of the gene whose output state we want to explain

            candidates=[ke for ke,v in d_inp_outp.items() if gene_outp in ke]
            # candidates: list of combinations of K+1 genes that may explain
            # gene_outp

            candidates=dict((ke, d_inp_outp[ke]) for ke in candidates)
            # candidates: dict, keys are possible combinations of K+1 genes that
            # that may explain gene_outp; values are the entropy of such combi
            # combination
        

        # 'a_B_B',for instance, will be renamed as 'B_B'
        # (only the capital letters: genes in "inp"). 
        # Thus,name is same as in d_inp, so we can search for intersections)
        # For instance, H(g1,G2) == H(G2), where the RHS is only genes in "inp"

            for key, value in candidates.items():
                key_noOutp=key.split("|")
                key_noOutp=key_noOutp[1:] 
                #key_noOutp: combination after removing gene_outp 
                key_noOutp="|".join(key_noOutp)
                candidates[key_noOutp] = candidates.pop(key)


            common = dict(set(candidates.iteritems()).\
            intersection(d_inp.iteritems()))  
            # common: dict, keys are combinations that are matches according to 
            # H(G)=H(g,G). Values are the entropy, so H(G)=H(g,G).


            # In case there is more than 1 solution
            solutions=[] #list of combinations whose input explain gene_outp  
            for keys, values in common.items():
                solutions.append(keys)
       
            regulations[gene_outp]=solutions
        
        for key, values in regulations.items():
            if values != []:
                solved.append(key)
                results.append(values)
        k+=1

    revealed={}
    for i in range(0,len(solved)):
        revealed[solved[i]]=results[i]

    #re-conver lower case to uppercase
    for key, value in revealed.items():
        revealed[key.upper()] = revealed.pop(key)
    

    return revealed


#create files 100 and 1000 time points.
#The files swilll be removed after Question 5
data= pd.read_csv('yeast_bin.tsv')

def some(x, n):
    return x.ix[random.sample(x.index, n)]
#https://stackoverflow.com/questions/15923826/random-row-\
#selection-in-pandas-dataframe

some_rows=some(data,100)
pd.DataFrame(some_rows).to_csv("file100_time_points", \
quoting=csv.QUOTE_NONE,index=False)

some_rows1000=some(data,1000)
pd.DataFrame(some_rows1000).to_csv("file1000_time_points", \
quoting=csv.QUOTE_NONE,index=False)




                                                                                                                                                                                                                
if __name__ == "__main__":


    # Data used in Liang et al, PSB 1998, in Fig. 1 and on pp. 23:
    
    inp = {'A':[0,0,0,0,1,1,1,1],'B':[0,0,1,1,0,0,1,1],'C':[0,1,0,1,0,1,0,1]}
    outp ={'A':[0,0,1,1,0,0,1,1],'B':[0,1,0,1,1,1,1,1],'C':[0,0,0,1,0,1,1,1]}



    #QUESTION 1
    print "QUESTION 1"
    print 
    print check_table(inp,outp)
    network = reveal(inp,outp)
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", \
              network[node]
    print
    print


    #QUESTION 2
    ready=read_tsv_file(open('yeast_bin.tsv'))
    inp=ready[0]
    outp=ready[1]



    #QUESTION 3
    print "QUESTION 3"
    print
    print "for k=1:" 
    ready=read_tsv_file(open('yeast_bin.tsv'))
    inp=ready[0]
    outp=ready[1]
    network = reveal(inp,outp,1)
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", \
              network[node]
    print("--- %s seconds ---" % round((time.time() - start_time),4))
    print 
    print "for k=2:" 
    print "This takes around 150 seconds"
    ready=read_tsv_file(open('yeast_bin.tsv'))
    inp=ready[0]
    outp=ready[1]
    network = reveal(inp,outp,2)
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", \
              network[node]
    print("--- %s seconds ---" % round((time.time() - start_time),4))
    print 
##    print "for k=3:" 
##    ready=read_tsv_file(open('yeast_bin.tsv'))
##    inp=ready[0]
##    outp=ready[1]
##    network = reveal(inp,outp,3)
##    print
##    for node in network.keys():
##        print "The expression of node", node, "is best explained by nodes", \
##              network[node]
##    print("--- %s seconds ---" % round((time.time() - start_time),4))
##    print 




  #QUESTION 5
    print "QUESTION 5"
    print
    print "for 1000 random points and k2:" 
    ready=read_tsv_file(open('file1000_time_points'))

    inp=ready[0]
    outp=ready[1]
    network = reveal(inp,outp,2)
    print network
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", \
              network[node]
    if network=={}:
        print "It was not possible to explain the expression of any gene"

    print 
    print "for 100 random points and k3:" 
    ready=read_tsv_file(open('file100_time_points'))

    inp=ready[0]
    outp=ready[1]
    network = reveal(inp,outp,3)
    print network
    print
    for node in network.keys():
        print "The expression of node", node, "is best explained by nodes", \
              network[node]
    print 
    os.remove("file100_time_points")
    os.remove("file1000_time_points")




"""
Notes:

    *We search for H(g1,G2) == H(G2), as in the equation just after figure5 in 
    [1]. Thus, the fucntion "mutual_information" is not used in the the reveal 
    algorithm. 


References: 

    [1] S. Liang, S. Fuhrman, R. Somogyi, REVEAL
        http://psb.stanford.edu/psb-online/proceedings/psb98/liang.pdf
"""



