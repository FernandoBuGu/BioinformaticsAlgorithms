#!/usr/bin/env python

"""
Author: Fernando Bueno Gutierrez
Student nr:890605143090
Script to: Train a profile Hidden Markov profile_states 
           (profile HMM) for multiple sequence alignment
"""

from sys import argv
from random import random
import math
import decimal
import collections
import numpy as np
import random
import pandas as pd 
import numpy as np


# Background amino acid probabilities
pa = { 'A':0.074, 'C':0.025, 'D':0.054, 'E':0.054, 'F':0.047, 'G':0.074, 'H':0.026, 'I':0.068, 'L':0.099, 'K':0.058, 'M':0.025, 'N':0.045, 'P':0.039, 'Q':0.034, 'R':0.052, 'S':0.057, 'T':0.051, 'V':0.073, 'W':0.013, 'Y':0.034 } 


class HMM():
    """HMM object to store an HMM model

    This object is designed to keep track of all HMM states, emissions, and 
    transitions. It may be used in your implementation, but may also be ignored,
    and replaced by a data structure of choice
    """
    # Emission probabilities for the match and insert states
    e_m   = []; e_i   = pa; 
    
    # Transition probabilities from/to matches, inserts and deletions
    t_mm  = []; t_mi  = []; t_md = [];
    t_im  = []; t_ii  = []
    t_dm  = []; t_dd  = []; 
    
    def __init__(self,nmatches):
        """Initialize HMM object with number of match states
        
        nmatches: int, number of match states
        """
    
        self.nmatches = nmatches
        
        self.e_m   = [dict(pa) for i in range(0,nmatches)]
        for i in range(0,nmatches):
            for j in pa.keys():
                self.e_m[i][j] = 0.0
        self.e_i   = pa;

        self.t_mm  = [0.0 for i in range(0,nmatches+1)]
        self.t_mi  = [0.0 for i in range(0,nmatches+1)]
        self.t_im  = [0.0 for i in range(0,nmatches+1)]
        self.t_ii  = [0.0 for i in range(0,nmatches+1)]
        self.t_md  = [0.0 for i in range(0,nmatches+1)]
        self.t_dm  = [0.0 for i in range(0,nmatches+1)]
        self.t_dd  = [0.0 for i in range(0,nmatches+1)]

def sample(events):
    """Return a key from dict based on the probabilities 

    events: dict of {key: probability}, sum of probabilities
    should be 1.0. 
    """
    k = events.keys()
    cum = [0 for i in k]

    cum[0] = events[k[0]]
    for i in range(1,len(events)):
        cum[i] = cum[i-1] + events[k[i]]
    # Should not be necessary, but for safety
    cum[len(cum)-1] = 1.0

    r = random()
    pick = ''
    i = 0
    while (pick == '') and (i < len(cum)):
        if r < cum[i]:
            pick = k[i]
        i = i + 1
    return pick

def parse_fasta(fasta_file):
    """
    Return a dictionary corresponding to a parsed fasta file 

    Arguments
    ---------
   	fasta_file:     string, filename of the fasta file to be parsed

    Returns
    ---------
    parsed_seqs:    dictionary, parsed fasta file. Each key is the header
                    of a sequence as it apperas in the fasta, and each value
                    is the corresponding sequence.
    """

    with open(fasta_file) as fasta_file:
        fasta_list = fasta_file.read().splitlines()

        parsed_seqs = {}
        for line in fasta_list:
            if line.startswith(">"):
                label = line[1:]
                parsed_seqs[label] = ""

            else:
                parsed_seqs[label] += line

    return parsed_seqs


def count_match(parsed_seqs):
    """
    Given a parsed fasta file, returns an integer corresponding to the number 
    of match states of a sequence alignment in a profile HMM.

    Return a list of strings corresponding to the states of a profile HMM.

    Arguments
    ---------
   	parsed_seqs:    dictionary, parsed fasta file. Each key is the header
                    of a sequence as it apperas in the fasta, and each value
                    is the corresponding sequence. 

    Returns
    ---------
    count_match_seq_int:      integer, number of match states of the HMM 

    profile_states:           list of strings corresponding to the states of 
                              the profile HMM. Each string is one position in 
                              the multiple seqeunce alignment. Four different 
                              states are considered: start ('i'),match ('M'),
                              non-match ("-"), and end ("!"). 'M' is assigned 
                              to those positions for which at least half of the
                              sequences have an amino-acid.
    """
    count_match_per_position= [0] * len(parsed_seqs.values()[0]) 
    #for each position we count the number of sequences for which there is an 
    #amino-acid in that position
    for key, value in parsed_seqs.items():
        for idx, position in enumerate(value):
            if position != "-":
                count_match_per_position[idx]+=1

    threshold=len(file_parsed)//2+1

    count_match_seq=[] 

    count_match_seq_int=0
    #count_match_seq_int counts the number of positions for which we assign 'M'

    for position in count_match_per_position:
        if position>=threshold:
            count_match_seq.append('M') #'M': Match state
            count_match_seq_int+=1
        else:
            count_match_seq.append('-') #'-': Non-match state
    
    start=['i']
    end=['!']
    profile_states= sum([start, count_match_seq, end], [])
    #profile_states will be used in other functions
    
    return count_match_seq_int, profile_states
		


def list_per_position(list_seqs):
    """
    Given a list of sequences, returns a list of lists where each inner list 
    corresponds to a position in the sequence alignment and it contains a 
    string for each aa found in that position.

    Arguments
    ---------
    list_seqs: list of strings, each element is a sequence, each letter an aa

    Returns
    ---------
    list_seqs_or: list of lists, each inner list corresponds to a position in 
                  the sequence alignment and it contains a string for each aa 
                  found in that position.
    """

    list_seqs_or = [[] for i in range(len(list_seqs[0]))]  
    #"list_seqs_or" will have one element per position. Assumming that all 
    #sequences have equal length, the first element is used to define the 
    #number of elements of "list_seqs_or"

    for seq in list_seqs:
        for i in range(len(seq)):
            list_seqs_or[i].append(seq[i])
        
    return list_seqs_or



def return_dictionary(list_position):
    """
    Given a list of strings corresponding to the aa found in a given
    position in the multiple sequence alignment, returns a dictionary where the
    keys are the aa and the values are the emissions of that aa for that 
    particular position.

    Arguments
    ---------
    list_position: list of strings, the strings are aa found in a given 
                   position. It is analogous to the inner lists of 
                   "list_seqs_or" 

    Returns
    ---------
    newD: dictionary, each keys is an aa and each value is the emission of 
          that aa(i.e the proportion of sequences that have that aa in that 
          position. 

    """

    d={}  
    total=0 #total: the number of sequences that have any aa in that position
    for letter in list_position:
        if letter != "-":
            if letter in d:
                d[letter] += 1
            else:
                d[letter] = 1
            total+=1

    newD = {key:give_proportion(value,total) for key, value in d.items()}

    return newD



def give_proportion(value,total):
    """
    General function. Given a count in one class and the total count in all 
    classes, returns the proportion of the first class within all classes.

    Arguments
    ---------
    value: integer, count of one class
    total: integer, total count of all classes
    
    Returns
    ---------
    proportion_rounded: float with two decimals, proportion of the first class 
                        within all classes.
    """ 

    value=float(value)
    proportion=value/total
    proportion_rounded=decimal.Decimal(proportion)
    proportion_rounded=round(proportion_rounded,2)
    return proportion_rounded


def return_emisions(list_seqs,profile_states):  
    """
    Given a list of sequences and a list with the states of the HMM, returns a
    list of dictionaries, where each dictionary corresponds to the emissions of
    each match position in the multiple sequences alignment. 

    Arguments
    ---------
    list_seqs: list of strings, each element is a sequence, each letter an aa
    profile_states: 

    profile_states: list of strings corresponding to the states of 
                    the profile HMM. Each string is one position in 
                    the multiple seqeunce alignment. Four different 
                    states are considered: start ('i'),match ('M'),
                    non-match ("-"), and end ("!"). 'M' is assigned 
                    to those positions for which at least half of the
                    sequences have an amino-acid.

    Returns
    ---------
    emissions: list of dictionaries, where each dictionary corresponds to the 
               emissions of each match position in the multiple sequences 
               alignment. 
    """
    
    list_seqs_or=list_per_position(list_seqs)

    emissions=[]
    for i in range(len(profile_states)):
        if profile_states[i] == 'M':   
        #Non-match states will have other emissions ("pa")
            dict_position=return_dictionary(list_seqs_or[i-1]) 
            #"dict_position": dictionary of a given position. "-1" because 
            #"list_seqs_or" starts in position "1" with respect to 
            #"profile_states" (since profile_states starts in "i")
            emissions.append(dict_position)

    return emissions




def transform_string(seq,profile_states):
    """
    Given a sequence and the profile states, return a list of strings 
    corresponding to the states of that sequence with respect to the 
    profile states.

    Arguments
    ---------
    seq:            string, sequence, each letter is an aa in the sequence
    profile_states: list of strings corresponding to the states of 
                    the profile HMM. Each string is one position in 
                    the multiple seqeunce alignment. Four different 
                    states are considered: start ('i'),match ('M'),
                    non-match ("-"), and end ("!"). 'M' is assigned 
                    to those positions for which at least half of the
                    sequences have an amino-acid.

    Returns
    ---------
    seq_states: list of strings corresponding to the states of a particular 
                sequence with respect to the profile HMM. Each string is one 
                position in the HMM. The number in the string is the number of 
                the state in the profile HMM, and the non-numeric characters
                are: start ('i'),match ('M'), non-match ("-"), and end ("!"). 
                'M' is assigned to those positions for which at least half of 
                the sequences have an amino-acid. (i.e. if a sequence has a 
                '5D', it means that the sequence has a deleteion in the state 5
                of the profile HMM). 
    """


    seq='i'+seq+'!'
    seq_states=[]
    S=0
    

    for i in range(len(profile_states)):

        if profile_states[i] == 'i' or profile_states[i] == '!':
            if  profile_states[i] == 'i':           
                seq_states.append('00i')
            if profile_states[i] == '!':
                seq_states.append(profile_states[i])

        else:
            if profile_states[i] == 'M' and seq[i] != '-':
                S+=1
                app=str(S)+'M'

            if profile_states[i]!='M' and seq[i] !='-':
                app=str(S)+'I'

            if profile_states[i]!='M' and seq[i] =='-':
                continue

            if profile_states[i]=='M' and seq[i] =='-':
                S+=1
                app=str(S)+'D'
            seq_states.append(app)

    return seq_states      


def transitions_states(seq):
    """
    Given the states, returns a list of strings with the transitions between
    states. 

    Arguments
    ---------
    seq:   List of str, each string represents a state. (i.e. 'I': start, 
          '1M': first match, '!': end)

    Returns
    ---------
    list_seqs: List of str, each string represents the transition form one 
               state to another. The original state apperad before "-", whereas 
               the "next" state appears after.  
               (i.e. ['00i-M', '01M-M', '02M-I'] means thatthe sequence goes 
               from start-state to the first match, then to the second match, 
               and then one insertion). 
    """

    list_seqs=[]
    for i in range(len(seq)-1):
        if seq[i+1] != '!':
            paste=seq[i]+"-"+seq[i+1][-1]
        else:
            paste=seq[i]+"-"+'!'
        list_seqs.append(paste)

    return list_seqs



def once_transformed(seq_states):

    """
    Given the transitions between states, returns a list of strings with the 
    count of sequences that such transitions. Thus, one count per transition. 

    Arguments
    ---------
    seq_states:  List of str, each string represents the transition form one  
                 state to another. The original state apperad before "-", 
                 whereas the "new" state appears after.  
                 (i.e. ['00i-M', '01M-M', '02M-I'] means that the sequence goes 
                 from start-state to the first match, then to the second match, 
                 and then one insertion).

    Returns
    ---------
    DT:   Dict, each key is a string corresponding to the transtion between
          one state and the next. Each value is an integer corresponding to the
          number of sequences in the data that have such transition. 
    """
    

    DT={} #DT: dict transitions
    for SEQ in seq_states:
        for pos in SEQ:
            if pos not in DT:
                DT[pos] = 1
            else:
                DT[pos] += 1
    return DT


def transitions_in_proportions(DT):

    """
    Given the count of transitions between states, returns a list of strings 
    with the transition probabilities. Thus, one probability per transition. 

    Arguments
    ---------
    DT:   Dict, each key is a string corresponding to the transtion between
          one state and the next. Each value is an integer corresponding to the
          number of sequences in the data that have such transition. 

    Returns
    ---------
    DT:   Dict, each key is a string corresponding to the transtion between
          one state and the next. Each value is a float corresponding to the
          transition probability. 
    """

    DT_auxiliar={} #DT_auxiliar: dict transitions auxiliar
    for key, values in DT.items():
        fromm=key.split("-")[0]                     
        if fromm not in DT_auxiliar:
            DT_auxiliar[fromm] = values
        else:
            DT_auxiliar[fromm] += values


    for i in range(len(DT)):

        fromm=DT.keys()[i]

        fromm2=fromm.split("-")[0]
 
        this=DT[fromm]
        total=DT_auxiliar[fromm2]

        DT[fromm]=give_proportion(this,total)
    return DT



def return_transitions(seqs,profile_states):
    
    """
    Given a string corresponding to the sequence of interest and a list of
    strings corresponding to the profile states of the HMM, returns the
    transition probabilities. 

    Arguments
    ---------
    seqs:           Str, sequence of interest. Each character is an amino-acid

    profile_states: List of str, each string represents the transition from one 
                    state to another. The original state apperad before "-", 
                    whereas the "new" state appears after. (i.e. 
                    ['00i-M', '01M-M', '02M-I'] means that the sequence goes 
                    from start-state to the first match, then to the second
                    match, and then one insertion).

    Returns
    ---------
            DT2:  Dict, each key is a string corresponding to the transtion 
                  between one state and the next. Each value is a float 
                  corresponding to the transition probability. 

    d_ordered:    collections.OrderedDict, as d2, but sorted from the start 
                  state to the end. Allows for eassier visualization. 
    """


    seqs_transformed=[]
    for seq in seqs: 
        seqs_transformed.append(transform_string(seq,profile_states))

    seqs_transformed2=[]
    for seq in seqs_transformed: 
        seqs_transformed2.append(transitions_states(seq))

    DT = once_transformed(seqs_transformed2)

    DT2 = transitions_in_proportions(DT)

    DT_sorted = collections.OrderedDict(sorted(DT2.items()))

    return DT2, DT_sorted


def sample_from_HMM(candidates):
    """
    Given the transition probabilities, samples the "next" state.
    
    Arguments
    ---------
    candidates: dict, each key is one of the possible transition from the same 
                state, each value is the transition probabilities.

    Returns
    ---------
    random.choice(candidates): char, corresponding to the sampled transition

    Example
    --------
    The for input {'3M-M': 0.88, '3M-D': 0.06, '3M-I': 0.06}, it is more
    likelly that '3M-M' is sampled, but 3M-D or 3M-I coudl also be sampled
    """

    lis = candidates.keys()
    times= candidates.values()
    times = [ int(100*x) for x in times ]
    times=tuple(times)
    candidates = sum(([x]*y for x,y in zip(lis, times)),[])
    return random.choice(candidates)



def sample_transition(number,state,dic_transitions):
    """
    Given the position we are interested in, the state at such position for
    a seqeunce of interest (M,I or D) and the transitions probabilities of the 
    HMM, returns a sampled transition for the input position ans state.

    Arguments
    ---------
    number:             int, the position in the HMM 

    state:              str (M,I or D), the state of the sequence of interest 
                        in the position "number"
           
    dic_transitions:    Dict, each key is a string corresponding to the  
                        transtion between one state and the next. Each value is  
                        a float corresponding to the transition probability. 

    Returns
    ---------
    sampled_transition: str, the sampled transition given the input
    """

    key=str(number)+state
    def slicedict(d, s):        #https://stackoverflow.com/questions/4558983/slicing-a-dictionary-by-keys-that-start-with-a-certain-string
        return {k:v for k,v in d.iteritems() if k.startswith(s)}
    candidates= slicedict(dic_transitions, key)
    sampled_transition=sample_from_HMM(candidates)
    return sampled_transition


def sample_aa(location,dic,insertion=False):
    
    """
    Given the location in the HMM map of states, the emission probabilities and 
    whether the state (in the aligned sequence we are interested in) is an 
    inserction, returns a sampled aa.

    Arguments
    ---------
    location:   int, location in the HMM map of states

    dic:        list of dict, where each dictionary corresponds to the 
                emissions of each match position in the multiple sequences
                alignment.

    insertion:  logical that is true if for the sequence of interset in the
                aa number "location" there is an insertion. In such case, 
                location will not increase  

    Returns
    ---------
    choice: int, sampled aa

    """
    if insertion==False:
        candidates=dic[location-1]
    else:
        candidates=dic
    choice=sample_from_HMM(candidates)
    return choice


def sample_sequence_from_HMM(dic_transitions,dic_emissions,pa):

    """
    Given the transitiions probabilities and the emission probabilities for
    the states and teh insertion states, samples a sequence from the HMM.

    Arguments
    ---------
    dic_transitions:    Dict, each key is a string corresponding to the  
                        transtion between one state and the next. Each value is  
                        a float corresponding to the transition probability. 

    emissions:          list of dict, where each dictionary corresponds to the 
                        emissions of each match position in the multiple
                        sequences alignment.

    pa:                 Dict, each key is a string corresponding to an aa. 
                        Each value is the background frequency of such aa. 

    Returns
    ---------
    seq:   str, each character is either an aa or a deletion ("-") from the 
           sampled sequence  
    """

    seq=[]
    location=0
    start='00i'
    next=sample_transition('00','i',dic_transitions)
    next=next.split("-")[1]
    if next == 'M' or next == 'D':
        location+=1
    while next != '!':
        if next == 'M':
            new_aa=sample_aa(location,dic_emissions)
            seq.append(new_aa)
        if next == 'I':
            new_aa=sample_aa(location,pa,insertion=True)
            seq.append(new_aa)
        if next == 'D':
            seq.append("-")

        next=sample_transition(location,next,dic_transitions)
        next=next.split("-")[1]
        if next == 'M' or next == 'D':
            location+=1
    seq="".join(seq)
    return seq


def fill_with_backgrounds(count,pa):
    """
    Returns list of dictionaries corresponding to the emissions using 
    background frequencies

    Arguments
    ---------  
    count: integer, indicates the number of match states

    pa: dict, background frequencies

    Returns
    ---------
    filles_with_backgrounds: list of dict, emissions using background 
                             frequencies. Each dictionary contains the emission 
                             probabilities of a match state, respectively. Each 
                             key in dictionary is an aa, and each value is the
                             probability of finding such aa in that position. 
    """
    
    filles_with_backgrounds=[]
    for i in range(0,count):
        filles_with_backgrounds.append(pa)
    return filles_with_backgrounds

    
           
if __name__ == "__main__":

    # implement main code here
    infile = 'test.fasta'

#QUESTION 1
    print "QUESTION 1"
    print
    filename = "test.fasta"
    file_opened = open(filename, "r")
    file_readed = file_opened.readlines()
    file_parsed = parse_fasta(filename)

    count=count_match(file_parsed)
    print "The number of match states needed is:"
    print count[0]
    print



#    #QUESTION 2
    print "QUESTION 2"
    print
    print "emissions in the insert states"
    print pa
    print
    print "emissions"
    print return_emisions(file_parsed.values(),count[1])
    print
    print "trasnsitions"
    print return_transitions(file_parsed.values(),count[1])[0]
    print


#    #QUESTION 3
    print "QUESTION 3"
    print
    print "emissions in PSSM format:"
    print
    dic_emissions = return_emisions(file_parsed.values(),count[1])
    df=pd.DataFrame(dic_emissions)
    df[np.isnan(df)] = 0
    print df
    print



#    #QUESTION 4
    print "QUESTION 4"
    print
    dic_emissions = return_emisions(file_parsed.values(),count[1])
    dic_transitions = return_transitions(file_parsed.values(),count[1])[0]
    print "samples generated with the HMM:"
    print
    for i in range(0,10):
        print sample_sequence_from_HMM(dic_transitions,dic_emissions,pa)
    print


#    #QUESTION 5

    filled_with_backgrounds=fill_with_backgrounds(count[0],pa)

    print "QUESTION 5"
    print
    print "emissions in the insert states"
    print pa
    print
    print "emissions using the background frecuencies"
    print filled_with_backgrounds
    print
    print "trasnsitions"
    print return_transitions(file_parsed.values(),count[1])[0]
    print

    print
    print "emissions in PSSM format:"
    print
    df=pd.DataFrame(filled_with_backgrounds)
    df[np.isnan(df)] = 0
    print df
    print

    print "samples generated with the HMM:"
    print
    dic_transitions = return_transitions(file_parsed.values(),count[1])[0]
    for i in range(0,10):
        print sample_sequence_from_HMM(dic_transitions,filled_with_backgrounds\
                ,pa)
    print

#    #QUESTION 6

    print "QUESTION 6"
    print
    filename = "test_large.fasta"
    file_opened = open(filename, "r")
    file_readed = file_opened.readlines()
    file_parsed = parse_fasta(filename)
    count=count_match(file_parsed)
    filled_with_backgrounds=fill_with_backgrounds(count[0],pa)
    print "emissions in the insert states"
    print pa
    print
    print "emissions"
    print return_emisions(file_parsed.values(),count[1])
    print
    print "trasnsitions"
    print return_transitions(file_parsed.values(),count[1])[0]
    print
    dic_emissions = return_emisions(file_parsed.values(),count[1])
    dic_transitions = return_transitions(file_parsed.values(),count[1])[0]
    print
    print "emissions in PSSM format:"
    print
    df=pd.DataFrame(dic_emissions)
    df[np.isnan(df)] = 0
    print df
    print
    print "samples generated with the HMM:"
    print
    for i in range(0,10):
        print sample_sequence_from_HMM(dic_transitions,dic_emissions,pa)
    print
    print "Now we see what happens when we change to background frecuencies:"
    print

    print "emissions in the insert states"
    print pa
    print
    print "emissions using the background frecuencies"
    print filled_with_backgrounds
    print
    print "trasnsitions"
    print return_transitions(file_parsed.values(),count[1])[0]
    print
    print "emissions in PSSM format:"
    print
    df=pd.DataFrame(filled_with_backgrounds)
    df[np.isnan(df)] = 0
    print df
    print
    print "samples generated with the HMM:"
    print
    dic_transitions = return_transitions(file_parsed.values(),count[1])[0]
    for i in range(0,10):
        print sample_sequence_from_HMM(dic_transitions,filled_with_backgrounds\
                ,pa)
    print





    

    













