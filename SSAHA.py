#!/usr/bin/env python
"""
Author: Fernando Bueno Gutierrez
Implementation of the SSAHA algorithm to find high-scoring short local alignments(indexing)
"""

# import statements
from sys import argv
import time
# implement your functions here
def generate_hash_table(seqlist, k):
    """
    :param seqlist: list containing DNA sequences
    :param kmer_length: integer specifying the length of the k-mers generated
    by the function
    :return: hash table containing the non-overlaping k-mers as keys and the
    positions as values.
    """
    hash_table = {}
    for index, seq in enumerate(seqlist):
        for kmer in range(0, len(seq)-k, k):
            kmerseq = seq[kmer:kmer + k]

            if kmerseq in hash_table:
                hash_table[kmerseq].append((index+1, kmer + 1))

            else:
                hash_table[kmerseq] = []
                hash_table[kmerseq].append((index+1, kmer + 1))

    return hash_table


def get_hits(seq, hash_table, k):
    """
    :param k: word size (should be the same as the one used for the hash table)
    :param seq: Query DNA sequence.
    :param hash_table: Hash table containing positions for certain k-mers.
    :return: list of hits of the query sequence in the database ("hash table")
    sorted by index and shift
    """
    hits = []
    for idx in range(0, (len(seq) - k + 1)):
        kmer = seq[idx:idx + k]
        if kmer in hash_table:
            for position in hash_table[kmer]:
                pos = []
                pos.append(position[0])
                pos.append(position[1] - idx)
                pos.append(position[1])
                hits.append(pos)
    hits = sorted(hits)
    return hits


def get_longest_match(hits):
    """
    :param hits: list containing the sorted hits of a query sequence in a
    database of sequences. Each hit should contain the sequence index, its
    shift and the offset.
    :return: 3 objects. Matches is a list containing the sequence of all the
    longest matches found in the database. seq_coords is a list containing the
    index and start and end position of the matches in the sequence database.
    query_coords is a list containing the start and end position of the longest
    match in the query sequence.
    """
    matchcount = []
    counter = 0
    for hit in range(0, len(hits) - 1):
        if hits[hit][0] == hits[hit + 1][0] and hits[hit][1] == hits[hit + 1][
            1]:
            counter += 1
            matchcount.append(counter)
        else:
            counter = 0
            matchcount.append(0)
    maxmatch = max(matchcount)
    idxs = [idx for idx, count in enumerate(matchcount) if count == maxmatch]
    seq_coords = []
    query_coords = []
    matches = []
    for idx in range(0, len(idxs)):
        whichseq = hits[idxs[idx]][0] - 1
        match_start = hits[idxs[idx] - maxmatch + 1][2] - 1
        match_end = hits[idxs[idx] + 1][2]
        query_start = hits[idxs[idx] - maxmatch + 1][2] - \
                      hits[idxs[idx] - maxmatch + 1][1]
        query_end = match_end - match_start + query_start

        seq_coords.append([whichseq, match_start, match_end])
        query_coords.append([query_start, query_end])
        match = query[query_start:query_end + 1]
        matches.append(match)
    return matches, seq_coords, query_coords


def print_alignment(query_seq, db_seqs, query_coords, seq_coords):
    """
    :param query_seq: DNA sequence that is wanted to compare with a database of
    sequences.
    :param db_seqs: list of sequences that are used as a database to compare
    the query with.
    :param query_coords: Start and end position of the longest match in the
    query sequence
    :param seq_coords: Index, start and end position of the longest match in
    the sequences.
    :return: Top, middle and bottom, three strings that when printed show the
    alignment between the query sequence and the longest match in the database
    sequences.
    """

    top = []
    bottom = []
    middle = []

    empty = seq_coords[0][1] - query_coords[0][0]
    if empty >= 0:
        top.append(db_seqs[seq_coords[0][0]])
        for i in range(empty):
            bottom.append(" ")
        bottom.append(query_seq)
    if empty < 0:
        bottom.append(query_seq)
        for i in range(-empty):
            top.append(" ")
        top.append(db_seqs[seq_coords[0][0]])
    for i in range (empty + len(query_seq), len(db_seqs[seq_coords[0][0]])):
        bottom.append(" ")
    top = "".join(top)
    bottom = "".join(bottom)

    for i in range(0,len(top)):
      if top[i] == bottom[i]:
          middle.append('|')

      elif top[i] == " " or bottom[i] == " ":
          middle.append(" ")

      else:
          middle.append('.')

    middle = "".join(middle)
    print top
    print middle
    print bottom
    return top, middle, bottom

def fasta_parse(filename):
    """
    :param filename: filename (or if in different than working directory, the
    complete path) of the fasta file that needs to be parsed.
    :return: two lists, seqs, which contains the DNA sequences in the
    fasta file, and names, which contains the names of the sequences.
    """
    fastaFile = open(filename, 'r')
    seq = []
    seqs = []
    names = []
    for line in fastaFile:
        if line[0] == ">":
            line = line.rstrip()
            names.append(line[1:])
            if seq != []:
                seq = "".join(seq)
                seqs.append(seq)
                seq = []
        else:
            line = line.rstrip()
            seq.append(line)
    seq = "".join(seq)
    seqs.append(seq)
    fastaFile.close()
    return names, seqs

def complement_seq(fwd_seqs):
    """
    Function that takes a DNA sequence as input and outputs its complementary.
    :param fwd_seqs: list of DNA sequence in the forward sense.
    :return: rvs_seqs: list of the same DNA sequences in the reverse sense.
    """
    rvs_seqs = []
    for seq_idx, seq in enumerate(fwd_seqs):    #for every sequence in list

        currentSeq = []
        rvs_seqs.append(currentSeq)
        for char_idx, char in enumerate(seq): #for every character in sequence

            if seq[char_idx] == "A":
                currentSeq.append("T")

            elif seq[char_idx] == "T":
                currentSeq.append("A")

            elif seq[char_idx] == "G":
                currentSeq.append("C")

            elif seq[char_idx] == "C":
                currentSeq.append("G")

        rvs_seqs[seq_idx] = currentSeq  #update list with currentSeq
        rvs_seqs[seq_idx] = rvs_seqs[seq_idx][::-1] #reverse the generated seq
        rvs_seqs[seq_idx] = "".join(rvs_seqs[seq_idx]) # convert list to string
    return rvs_seqs

if __name__ == "__main__":

    # the code below should produce the results necessary to answer the questions
    # in other words, if we run your code, we should see the data that you used
    # to answer the questions
    start_time = time.time()


    #QUESTION 1:

    print "\n QUESTION 1 \n   "
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    seqs = [s1, s2, s3]



    hash = generate_hash_table(seqs, 2)
    print "The generated hash table is: \n"
    for key in sorted(hash.keys()):
        print "%s: %s" % (key, hash[key])

    #QUESTION 2:

    print "\n QUESTION 2 \n"

    query = 'TGCAACAT'
    hits = get_hits(query, hash, 2)
    print("There is a total of " + str(len(hits)) + " hits")
    print "The first and last sorted hits are the following:"
    print (hits)[0], (hits)[-1]

    #QUESTION 3:

    print "\n QUESTION 3 \n"
    maxmatch, seq_coords, query_coords = get_longest_match(hits)

    print "The longest alignment is the following: \n"
    print_alignment(query, seqs, query_coords, seq_coords)

    #QUESTION 4:

    print "\n QUESTION 4 \n"

    #fastaFile = argv[1]
    fastaFile = "TAIR10.fasta.txt"

    names, seqs = fasta_parse(fastaFile)
    total_length = 0
    for i in range (0,len(names)):
        print ("The length of the sequence of " + names[i] + " is of " +
                str(len(seqs[i])) + " bp.")
        total_length += len(seqs[i])
    print ("\nThe total length of the provided Arabidopsis thaliana genome is"+
            " of " + str(total_length) + " bp. \nThere is a total of " +
           str(len(names)) + " sequences.")

    print (time.time() - start_time)


    #QUESTION 5:

    print "\n QUESTION 5 \n"

    hash = generate_hash_table(seqs, 15)
    hash_length = len(hash)
    print ("There is a total of " + str(hash_length) + " unique 15-mers in " +
           "the provided \nArabidopsis thaliana genome.")

    #QUESTION 6:

    print "\n QUESTION 6 \n"

    fastaQuery = "athal_query.fasta.txt"
    #fastaQuery = argv[2]
    query_name, query_seqs = fasta_parse(fastaQuery)

    n= -1
    for query in query_seqs:
        n += +1
        hits = get_hits(query, hash, 15)
        maxmatch, seq_coords, query_coords = get_longest_match(hits)
        for match in maxmatch:
            if len(match) > 15:
                print("The longest match for query " + query_name[n] +
                      " is found in Chromosome " + str(seq_coords[0][0] + 1) +
                      " from bp " + str(seq_coords[0][1]) + " to bp " +
                      str(seq_coords[0][2]) + " so a total stretch of " +
                      str(seq_coords[0][2] - seq_coords[0][1]) + " bp.")
                print("\nThe matching sequence is the following:")
                print(match)
        if max(map(len,maxmatch)) < 15:
            print("\nNone of the 15-mers of " + query_name[n]+
                  " match with the 15-mers in the hash table")


    #OPTIONAL PART:

    query_seqs = complement_seq(query_seqs)

    n= -1
    for query in query_seqs:
        n += +1
        hits = get_hits(query, hash, 15)
        maxmatch, seq_coords, query_coords = get_longest_match(hits)
        for match in maxmatch:
            if len(match) > 15:
                print("The longest match for query " + query_name[n] +
                      " in the reversed direction is found in Chromosome " +
                      str(seq_coords[0][0] + 1) +
                      " from bp " + str(seq_coords[0][1]) + " to bp " +
                      str(seq_coords[0][2]) + " so a total stretch of " +
                      str(seq_coords[0][2] - seq_coords[0][1]) + " bp.")
                print("\nThe matching sequence is the following:")
                print(match)
        if max(map(len,maxmatch)) < 15:
            print("\nNone of the 15-mers of " + query_name[n]+
                  " match with the 15-mers in the hash table")





    print (time.time() - start_time)