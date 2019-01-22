#!/usr/local/bin/python



"""
#Author: Fernando Bueno

#Description: Implementation of Needleman-Wunsch algorithm (Dynamic Programming) for global sequence alignment
	#Uses BLOSUM 62 substitution matrix
	#Control end-gap penalties
	#Computes: percentage of identity between two aligned sequences
"""


library("Biostrings")
data(BLOSUM62)


intit_matrix = function(gap_penalty,seq1,seq2) 
{
    """
    #Given two sequences of a.a. and a gap penalty at the end of sequences, it returns the initial matrix of Needleman-Wunsch
    """

	seq1_ch<-c("0","0")
	seq1_ch<-append(seq1_ch,unlist(strsplit(seq1,"")))
	seq2_ch<-c("0","0")
	seq2_ch<-append(seq2_ch,unlist(strsplit(seq2,"")))
	D<-data.frame(matrix(NA, nrow = length(seq2_ch), ncol = length(seq1_ch)))
	D[1,]<-c(seq1_ch)
	D[,1]<-c(seq2_ch)
	w=0
	for(i in 2:length(seq1_ch)) {
		if(exists("w"))
		{
			D[2,i]<-w
		} else {
			D[2,i]<-0
		}
		w=w-gap_penalty
	}
	w=0
	for(j in 2:length(seq2_ch)) {
		if(exists("w"))
		{
			D[j,2]<-w
		} else {
			D[j,2]<-0
		}
		w=w-gap_penalty
	}
	return(D)
}



fill_matrix = function(D,gap_penalty)
{
	#Given the initial matrix of Needleman-Wunsch and the gap penalty it fills the Needleman-Wunsch for those sequences

	seq1_ch<-c("0","0")
	seq1_ch<-append(seq1_ch,unlist(strsplit(seq1,"")))
	seq2_ch<-c("0","0")
	seq2_ch<-append(seq2_ch,unlist(strsplit(seq2,"")))
	for(j in 3:length(seq2_ch)) {
		for(i in 3:length(seq1_ch)) {
		H=D[1,i]			
		V=D[j,1]
		S=BLOSUM62[colnames(BLOSUM62)==as.character(V),rownames(BLOSUM62)==as.character(H)]
		a=as.numeric(D[j-1,i-1])+S
		b=as.numeric(D[j-1,i])-gap_penalty
		c=as.numeric(D[j,i-1])-gap_penalty
		D[j,i]<-max(a,b,c)
		}
	}
	return(D)
}





intit_matrix_back = function(seq1,seq2) 
{
	#Given two sequences of a.a., it returns initial traceback-matrix of Needleman-Wunsch
	seq1_ch<-c("0","0")
	seq1_ch<-append(seq1_ch,unlist(strsplit(seq1,"")))
	seq2_ch<-c("0","0")
	seq2_ch<-append(seq2_ch,unlist(strsplit(seq2,"")))
	D_back<-data.frame(matrix(NA, nrow = length(seq2_ch), ncol = length(seq1_ch)))
	D_back[1,]<-c(seq1_ch)
	D_back[,1]<-c(seq2_ch)
	D_back[2,2]<-"done"
	for(i in 2:length(seq1_ch)) {
		D_back[2,i]<-"left"
	}
	for(j in 2:length(seq2_ch)) {
		D_back[j,2]<-"up"
	}
	D_back[2,2]<-"done"
	return(D_back)
}


intit_matrix_back_end = function(seq1,seq2) 
{
	#Given two sequences of a.a., it returns initial traceback-matrix of Needleman-Wunsch for situations in which no penalty is put at the ends of sequences
	seq1_ch<-c("0","0")
	seq1_ch<-append(seq1_ch,unlist(strsplit(seq1,"")))
	seq2_ch<-c("0","0")
	seq2_ch<-append(seq2_ch,unlist(strsplit(seq2,"")))
	D_back<-data.frame(matrix(NA, nrow = length(seq2_ch), ncol = length(seq1_ch)))
	D_back[1,]<-c(seq1_ch)
	D_back[,1]<-c(seq2_ch)
	D_back[2,2]<-"done"
	for(i in 2:length(seq1_ch)) {
		D_back[2,i]<-"diag"
	}
	for(j in 2:length(seq2_ch)) {
		D_back[j,2]<-"diag"
	}
	D_back[2,2]<-"done"
	return(D_back)
}



fill_matrix_back = function(D_back,D,gap_penalty)
{
	#Given the initial traceback-matrix, the initial matrix and the gap_penalty, it fits the traceback-matrix for those sequences
	seq1_ch<-c("0","0")
	seq1_ch<-append(seq1_ch,unlist(strsplit(seq1,"")))
	seq2_ch<-c("0","0")
	seq2_ch<-append(seq2_ch,unlist(strsplit(seq2,"")))
	for(j in 3:length(seq2_ch)) {
		for(i in 3:length(seq1_ch)) {
		H=D[1,i]			
		V=D[j,1]
		S=BLOSUM62[colnames(BLOSUM62)==as.character(V),rownames(BLOSUM62)==as.character(H)]
		a=as.numeric(D[j-1,i-1])+S
		b=as.numeric(D[j-1,i])-gap_penalty
		c=as.numeric(D[j,i-1])-gap_penalty
		D[j,i]<-max(a,b,c)
		if(a==max(a,b,c))
		{
			D_back[j,i]<-"diag"
			next
		}
		if(b==max(a,b,c))
		{
			D_back[j,i]<-"up"
			next
		}
		if(c==max(a,b,c)) 
		{
			D_back[j,i]<-"left"
		}
		}
	}
	return(D_back)
}


provide_seq = function(DFB) {
	#Given the traceback-matrix of Needleman-Wunsch for two sequences of a.a., it returns the best possible alignment for those sequences
	seq1<-c()
	seq2<-c()
	move=DFB[nrow(DFB),ncol(DFB)]
	u=0
	v=0
	for(i in 1:9999) {
		if(move=="diag"){
			coord<-c(nrow(DFB)-u,ncol(DFB)-v)
			seq1<-append(seq1,DFB[1,coord[2]])
			seq2<-append(seq2,DFB[coord[1],1])
			move=DFB[nrow(DFB)-1-u,ncol(DFB)-1-v]
			u=u+1
			v=v+1
		}
		if(move=="up"){
			coord<-c(nrow(DFB)-u,ncol(DFB)-v)
			seq1<-append(seq1,"_")
			seq2<-append(seq2,DFB[coord[1],1])
			move=DFB[nrow(DFB)-1-u,ncol(DFB)-v]
			u=u+1
		}
		if(move=="left"){
			coord<-c(nrow(DFB)-u,ncol(DFB)-v)
			seq1<-append(seq1,DFB[1,coord[2]])
			seq2<-append(seq2,"_")
			move=DFB[nrow(DFB)-u,ncol(DFB)-1-v]
			v=v+1
		}
		if(move=="done"){
			break
		}
		i=i
	}
	lista<-list("alig1"=paste(rev(seq1),collapse=""),"alig2"=paste(rev(seq2),collapse=""))
	return(lista)
}




provide_seq_end = function(DFB,D_fill) {
	#Given the traceback-matrix of Needleman-Wunsch for two sequences of a.a., it returns the best possible alignment for those sequences, in situations in which end-sequence-gap penalty is different than the gap penalty in between the sequences
	seq1<-c()
	seq2<-c()
	D_fill_m<-D_fill[-1,-1]
	h<-sapply(D_fill_m, as.numeric)
	highest<-which(h == max(h), arr.ind = TRUE)+1
	move=DFB[nrow(DFB),ncol(DFB)]
	u=0
	v=0
	for(i in 1:9999) {
		if(move=="diag"){
			if(u==0 & v==0)
			{
				coord<-highest
			} else {
				coord<-c(nrow(DFB)-u,ncol(DFB)-v)
			}
			seq1<-append(seq1,DFB[1,coord[2]])
			seq2<-append(seq2,DFB[coord[1],1])
			move=DFB[nrow(DFB)-1-u,ncol(DFB)-1-v]
			u=u+1
			v=v+1
		}
		if(move=="up"){
			if(u==0 & v==0)
			{
				coord<-highest
				u=u+1
				v=v+1
			} else {
				coord<-c(nrow(DFB)-u,ncol(DFB)-v)
			} 
			seq1<-append(seq1,"_")
			seq2<-append(seq2,DFB[coord[1],1])
			move=DFB[nrow(DFB)-u,ncol(DFB)-v]
			u=u+1
		}
		if(move=="left"){
			if(u==0 & v==0)
			{
				coord<-highest
				move=DFB[nrow(DFB)-u,ncol(DFB)-v]
				v=v+1
				u=u+1
			} else {
				coord<-c(nrow(DFB)-u,ncol(DFB)-v)
				move=DFB[nrow(DFB)-u,ncol(DFB)-v]
			} 
			seq1<-append(seq1,DFB[1,coord[2]])
			seq2<-append(seq2,"_")
			v=v+1
		}
		if(move=="done"){
			break
		}
		i=i
	}
	lista<-list("alig1"=paste(rev(seq1),collapse=""),"alig2"=paste(rev(seq2),collapse=""))
	return(lista)
}




score = function(res1,res2){
	#given two a.a. residues, it returns the alignment score using the Blosum 62 score matrix
	S=BLOSUM62[colnames(BLOSUM62)==as.character(res1),rownames(BLOSUM62)==as.character(res2)]
	return(S)
	}


score_al = function(alig1,alig2,gap_penalty) {
	#given two aligned sequences and the gap_penalty, it returns the score of the alignment
	scores=c()
	alig1<-unlist(strsplit(alig1,""))
	alig2<-unlist(strsplit(alig2,""))
	for(i in 1:length(alig1)){
		if(alig1[i] == "_" | alig2[i] == "_")
		{
			scores<-append(scores,-gap_penalty)
		} else {
			scores<-append(scores,score(alig1[i],alig2[i]))
		}
	}
	return(sum(scores))
}


score_al_end_noP = function(alig1,alig2,gap_penalty,pen_gap_end) {
	#given two aligned sequences and the gap_penalty, it returns the score of the alignment, in situations in which end-sequence-gap penalty is different than the gap penalty in between the sequences
	scores=c()
	alig1<-unlist(strsplit(alig1,""))
	alig2<-unlist(strsplit(alig2,""))
	for(i in 1:length(alig1)){
		if(alig1[i] == "_" | alig2[i] == "_")
		{
			b<-alig2[i:length(alig2)]
			if(any(b != "_")){
				scores<-append(scores,-gap_penalty)
			} else {
				scores<-append(scores,-pen_gap_end)
			}
		} else {
			scores<-append(scores,score(alig1[i],alig2[i]))
		}
	}
	return(sum(scores))
}



per_identy = function(alig1,alig2)
{	
	#Given two sequences it returns the percentage of identity for thos sequences
	suma<-c()
	alig1<-unlist(strsplit(alig1,""))
	alig2<-unlist(strsplit(alig2,""))
	for(j in 1:length(alig2)) {
		for(i in 1:length(alig1)) {
			if(alig1[i] == alig2[j])
			{	
				suma<-append(suma,1)
			}
		}
	}
	return(sum(suma)/length(alig1)*100)
}
			





#Call your functions here
seq2<-"THISLINE"			#Naming the sequences like this, matrices will be like in book Ch5.2	
seq1<-"ISALIGNED"


#In any case
D<-intit_matrix(3,seq1,seq2)
D_fill<-fill_matrix(D,3)


#In situations where the end gap penalty == gap penalty in between the sequences
D_back<-intit_matrix_back(seq1,seq2)
DFB<-fill_matrix_back(D_back,D,3)
best_alig<-provide_seq(DFB)
score_al(best_alig$alig1,best_alig$alig2,3)
PI<-per_identy(best_alig$alig1,best_alig$alig2)

#In situations where the end gap penalty != gap penalty in between the sequences
D_back_end<-intit_matrix_back_end(seq1,seq2)
DFB_end<-fill_matrix_back(D_back_end,D,8)
best_alig_end<-provide_seq_end(DFB_end,D_fill)
score_al_end_noP(best_alig_end$alig1,best_alig_end$alig2,8,0)






