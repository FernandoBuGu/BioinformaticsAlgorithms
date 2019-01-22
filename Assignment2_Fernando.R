#!/usr/bin/env R

######################################
#Author: Fernando Bueno Gutierrez
#Student number: 890605143090
#Implementation of the SSAHA algorithm
######################################
#Input arguments:






# import statements
x<-c("stringr","data.table","reshape")
#lapply(x, install.packages)
lapply(x, require, character.only = TRUE)


# implement your functions here

kmers_overl<-function(listQ,k=2)
{
#Given a sequence it returns the overlaping kmers of size "k"
	#listQ: A character or a list of characters with letters corresponding to DNA/RNA bases or amino acids from the query sequence
	#k: kmers size
	#returns a character containing the kmers 

	kmers<-c()
	for(Q in 1:length(listQ)){
		seq<-listQ[Q]
		for(i in 1:nchar(seq)){
			kmer<-str_sub(seq, start=i, end=i+(k-1))
			if(nchar(kmer)==k) {kmers<-append(kmers,kmer)}
		}
	}
	return(kmers)
}



kmers_NONoverl<-function(seq,k=15)
{
#Takes a sequence and returns the non-overlaping kmers of size "k"
	#seq: A character with letters corresponding to DNA/RNA bases or amino acids
	#k: kmers size
	#returns a character containing the kmers 

	kmers<-c()
	for(i in seq(1,nchar(seq),by=k)) {
		kmer<-str_sub(seq, start=i, end=i+(k-1))
		if(nchar(kmer)==k) {kmers<-append(kmers,kmer)}
	}
	return(kmers)
}



return_L<-function(D,k=15)
{
#Takes the search sequences and returns list containing the positions of each kmer name
	#D: a list containing one or more search sequences
	#k: kmers size	
	#returns a list where each element is a kmer name and each subelement is the position at which it appears. 
		#output format: tupple: (search_sequence,base-pair)

	kmers_list<-lapply(D,kmers_NONoverl,k=k)
	kmer_name<-c()
	position<-c()
	for(seq in 1:length(kmers_list)){
		k_mers_seq<-kmers_list[seq]
		for(i in 1:length(unlist(k_mers_seq))) {
				kmer_name<-append(kmer_name,unlist(k_mers_seq)[i])
				position<-append(position,paste("(",seq,",",i*k-1,")",sep=""))
		}
	}	
	names(position)<-kmer_name
	position<-position[order(names(position))]

	List_kmers<-list()
	for(i in 1:length(position)) {
		kmer<-as.character(names(position))[i]
		List_kmers[[kmer]]<-position[names(position) == kmer]}

	return(List_kmers)
}



print_table<-function(D,k=15)
{
#Given a list of search sequences, it prints table indicating the positions of the different kmer names in the search sequences. 
	#D: a list containing one or more search sequences
	#k: kmers size
	#prints a matrix in the same format as table 1 in SSAHA paper.	

	List_kmers<-return_L(D,k)
	max.length <- max(sapply(List_kmers, length))
	List_kmers <- lapply(List_kmers, function(v) { c(v, rep(NA, max.length-length(v)))})
	List_kmers<-do.call(rbind, List_kmers)
	List_kmers[is.na(List_kmers)] <- ""
	colnames(List_kmers) <- rep("",dim(List_kmers)[2])
	print(noquote(List_kmers))
}



return_table<-function(D,k=15)
{
#As avobe but it returns object rather than printing it

	List_kmers<-return_L(D,k)
	max.length <- max(sapply(List_kmers, length))
	List_kmers <- lapply(List_kmers, function(v) { c(v, rep(NA, max.length-length(v)))})
	List_kmers<-do.call(rbind, List_kmers)
	colnames(List_kmers) <- rep("",dim(List_kmers)[2])
	return(List_kmers)
}


kmers_overl<-function(seq,k=2)
{
#Given a sequence it returns the overlaping kmers of size "k"
	#seq: A character with letters corresponding to DNA/RNA bases or amino acids
	#k: kmers size
	#returns a character containing the kmers 

	kmers<-c()
	for(i in 1:nchar(seq)){
		kmer<-str_sub(seq, start=i, end=i+(k-1))
		if(nchar(kmer)==k) {kmers<-append(kmers,kmer)}
	}
	return(kmers)
}


count_matches<-function(D,query,k)
{
#Given a list of search sequences and the query it returns a table indicating count of overlaping kmers of querry in search sequences
	#D: a list containing one or more search sequences
	#query: a character corresponding to the sequence to be aligned
	#returns a list containing a table indicating count of overlaping kmers, and a matrix that is the transpose of table in 
		#previous functions

	transpose_list<-t(return_table(D,2))
	kmers_search<-colnames(transpose_list)
	kmers_query<-unlist(kmers_overl(query))

	bp<-c() #base-pair index at which the kmer appears
	match<-c()	#kmer name of querry that was found in the search sequences
	for(i in 1:length(kmers_search)){
		if(kmers_search[i] %in% kmers_query){
			t<-which(kmers_query == kmers_search[i])
			bp<-append(bp,t-1)
			match<-append(match,kmers_query[t])
		}
	}
	names(bp)<-match
	bp<-sort(bp)
	count<-table(names(bp))
	output<-list("count"=count,"transpose_list"=transpose_list,"base_pairs"=bp)
	return(output)
}



return_list<-function(transpose_list,count,bp)
{
# Intermediate function: Given the output from previous function it returns input required for next
	#transpose_list: a matrix that is the transpose of table in previous functions
	#count: a table indicating count of overlaping kmers
	#returns a dataframe indicating the position at which each kmer was found
	
	another<-data.frame(matrix(nrow = dim(transpose_list)[1]))
	to<-dim(g$transpose_list)[2]
	for(i in 1:to){
		if(!colnames(transpose_list)[i] %in% names(count)){
			transpose_list[,i]<-NA
		} else {
			this<-which(names(count) == colnames(transpose_list)[i])
			repQ<-count[this]
			while(repQ != 1){
				length(another[,1]) = length((transpose_list)[,i])
				another<-cbind(another,as.factor((transpose_list)[,i]))
				j<-dim(another)[2]
				colnames(another)[j]<-colnames(transpose_list)[i]
				repQ<-repQ-1
			}
		}
	}
	another[,1]<-NULL
	another<-as.matrix(another)
	transpose_list<-cbind(transpose_list,another)
	transpose_list <- transpose_list[,colSums(is.na(transpose_list))<nrow(transpose_list)]
	list <- names(bp)
	transpose_list<-transpose_list[,list]
	list<-t(transpose_list)
	list<-as.data.frame(as.table(list))
	list[,2]<-NULL
	colnames(list)[1]<-"rn"
	return(list)
}


return_M<-function(list,bp)
{
# Intermediate function: Given the output from previous function it returns input required for next


	library(data.table)
	bp_df<-as.data.frame(t(data.frame(as.list(bp))))	
	bp_df<-setDT(bp_df, keep.rownames = TRUE)[]
	nam<-gsub("\\..*","",bp_df$rn)
	bp_df$rn<-nam

	library(reshape)
	DD<-melt.data.frame(list, "rn")
	DD$variable<-NULL
	DD<-DD[ order(DD[,1]), ]
	DD<-DD[complete.cases(DD[,2]),]	
	M<-merge(DD,bp_df,by="rn") #this is good
	M<-M[ order(M[,3]), ]

	library(stringr)
	a<-str_split_fixed(M$value, ",", 2)
	nam<-as.numeric(gsub(")","",a)[,2])
	M<-cbind(M,nam)
	M$dif<-M$nam-M$V1
	M$H<-paste(a[,1],",",M$dif,",",a[,2],sep="")
	M<-M[ order(M[,3]), ]
	return(M)
}



print_table2<-function(M)
{
#Given output from previous, it returns table2


	M<-unique(M)
	a<-str_split_fixed(M$H, ",", 3)
	a1<-a[,1]
	a1<-as.numeric(gsub("[^0-9]", "", a[,1]))
	a2<-a[,2]
	a2<-as.numeric(gsub("[^0-9]", "", a[,2]))

	M<-cbind(M,a1)
	M<-cbind(M,a2)
	table2<-M
	table2<-table2$H[order(table2$a1,table2$a2)]
	table2<-cbind(table2,M)
	table2$nam<-NULL
	table2$dif<-NULL
	table2$a1<-NULL
	table2$a2<-NULL
	table2<-table2[,c(4,2,3,5,1)]
	colnames(table2)<-c("t","W","position","H","M")

	for(i in 2:dim(table2)[1]){
		if(table2[i,2]==table2[i-1,2]){
			table2[i,1]<-""}
	}
	table2$W<-as.character(table2$W)
	for(i in 2:dim(table2)[1]){
		if(table2[i,1]==""){
			table2[i,2]<-""}
	}
	return(table2)
}


report_max_match<-function(M,k=2)
{
#Given M, it reports the maximum match between a query sequence and the search sequences

	a<-str_split_fixed(unique(M$H), ",", 3)
	a1<-as.numeric(gsub("[^0-9]", "", a[,1]))
	a2<-as.numeric(a[,2])
	aj<-paste(a1,a2,sep=",")
	aj<-table(aj)
	maxi<-names(aj)[aj==max(aj)]
	maxi2<-as.numeric(str_split_fixed(maxi, ",", 2))
	length_al<-max(aj)*k
	printy<-paste("the maximum match is with search ",maxi2[1]," at position ",maxi2[2], " of search ", maxi2[1], " and continuing for ", length_al, " bases")
	l<-list("length_al"=length_al,"search"=maxi2[1],"from"=maxi2[2],"print"=printy)
	return(l)
}


print_alignment<-function(alig_seqs,D,query)
{
#Prints the best alignment between a seqrch sequence and the query

	S=unlist(D[alig_seqs$search])
	query=unlist(query)
	before<-paste(rep(" ",alig_seqs$from-2),collapse="")
	query<-paste(before,query,collapse="")
	later<-nchar(S)-nchar(query)
	later<-paste(rep(" ",later-1),collapse="")
	query<-paste(query,later,collapse="")
	interm<-paste(rep(" ",nchar(S)),collapse="")

	q<-str_split_fixed(query, "",nchar(query))
	S<-str_split_fixed(S, "",nchar(S))
	intT<-str_split_fixed(interm, "",nchar(interm))

	m<-matrix(ncol=length(S))
	m<-rbind(m,q)
	m<-rbind(m,intT)
	m<-rbind(m,S)
	m<-m[2:4,]

	for(i in 1:44){
		if(m[1,i]==m[3,i]){
			m[2,i]<-as.character("|")}
		if(m[1,i]!=m[3,i] & m[1,i]!=" "){
			m[2,i]<-as.character(".")}
	}

	que<-noquote(paste0(m[1,],collapse=""))
	int<-noquote(paste0(m[2,],collapse=""))
	Se<-noquote(paste0(m[3,],collapse=""))

	out<-capture.output(que,int,Se)
	print(noquote(out),row.names = FALSE)
}


answer_fasta<-function(link)
{
	thepage = readLines(as.character(link))
	seqs<-c()
	for(i in 1:length(thepage)){
		a<-substring(thepage[i], 1, 1)
		if(a==">"){
		seqs<-append(seqs,i)
		}
	}
	seqs<-append(seqs,length(thepage))
	length<-nchar(thepage[2])
	L<-c()
	for(i in 2:length(seqs)){
		L<-append(L,seqs[i]-seqs[i-1]-1)
	}
	L<-append(L,length(thepage))
	lengths<-L*length
	total_length<-round(sum(lengths)/1000000,2)
	printy<-paste("# seqs is",length(seqs)-1,".","Total length is", total_length, "Mbp")
	l<-list("print"=printy,"positions"=seqs,"page"=thepage)
}


parse_fasta<-function(positions,thepage)
{
	headers<-c()
	sequences<-list()
	le<-as.numeric(length(positions)-1)
		for(i in 1:le){
			header<-positions[i]
			headers<-append(headers,thepage[header])
			start<-positions[i]+1
			end<-positions[i+1]-1
			seq<-thepage[start:end]
			bases<-toString(seq)
			bases<-gsub(", ", "", bases)
			bases<-gsub("[^A,C,T,G]", "", bases)
			sequence<-c()
			sequence<-append(sequence,bases)
			sequences<-append(sequences,list(sequence))
		}
	lista<-list("headers"=headers,"sequs"=sequences)
	return(lista)
}



count_kmers_NONoverl<-function(list_seqs,k=15)
{
#Takes a sequence and returns the non-overlaping kmers of size "k"
	#seq: A character with letters corresponding to DNA/RNA bases or amino acids
	#k: kmers size
	#returns a character containing the kmers 

	kmers_count<-0
	for(S in 1:length(list_seqs)){
		seq=list_seqs[S]
		for(i in seq(1,nchar(seq),by=k)) {
			kmer<-str_sub(seq, start=i, end=i+(k-1))
			if(nchar(kmer)==k & !kmer %in% kmers_count) {kmers_count<-kmers_count+1}
		}
	}
	return(kmers_count)
}









# To see answers, please execute the following commands in order   
    query = 'TGCAACAT'
    s1 = 'GTGACGTCACTCTGAGGATCCCCTGGGTGTGG'
    s2 = 'GTCAACTGCAACATGAGGAACATCGACAGGCCCAAGGTCTTCCT'
    s3 = 'GGATCCCCTGTCCTCTCTGTCACATA'
    D<- list(s1,s2,s3)

#question 1
KMERS<-kmers_overl(query,2)
print_table(D,2)


#question 2
g<-count_matches(D,query,2)
p<-return_list(g$transpose_list,g$count,g$base_pairs)
M<-return_M(p,g$base_pairs)
table2<-print_table2(M)
print(table2,row.names = FALSE)


#question 3
out<-report_max_match(M,2)
print(out$print,row.names = FALSE)
print_alignment(out,D,query)


#question 4
link="http://www.bioinformatics.nl/courses/BIF-31306/TAIR10.fasta"
h<-answer_fasta(as.character(link))
print(h$print)
parsed<-parse_fasta(h$position,h$page)		


#question 5
kmers_count_list<-count_kmers_NONoverl(parsed$sequs,k=15)
paste<-paste("the total # of keys is",kmers_count_list)
print(paste)


#question 6
link="http://www.bioinformatics.nl/courses/BIF-31306/athal_query.fasta"
queries<-answer_fasta(as.character(link))
queries_parsed<-parse_fasta(queries$position,queries$page)
KMERS<-kmers_overl(queries_parsed$sequs,15)

	#To reduce computational time, modify function to consider only kmers taht are in query
	kmers_NONoverl<-function(seq,k=15)
	{
	#Takes a sequence and returns the non-overlaping kmers of size "k"
		#seq: A character with letters corresponding to DNA/RNA bases or amino acids
		#k: kmers size
		#returns a character containing the kmers 

		kmers<-c()
		for(i in seq(1,nchar(seq),by=k)) {
			kmer<-str_sub(seq, start=i, end=i+(k-1))
			if(nchar(kmer)==k & kmer %in% KMERS) {kmers<-append(kmers,kmer)}
		}
		return(kmers)
	}

#Something I do wrong because it takes more than 1 hour. I did not manage to see the results. Commands would be:

#		for(queri in 1:length(queries_parsed$sequs)){
#			query=queries_parsed$sequs[queri]
#			g<-count_matches(parsed$sequs,query,2)
#			p<-return_list(g$transpose_list,g$count,g$base_pairs)
#			M<-return_M(p,g$base_pairs)
#			table2<-print_table2(M)
#			print(table2,row.names = FALSE)
#			out<-report_max_match(M,2)
#			P<-paste("for query", queri, "result is:")
#			print(P)
#			print(out$print,row.names = FALSE)
#		}

