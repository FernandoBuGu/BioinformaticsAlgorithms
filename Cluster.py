#!/usr/bin/env python
"""
Author: Fernando Bueno Gutierrez
Implementation of the Lloyd’s k-means clustering algorithm to identify proteins that have similar values in the different conditions


# import statements
import random   
import math
import numpy
from itertools import combinations



def csv_parser(lines):

    """Returns list of point coordinates as 
	[[x1,x2,...,xd],[y1,y2,...,yd],...,[z1,z2,...,zd]]
	where "d" is the number of dimensions

    Arguments
    ---------    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Returns
    ---------
    data_points: List of lists with the coordinates of the points.
	The first deparation are the points, then the coordinates.

    Notes
    ---------
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 

    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(map(float, items[1:])) #first item is the label
        except ValueError: #must be the header
            continue
    return data_points



def euc_distance(point1,point2):
	"""
    Returns an integer corresponding to the euclidean distance bt two points

    Arguments
    ---------
	point1: list of int, each int is a coordinate of a point in a d-dimensional 
            space, as follows: [x1,x2,...,xd]. 
	point2: list of int, each int is a coordinate of another point in a 
            d-dimensional space, as follows: [y1,y2,...,yd]. 

    Returns
    ---------
	euc_distance: int, euclidean distance bt the two points	 
	"""

	summing=[]  #summing: list of square differences
	for idx, dim in enumerate(point1):
		summing.append((point1[idx]-point2[idx])**2)
	euc_distance=math.sqrt(sum(summing))
	return euc_distance



def chosse_cent(data_points,k):
	"""
    Returns a list of integers corresponding to the indices of the points 
	that will be the centroids in the first iteration of k-means. 
	These indices refer to "data_points".

    Arguments
    ---------
	data_points: List of lists with the coordinates of the points.
	             The first separation are the points, then the coordinates.

	k:           Int, number of clusters that we want to create. More 
	             specifically, here "k" refers to the number of centroids that 
	             we want to consider in the first iteration

    Returns
    ---------
	centroids_1st: list of int, "data_points" indices of the points that will 
                   be the centroids in the first iteration of k-means. 
	"""
	centroids_1st=[]
	#random.seed(3)
	while len(centroids_1st)<k:
		index = random.randint(0,len(data_points)-1)
        #To avoid taking any point more than once
		if index not in centroids_1st: 
			centroids_1st.append(index) 
	return centroids_1st



def fun_points_iterate(centroids_1st,data_points):
	"""
    Returns a list of integers corresponding to the indices of the points 
	over which the algorithm will iterate in the first iteration of k-means. 
	These indices refer to "data_points". The points will be those that were 
	not chosen as centroids in the function "chosse_cent". 

    Arguments
    ---------
	centroids_1st: list of int, "data_points" indices of the points that will 
                   be the centroids in the first iteration of k-means. 

	data_points:   list of int, coordinates of the points.
	               The first separation are the points, then the coordinates.

    Returns
    ---------
	points_iterate: list of int, "data_points" indices of the points over which 
                    we will iterate in the first iteration of k-means. 
	"""

	points_iterate=[]
	for idx, point in enumerate(data_points):
		if idx not in centroids_1st:
			points_iterate.append(idx)
	return points_iterate


def assign_to_closest(centroids,data_points,points_iterate,first=False): 
	"""
    Returns a list of integers of length equal to the #points we iterate 
	over. The integers correspond to the cluster to which each point has been 
	assigned.

    Arguments
    ---------
	centroids: If first is True: 
                    list of int, "data_points" indices of the points that will 
	                be the centroids in the first iteration of k-means.
	                This objects is called "centroids_1st" outside this 
                    function.

               If first is False: 
	                list of lists, coordinates of the new centroids.
	                The outer lists are the centromers. The inner list is a 
	                list of integers of length d, containing the coordinates of 
	                the centroid in a d-dimensional space. This objects is  
	                called "new centroids" outside the function.

	data_points:    list of lists with the coordinates of the points.
	                The inner lists are points, integers are coordinates.

	points_iterate: list of int, "data_points" indices of the points over which 
	                the algorithm will iterate in the iterations of k-means. 

	first:          logical, True if we are in the first iteration of k-means.


    Returns
    ---------
	assign_step1: If first is True: 
	                    list of int of length equal to the #points the 
	                    algorithm  iterates over. The integers correspond to  
	                    the cluster in which each point has been assigned 
                        (i.e. "[1, 0, 1]" means that the first and last point 
	                    were assigned to one cluster and the second point was 
	                    assigned to a differenet cluster).

                  If first is False: 
	                    list of int of length equal to the #points the 
	                    algorithm  iterates over. The integers correspond to  
	                    the cluster in which each point has been assigned (i.e. 
	                    "[1, 0, 1]" means that the first and last point  were  
	                    assigned to one cluster and the second point was 
                        assigned to a differenet cluster). 
	"""



	clusters = []
	for i in range(len(centroids)):
		clusters.append(i)

	C=0
	assign_step1=[]	
	centr_used=[] #centroids already used			
	for point in points_iterate:
		point=data_points[point]
		Distance_p_c=[] #Distance_p_c: distances from point to centroidS
		for cent in centroids:
			if first:
				centr=data_points[cent]
			else:
				centr=cent #cent=centr=centroids
			Distance_p_c.append(euc_distance(point,centr))
			if len(Distance_p_c)==len(centroids): 
                # If the distance to all entroids has been appended
				centr=[i for i in range(len(Distance_p_c)) if \
                    Distance_p_c[i] == min(Distance_p_c)]
				centr="".join(map(str, centr))
				centr=int(centr)
                # Answer to question 7: unique labelling.
				if centr not in centr_used: 
					centr=clusters[C]
					centr_used=centr_used+list([clusters[C]])
					C+=1
				assign_step1.append(centr)




	#if a cluster is empty, allocate a random point in the empty cluster
	#If this part was ignored, then we would end up with a lower number of 
	#clusters. If the algorithm is run untill convergence, this random choice 
	#will not matter (see Question 7)


	for idx, cluster in enumerate(clusters):
		if cluster not in assign_step1:
			empty_cluster=cluster
			assign_to_empty=random.randint(0,len(assign_step1)-1)
			assign_step1[assign_to_empty]=empty_cluster



	return assign_step1



def assign_to_closest_step2(assign_1st_step1,data_points,\
	points_iterate):
	"""
    Returns a list of lists of length equal to the #clusters defined. 
	The inner lists correspond to the clusters and the integers correspond
	to the "data_points" indices of the points that correspond to that cluster 

    Arguments
    ---------
	assign_1st_step1: list of int of length equal to the #points we 
	                  iterate over. The integers correspond to the cluster in 
	                  which each point  has been assigned (i.e. "[1, 0, 1]" 
	                  means that the first and last point were assigned to one 
	                  cluster and the second point was assigned to a  
                      differenet cluster).

	data_points:      List of lists, inner lists are points, integer are  
	                  coordinates.

	points_iterate:   list of int, data_points" indices of the points over 
                      which we will iterate in the first iteration of k-means.

    Returns
    ---------	 
	assign_1st_step2: list of lists of length equal to the #clusters defined. 
	                  The inner lists correspond to the clusters and the 
	                  integers correspond to the "data_points" indices of the 
	                  points that correspond to that cluster 
	                   (i.e "[[0, 2], [3]]" means that the 1st and 3rd number  
	                  ("[0]" and "[2]")  of "data_points" correspond to the 
	                  same cluster, whereas the 4st number ("[3]") corresponds 
                      to a different cluster. The 2nd number ("[1]") was        
                      defined as a centroid).
	"""

	used_clusters=[]# clusters for which at least one point has been added
	for idx, code in enumerate(assign_1st_step1):
		if assign_1st_step1[idx] not in used_clusters:
			used_clusters.append(assign_1st_step1[idx])

	list.sort(used_clusters)	#Question 7, unique labelling

	assign_1st_step2=[] #Assign in first iteration, second step
	for idx, clust in enumerate(used_clusters):
		list_within=[i for i in range(len(assign_1st_step1)) if \
		assign_1st_step1[i] == used_clusters[idx]]

		
		lis_clust_small=[] 
        #lis_clust_small: Indexes of the points within a cluster
		for ele in list_within: #list_within: list of points within a cluster
			list_within_index=points_iterate[ele]
			lis_clust_small.append(list_within_index)
			if len(lis_clust_small) == len(list_within):
				assign_1st_step2.append(lis_clust_small)


	return assign_1st_step2


def fun_points_iterate_after1st(data_points):
	"""
    Returns a list of integers corresponding to the indices of the points 
	over which we will iterate in the second and next iterations of k-means. 
	These indices refer to "data_points". The points will be all in 
	"data_points". 

    Arguments
    ---------
	data_points: List of lists with the coordinates of the points.
	             The first deparation are the points, then the coordinates.

    Returns
    ---------	
	points_iterate: List of int, "data_points" indices of the points over which 
                    the algorithm will iterate in the first iteration of 
	                k-means. 
	"""
	points_iterate=[]
	for idx, point in enumerate(data_points):
		points_iterate.append(idx)
	return points_iterate


def new_centroids2(data_points,assign_1st_step2):
    """
    Given the coordinates of the points and their assigination to clusters,
    returns the coordinates of the new centroids.

    Arguments
    ---------
    data_points:      List of lists with the coordinates of the points.
                      The inner lists are points, the integers are coordinates.

    assign_1st_step2: List of lists of length equal to the #clusters defined. 
                      The inner lists are clusters and the integers are the 
                      "data_points" indices of the points that correspond to 
                      that cluster

    Returns
    ---------	
    new_centroids: list of lists of length equal to the number of clusters.
                   Inside the inner list, there are "d" integers corresponding  
                   to the coorinates of the new centroid for that cluster in a 
                   d-dimensional space. The coordiantes were computed by taking 
                   the mean of the coordiantes of all points in the cluster.
    """
    clusters=[]
    for clusts in assign_1st_step2:
	    L=[]
	    clusters.append(L)



    for idx, cluster in enumerate(assign_1st_step2):
	    for point in cluster:
		    clusters[idx].append(data_points[point])

    leng=len(data_points[0])
    dims=[]
    for Dim in range(leng):
	    L=[]
	    dims.append(L)


    Coor=[] #Coor: list of lists of lists with the coordinates of the points of 
            #each cluster
    for cluster in clusters:
	    dims=[]
	    for Dim in range(leng):
		    L=[]
		    dims.append(L)
	    for points in cluster:
		    for ide, dim in enumerate(points):
			    dims[ide].append(dim)
	    Coor.append(dims)

    new_centroids=[]
    for cluster in Coor:
	    centroid=[]
	    for dim in cluster:
       #av: average coordinate among points in a cluster, in a  given dimension
		    av=round(numpy.mean(dim),4) 
		    centroid.append(av)
	    new_centroids.append(centroid)

    return new_centroids


def final_location(points_iterate_after1st,new_centroids,data_points,\
	tolerance=0):
	"""	
	Returns a list of lists, where the outer lists refer to the clusters and 
	the inner lists contain the "data_points"-indices of the points that were 
	assigned to the cluster. It also returnsa an integer corresponding to the 
	number of iterations. 

    Arguments
    ---------
	points_iterate_after1st: list of int "data_points" indices of the points 
                             over which we iterate in the first iteration of 
                             k-means. 

	new_centroids:           list of lists containing the coordinates of the 
	                         new centroids. The outer lists are the centromers. 
                             The inner list is a list of integers of length d, 
                             containing the coordinates of the centroid in a 
	                         d-dimensional space. 

	data_points:             List of lists with the coordinates of the points.
	                         The first deparation are the points, then the 
                             coordinates.

	tolerance:               int, minimum euclidean distance between the
	                         centroids of the previous iteration and the new 
                             centroids that is required before k-means 
                             converges. For k centroids, "tolerance" refers to 
                             the mean of the euclidean distances between each 
                             centroid in the current interation, and the 
                             corresponding centroid in the previous iteration.

    Returns
    ---------	
	cluster_solultion:  list of lists, where the outer lists refer to the 
	                    clusters and the inner lists contain the "data_points"
                        indices of the points that were assigned to the 
                        cluster.

	noit:               int, number of iterations.

	new_centroids:      list of lists containing the coordinates of the new
                        centroids. The outer lists are the centromers. 
                        The inner list is a list of integers of length d, 
                        containing the coordinates of the centroid in a 
                        d-dimensional space.

	data_points:        List of lists with the coordinates of the points.
	                    The first deparation are the points, then the 
                        coordinates.
	"""

	it=1 #number of iterations
	convergence=False	#Whether the algorithm has converged
	if tolerance == 0:
		while not convergence:
			assign_next_step1=assign_to_closest(new_centroids,\
				data_points,points_iterate_after1st)
			cluster_solultion=assign_to_closest_step2(assign_next_step1,\
				data_points,points_iterate_after1st)
			previous_centroids=new_centroids
			new_centroids=new_centroids2(data_points,cluster_solultion)
			it+=1
			if previous_centroids == new_centroids:
				convergence = True

		noit= it
		return cluster_solultion, noit, new_centroids, data_points

	euc_dista_all=99999 #symbolic value to make sure that it starts
	if tolerance != 0:
		while euc_dista_all > tolerance: 
            #Then the value of euc_dista_all is reduced and will be of same 
            #magnitude as "tolerance".  
			euc_dista_all_list=[]
			assign_next_step1=assign_to_closest(new_centroids,\
				data_points,points_iterate_after1st)
			cluster_solultion=assign_to_closest_step2(assign_next_step1,\
				data_points,points_iterate_after1st)
			previous_centroids=new_centroids
			new_centroids=new_centroids2(data_points,cluster_solultion)

			for idx, new_centroid in enumerate(new_centroids):
				euc_dista=abs(euc_distance(previous_centroids[idx],\
					new_centroids[idx]))
				euc_dista_all_list.append(euc_dista)
			euc_dista_all=numpy.mean(euc_dista_all_list)
			it+=1
			if previous_centroids == new_centroids:
				convergence = True

		noit= it
		return cluster_solultion, noit, new_centroids, data_points



def run_Lloyd_algorithm(csv_file,k=3,tolerance=0):
	"""
    Calls all required functions to run the LLoyd's algorithm. It returns a 
	list of lists, where the outer lists refer to the clusters and the inner 
	lists contain the "data_points"-indices of the points that were 
	assigned to the cluster. It also returnsa an integer corresponding to the 
	number of iterations. 

    Arguments
    ---------	
	csv_file: path of the input file in .csv format.
	
	k:        int, number of clusters in which we want to separate the points. 

    Returns
    ---------
	solution: objects "cluster_solultion" and "noit".

              cluster_solultion:  list of lists, where the outer lists refer to  
	                              the clusters and the inner lists contain the 
                                  "data_points" indices of the points that were 
                                  assigned to the cluster.

	                       noit:  int, number of iterations.
	"""

	data_points = csv_parser(open(csv_file))
	centroids_1st=chosse_cent(data_points,k)
	points_iterate=fun_points_iterate(centroids_1st,data_points)
	assign_1st_step1=assign_to_closest(centroids_1st,data_points,\
		points_iterate,first=True)
	assign_1st_step2=assign_to_closest_step2(assign_1st_step1,data_points,\
		points_iterate)
	new_centroids=new_centroids2(data_points,assign_1st_step2)
	points_iterate_after1st=fun_points_iterate_after1st(data_points)
	solution=final_location(points_iterate_after1st,new_centroids,\
		data_points,tolerance)


	return solution


def choose_k(data_points, cluster_solultion, new_centroids):
	"""
	Returns a string specifying what are the maximum distances within and 
	between clusters. By running several times for different k values, this
	information allows to choose a proper value of k. 

    Arguments
    ---------
	data_points:        list of lists with the coordinates of the points.
	                    The first deparation are the points, then the 
                        coordinates.

	cluster_solultion:  list of lists. The outer lists refer to the 
	                    clusters and the inner lists contain the "data_points"
                        indices of the points that were assigned to the 
                        cluster.

	new_centroids:      list of lists containing the coordinates of the 
	                    new centroids. The outer lists are the centromers.
                        The inner list is a list of integers of length d, 
                        containing the coordinates of the centroid in a 
	                    d-dimensional space. 

    Returns
    ---------
	within:   str, maximum istance within clusters. This distance corresponds 
              to  the maximum euclidean distance between a point and its 
              corresponding centroid. 

	between:  str, maximum istance between clusters. Corresponds to the maximum 
	          euclidean distance between the centroids of two clusters. 
	"""

	euc_distances_clusters=[]
	for idx, cluster in enumerate(cluster_solultion):
		euc_distances=[]
		for i, point in enumerate(cluster):
			ed=euc_distance(data_points[cluster[i]],new_centroids[idx])
				#ed: euclidean distance
			euc_distances.append(ed)
		euc_distances_clusters.append(euc_distances)

	within=[]
	for cluster in euc_distances_clusters:
		maxi=max(cluster)
		within.append(maxi)
	within_max=round(max(within),4)
	within="The maximum distances within clusters is {}".format(within_max) 

	between=[]
	for combo in combinations(new_centroids, 2):
		ed=euc_distance(combo[0],combo[1])
		between.append(ed)
	between=round(max(between),4)
	between="The maximum distances between clusters is {}".format(between) 

	return within, between



if __name__ == "__main__":


    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.


    #QUESTION 1
    print "\nQUESTION 1 (a) \n"
    for i in range(0,10):
        FINAL1=run_Lloyd_algorithm\
        ('2dtest.csv',3)
        print "Run {}".format(i), FINAL1[0]


    print "\nQUESTION 1 (b) \n"
    nit=[]
    for i in range(0,1000):
        FINAL1=run_Lloyd_algorithm\
        ('2dtest.csv',3)
        nit.append(FINAL1[1])
    print "over 1000 runs, the mean of no.iterations needed to converge was:"
    print round(numpy.mean(nit),2)
    print
    print "and the standrad deviation was:"
    print round(numpy.std(nit),2)




    #QUESTION 2
    print "\nQUESTION 2 (a) \n"
    for i in range(1,7):
        FINAL1=run_Lloyd_algorithm\
        ('LargeSet_1.csv',i)
        print "For k={}".format(i), FINAL1[0]


    print "\nQUESTION 2 (b) \n"
    nit=[]
    for i in range(1,7):
        FINAL1=run_Lloyd_algorithm\
        ('LargeSet_1.csv',i)
        nit.append(FINAL1[1])
    print "over 6 k values (1:6), the mean no.iterations was:"
    print round(numpy.mean(nit),2)
    print
    print "and the standrad deviation was:"
    print round(numpy.std(nit),2)


    #QUESTION 3
    Set2_k2=run_Lloyd_algorithm\
    ('LargeSet_2.csv',2)
    print "\nQUESTION 3 \n"
    print "The algorithm tries to recover the correct cluster structure:"
    print
    print Set2_k2[0]
    print
    print "But fails for a few points"


    #QUESTION 7
    print "\nQUESTION 7 \n"
    data_points = csv_parser(open(\
	    '2dtest.csv'))
    centroids_1st=chosse_cent(data_points,3)
    points_iterate=fun_points_iterate(centroids_1st,data_points)
    assign_1st_step1=assign_to_closest(centroids_1st,data_points,\
	    points_iterate,first=True)
    print "The labels have been permuted. For example:"
    print "Run1", assign_1st_step1
    print "Run2", assign_1st_step1
    centroids_1st=chosse_cent(data_points,3)
    points_iterate=fun_points_iterate(centroids_1st,data_points)
    assign_1st_step1=assign_to_closest(centroids_1st,data_points,\
	    points_iterate,first=True)
    print "Run3", assign_1st_step1
    centroids_1st=chosse_cent(data_points,3)
    points_iterate=fun_points_iterate(centroids_1st,data_points)
    assign_1st_step1=assign_to_closest(centroids_1st,data_points,\
	    points_iterate,first=True)
    print "Run4", assign_1st_step1
    print
    print "The cluster of the firts point will always be labeled as '0', the" 
    print "cluster of the first point after this whose membership is different"  
    print "will always be labeled as '1', and so on. Thus this situation:"
    print "[0, 0, 0, 0, 0, 1, 1, 1, 1, 1] or [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]"
    print "will not occur, because both lists will always be labeled as the"
    print "first list" 



    #QUESTION 8
    FINAL1=run_Lloyd_algorithm\
    ('2dtest.csv',3)

    print "\nQUESTION 8 \n"
    print FINAL1[0]
    print(choose_k(FINAL1[3], FINAL1[0],FINAL1[2]))


    #QUESTION 9
    print "\nQUESTION 9 \n"
    FINAL1=run_Lloyd_algorithm\
	    ('2dtest.csv',5,0.01)
    print "tolerance:0.01", FINAL1[0]
    print "tolerance:0.01", FINAL1[1], "\n"
    FINAL1=run_Lloyd_algorithm\
	    ('2dtest.csv',5,0.2)
    print "tolerance:0.2", FINAL1[0]
    print "tolerance:0.2", FINAL1[1], "\n"




