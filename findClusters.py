import numpy as np
import matplotlib.pyplot as plt
import bisect
import math
import jp3fctn as jp
import itertools as itt

r = 4 # num rows
c = 4 # num col
N = r*c # number of sites
orderexp = 8 # order of expansion

# find all clusters
clusters, numclusters = jp.genAllClus(N,orderexp) # all clusters and their orders

# find all liked clusters
linkedclusters, numlinclus = jp.onlyLinkedClus(clusters, numclusters, r, c, orderexp)

# just an example
a = linkedclusters[40]

# find all symmetrically distinct clusters
printoutput = False
if printoutput:
    print("Symmetrically Distinct")
#symdistinct, numsymclus = jp.onlySymDist(linkedclusters, numlinclus, r, c, printoutput)

# finds topologically distinct clusters
# this is already saved, so we don't need to re-run this

zerc = np.zeros(c)
zerr = np.zeros(r)
Lvalues = []

# for each order
for order in symdistinct:
    orderclus = {}
    # for each symmetrically distinct cluster of that order
    for symclus,val in order.items():
        # find all clusters cg[i] related by a symmetry
        cg = jp.genclusters(linkedclusters[symclus])
        distinct = []
        
        # check if each cluster is distinct under translation
        for clus in cg:
            newclus = True
            
            # does clus equal any old array
            for comp in distinct:
                if np.array_equal(comp, clus): 
                    newclus = False
                    break
                    
            # try to find translation vector
            # translate to the right
            # column is blank
            currcol = c - 1
            while np.array_equal(clus[:,currcol], zerr) and newclus:
                # temp = translated clus
                temp = np.roll(clus,c-currcol,axis=1)

                # check whether translated clus is in distinct
                for comp in distinct:
                    if np.array_equal(comp, temp): 
                        newclus = False
                        break

                currcol = currcol - 1

            # translate to the left
            # column is blank
            currcol = 0
            while np.array_equal(clus[:,currcol], zerr) and newclus:
                # temp = translated clus
                temp = np.roll(clus,-currcol-1,axis=1)

                # check whether translated clus is in distinct
                for comp in distinct:
                    if np.array_equal(comp, temp): 
                        newclus = False
                        break

                currcol = currcol + 1

            # translate up and left/right
            # row is blank
            currow = 0
            # if next row only 0, then translate up
            while np.array_equal(clus[currow,:], zerc) and newclus:

                # temp = translated clus
                temp = np.roll(clus,-currow-1,axis=0)
                    #for k in range(currcol+1):
                        #temp[r-1-k,:] = 0

                # check whether translated clus is in distinct
                for comp in distinct:
                    if np.array_equal(comp, temp): 
                        newclus = False
                        break

                # try to find translation vector
                # translate to the right
                # column is blank
                currcol = c - 1
                while np.array_equal(temp[:,currcol], zerr) and newclus:
                    # temp = translated clus
                    temp2 = np.roll(temp,c-currcol,axis=1)


                    # check whether clus is in distinct
                    for comp in distinct:
                        if np.array_equal(comp, temp2): 
                            newclus = False
                            break

                    currcol = currcol - 1

                # translate to the left
                # column is blank
                currcol = 0
                while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                    temp2 = np.roll(temp,-currcol-1,axis=1)
                        #for k in range(currcol+1):
                        #    temp2[:,c-1-k] = 0

                    # check whether translated clus is in distinct
                    for comp in distinct:
                        if np.array_equal(comp, temp2): 
                            newclus = False
                            break
                        
                    currcol = currcol + 1

                currow = currow + 1

            # translate down and left/right
            # row is blank
            currow = r - 1
            while np.array_equal(clus[currow,:], zerc) and newclus:
                    # temp = translated clus
                temp = np.roll(clus,r-currow,axis=0)
                    #temp[range(r-currow),:] = 0

                # check whether clus is in distinct
                for comp in distinct:
                    if np.array_equal(comp, temp): 
                        newclus = False
                        break

                # try to find translation vector
                # translate to the right
                # column is blank
                currcol = c - 1
                while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                    temp2 = np.roll(temp,c-currcol,axis=1)

                    # check whether clus is in distinct
                    for comp in distinct:
                        if np.array_equal(comp, temp2): 
                            newclus = False
                            break
                        
                    currcol = currcol - 1

                # translate to the left
                # column is blank
                currcol = 0
                while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                    temp2 = np.roll(temp,-currcol-1,axis=1)
                    #    for k in range(currcol+1):
                     #       temp2[:,c-1-k] = 0

                    # check whether translated clus is in distinct
                    for comp in distinct:
                        if np.array_equal(comp, temp2): 
                            newclus = False
                            break
                    
                    currcol = currcol + 1            

                currow = currow - 1

            # if not (i.e. if it is new)
            if newclus:
                distinct.append(clus)
        orderclus[symclus] = len(distinct)
    Lvalues.append(orderclus)
    total = 0
    for a,b in orderclus.items():
        total = total + b
    print(total)
print(Lvalues)


# find all topologically distinct clusters
if printoutput: 
    print()
    print("Topologically Distinct")
topodistinct = jp.onlyTopologicallyDistinct(Lvalues, linkedclusters, r, c, printoutput)

#### this finds Y, the multiplicity matrix ####
# this has also been saved

clusarray = topodistinct
totalnumclus = 0
toPrint = False
ydict={}
for toporder in clusarray:
    for comp, val in toporder.items():
        ydict[comp]=totalnumclus
        totalnumclus = totalnumclus +1
Y = np.zeros((totalnumclus, totalnumclus))

# generate all possible connected subclusters
# for each cluster
subclusters = []
numSubClusters = []

for clusorder in range(orderexp):
    for clus, val in clusarray[clusorder].items():
        clus = int(clus)

    #  find indixes of non-zero points in cluster 
        # first find any nonzero entry
        found = False
        rnonzero = 0
        cnonzero = 0

        for i in range(r):
            if found:
                break;
            for j in range(c):
                if linkedclusters[clus][i][j] == 1:
                    rnonzero = i
                    cnonzero = j
                    found = True
                    break;
                    
        # find array of indices of nonzero 
        indnonzero = []
        vis = np.zeros((r,c))
        indnonzero = jp.siteIndices(rnonzero, cnonzero, r, c, linkedclusters[clus], indnonzero, vis)

        i = 1
        #for i in range(number of points in cluster)
        while (i < (clusorder + 2)):
            if toPrint: print(i)
            # choose i of those points to be one and store the arrays
            combs = itt.combinations(indnonzero, i)
            for points in combs:
                arr = np.zeros((r,c))
                for pair in points:
                    arr[pair[0]][pair[1]] = 1
                # check if connected
                vis = np.zeros((r,c))
                arr2 = np.zeros((r,c))
                jp.maxCon2d(pair[0], pair[1], r, c, arr, arr2, vis)
                if np.array_equal(arr, arr2):
                    # then arr is connected
                    # note: we have found a subcluster s with i nonzero elements
                    # s is a subcluster of clus
                    if toPrint: print(arr)
                    if toPrint: print()
                    newclus = True 
                    adjclus = jp.adjmatrix(arr, i-1, r, c) # check i - 1

                    # for each matrix in list of topologically distinct matrices of current order
                    for comp, val in clusarray[i-1].items():
                        if not newclus:
                            break
                        orig = linkedclusters[comp]
                        adjcomp = jp.adjmatrix(orig, i-1, r, c)
                        permutations = jp.generate_permutations(adjcomp, len(adjcomp))

                        # check if tc has same adjacency matrix as arr
                        for perm in permutations:
                            if np.array_equal(perm,adjclus):
                                newclus = False
                                # if equivalent, update Yij
                                Y[ydict[clus]][ydict[comp]] += 1
                                break

            i = i + 1