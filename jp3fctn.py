# don't forget to cite: http://www.mafagafogigante.org/heaps-algorithm/
import numpy as np
import matplotlib.pyplot as plt
import bisect
import math
import itertools as itt

# constants
rsig = np.array([[0,1],[0,0]])
lsig = np.array([[0,0],[1,0]])
zsig = np.array([[1,0],[0,-1]])

###
# row(site,c): rerurns row
# col(site,c): returns column
# site(row,col,c): returns site number
# outOfBounds(row, col, rtot, ctot)
# maxCon1d(row, col, rtot, ctot, arr1, arr2, vis)
# maxCon2d(row, col, rtot, ctot, arr1, arr2, vis)
# siteIndices(row, col, rtot, ctot, arr1, arr2, vis): returns array of (row, col) of each occupied site



# return row of a site on grid
def row(site,c):
    return(site//c)

# return column of a site on grid
def col(site,c):
    return(site%c)

# return site number of a site on grid
def site(row,col,c):
    return(row*c+col)

# is (row, col) out of the lattice?
def outOfBounds(row, col, rtot, ctot):
    if row<0:
        return True
    if row > (rtot-1):
        return True
    if col<0:
        return True
    if col > (ctot-1):
        return True
    return False

# maximally connected subgraph of  1d arr1
def maxCon1d(row, col, rtot, ctot, arr1, arr2, vis):
    if outOfBounds(row,col,rtot,ctot):
        return 
    s = site(row,col,ctot)
    if vis[s] == 1:
        return
    vis[s] = 1
    if arr1[s] == 1:
        arr2[s] = 1
        maxCon1d(row-1,col,rtot,ctot,arr1,arr2, vis)
        maxCon1d(row,col+1,rtot,ctot,arr1,arr2, vis)
        maxCon1d(row+1,col,rtot,ctot,arr1,arr2, vis)
        maxCon1d(row,col-1,rtot,ctot,arr1,arr2, vis)
    return

# returns maximally connected subarray of 2d lattice
def maxCon2d(row, col, rtot, ctot, arr1, arr2, vis):
    if outOfBounds(row,col,rtot,ctot):
        return
    s = site(row,col,ctot)
    if vis[row][col] == 1:
        return
    vis[row][col] = 1
    if arr1[row][col] == 1:
        arr2[row][col] = 1 
        maxCon2d(row-1,col,rtot,ctot,arr1,arr2, vis)
        maxCon2d(row,col+1,rtot,ctot,arr1,arr2, vis)
        maxCon2d(row+1,col,rtot,ctot,arr1,arr2, vis)
        maxCon2d(row,col-1,rtot,ctot,arr1,arr2, vis)
    return arr2


# returns array of site numbers of points on 2d lattice
def siteIndices(row, col, rtot, ctot, arr1, arr2, vis):
    if outOfBounds(row,col,rtot,ctot):
        return
    s = site(row,col,ctot)
    if vis[row][col] == 1:
        return
    vis[row][col] = 1
    if arr1[row][col] == 1:
        arr2.append((row,col))
        siteIndices(row-1,col,rtot,ctot,arr1,arr2, vis)
        siteIndices(row,col+1,rtot,ctot,arr1,arr2, vis)
        siteIndices(row+1,col,rtot,ctot,arr1,arr2, vis)
        siteIndices(row,col-1,rtot,ctot,arr1,arr2, vis)
    return arr2

# finds every permutation of possible clusters
def genAllClus(N, order):
    clusters = [] # all clusters
    numclusters = np.zeros(order) # keep track of each order
    allSites = np.array(range(N))
    i = 1
    while i < order + 1:
        # choose i of those points to be one and store the arrays
        combs = itt.combinations(allSites, i)
        for points in combs:
            arr = np.zeros(N)
            for idx in points:
                arr[idx] = 1
            clusters.append(arr)
        numclusters[i-1] = len(clusters)
        i = i + 1
    return clusters, numclusters

# finds every permutation of possible subclusters of 2d array
def genAllSubClus(indnonzero, r, c):
    subClusters = [] # all clusters
    numclusters = np.zeros(order) # keep track of each order
    order = len(indnonzero)
    i = 1
    while i < order + 1: 
        # choose i of those points to be one and store te arrays
        combs = itt.combinations(indnonzero, i)
        for points in combs:
            arr = np.zeros(N)
            for idx in points:
                row = indnonzero[idx][0]
                col = indnonzero[idx][1]
                arr[row][col] = 1
            subClusters.append(arr)
        numclusters[i-1] = len(subClusters)
        i = i + 1
    return subClusters, numclusters


# takes an array of 1d clusters and returns array of all 2d linked clusters
def onlyLinkedClus(clusters, numclusters, r, c, order):
    numlinclus = np.zeros(order) # number of total linked clusters up to each order                      
    linkedclusters=[]
    current = 0
    for i in range(len(clusters)):
        arr1 = clusters[i]
        arr2 = np.zeros(arr1.size)
        for j in range(arr1.size):
            if arr1[j] == 1:
                nonzero = j
                break
        vis = np.zeros(arr1.size)
        maxCon1d(row(nonzero, c), col(nonzero,c) , r, c, arr1, arr2, vis)
        if np.array_equal(arr1, arr2):
            linkedclusters.append(np.reshape(arr2,(r,c))) # changes array from 1d to 2d

        # index of linked cluster in each incremented order of expansion
        if i + 2 > numclusters[current]:
            numlinclus[current] = len(linkedclusters)
            current = current + 1
    return linkedclusters, numlinclus

# takes array of linked clusters and outputs symmetrically distinct clusters
def onlySymDist(linkedclusters, numlinclus, r, c, printoutput):
    distinct = [] # symmetrically distinct
    j = 0 # count which cluster we're looking at
    zerc = np.zeros(c)
    zerr = np.zeros(r)
    current = 0
    numsymclus = np.zeros(len(numlinclus))

    # for each order
    for i in range(len(numlinclus)):
        orderclus = {j: 1} # add first cluster in each order
        j = j + 1

        # for each cluster in linked clusters
        while j < numlinclus[i]:
            # generate all symmetrically related clusters
            cg = genclusters(linkedclusters[j])

            newclus = True # is linkedclusters[j] a distinct cluster?

            # for each generated cluster
            for clus in cg:
                if not newclus:
                    break;

                # check whether it is in distinct
                newclus, orderclus = cluscomp(clus, orderclus, linkedclusters)

                # try to find translation vector
                # translate to the right
                # column is blank
                currcol = c - 1
                while np.array_equal(clus[:,currcol], zerr) and newclus:
                    # temp = translated clus
                    temp = np.roll(clus,c-currcol,axis=1)

                    # check whether clus is in distinct
                    newclus, orderclus = cluscomp(temp, orderclus, linkedclusters)

                    currcol = currcol - 1

                # translate to the left
                # column is blank
                currcol = 0
                while np.array_equal(clus[:,currcol], zerr) and newclus:
                    # temp = translated clus
                    temp = np.roll(clus,-currcol-1,axis=1)

                    # check whether translated clus is in distinct
                    newclus, orderclus = cluscomp(temp, orderclus, linkedclusters)

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
                    newclus, orderclus = cluscomp(temp, orderclus, linkedclusters)

                    # try to find translation vector
                    # translate to the right
                    # column is blank
                    currcol = c - 1
                    while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                        temp2 = np.roll(temp,c-currcol,axis=1)


                        # check whether clus is in distinct
                        newclus, orderclus = cluscomp(temp2, orderclus, linkedclusters)

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
                        newclus, orderclus = cluscomp(temp2, orderclus, linkedclusters)

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
                    newclus, orderclus = cluscomp(temp, orderclus, linkedclusters)

                    # try to find translation vector
                    # translate to the right
                    # column is blank
                    currcol = c - 1
                    while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                        temp2 = np.roll(temp,c-currcol,axis=1)

                        # check whether clus is in distinct
                        newclus, orderclus = cluscomp(temp2, orderclus, linkedclusters)

                        currcol = currcol - 1

                    # translate to the left
                    # column is blank
                    currcol = 0
                    while np.array_equal(temp[:,currcol], zerr) and newclus:
                        # temp = translated clus
                        temp2 = np.roll(temp,-currcol-1,axis=1)
                        for k in range(currcol+1):
                            temp2[:,c-1-k] = 0

                        # check whether translated clus is in distinct
                        newclus, orderclus = cluscomp(temp2, orderclus, linkedclusters)

                        currcol = currcol + 1            

                    currow = currow - 1

            # if not (i.e. if it is new)
            if newclus:
                orderclus[j] = 1
            j = j + 1
        distinct.append(orderclus)
        if printoutput:
            print(distinct[i])
        if j + 2 > numlinclus[current]:
            numsymclus[current] = len(orderclus)
            current = current + 1
    for i in range(len(numsymclus) - 1):
        numsymclus[i+1] += numsymclus[i]
    return distinct, numsymclus

##### symmetries #####
def identity(arr1):
    return arr1

def rotate90(arr1):
    return np.rot90(arr1)

def rotate180(arr1):
    return np.rot90(arr1,2)

def rotate270(arr1):
    return np.rot90(arr1,3)

# reflect along x = 0
def refx0(arr1):
    return np.fliplr(arr1)

# reflect along y = 0
def refy0(arr1):
    return np.flipud(arr1)

# reflect along y=x
def refyx(arr1):
    (r,c)=arr1.shape
    # flatten array, reverse it, reshape it, transpose it
    return np.transpose(np.reshape(list(reversed(arr1.flatten())),(r,c)))

# reflect along y=-x
def refynx(arr1):
    return np.transpose(arr1) 
#### end symmetries ####

#format output of array
def printf(arr1):
    for arr in arr1:
        for ele in arr:
            print(str(int(ele)), end=" ")
        print()
        
# compare clus to all clusters in distinct array of the same order
def cluscomp(clus, orderclus, linkedclusters):
    for comp, val in orderclus.items():
        if np.array_equal(clus,linkedclusters[comp]):
            # if yes, increase multiplicity by 1
            orderclus[comp] = val + 1
            return False, orderclus
    return True, orderclus

# generate related clusters
def genclusters(arr1):
    cg = [] # clusters generated
    cg.append(identity(arr1))
    cg.append(rotate90(arr1))
    cg.append(rotate180(arr1))
    cg.append(rotate270(arr1))
    cg.append(refx0(arr1))
    cg.append(refy0(arr1))
    cg.append(refyx(arr1))
    cg.append(refynx(arr1))
    return cg

# generates adjacency matrix
def largeadjmatrix(clus, r, c):
    vis = np.zeros((r,c)) 
    adj = np.zeros((r*c,r*c))
    for i in range(r):
        for j in range(c):
            if clus[i][j] == 1:
                rnonzero = i
                cnonzero = j
                break
    neighborsite(rnonzero,cnonzero,r,c,clus,adj,site(rnonzero,cnonzero,c),vis)
    for i in range(r*c):
        adj[i][i] = 0
    return adj

def neighborsite(row, col, rtot, ctot, clus, adj, prev, vis):
    if outOfBounds(row, col, rtot, ctot):
        return
    if vis[row][col] == 1:
        return
    vis[row][col] = 1
    if clus[row][col] == 1:
        scur = site(row,col,c)
        adj[scur][prev] = 1
        adj[prev][scur] = 1
        neighborsite(row-1,col,rtot,ctot,clus,adj, scur, vis)
        neighborsite(row,col+1,rtot,ctot,clus,adj, scur, vis)
        neighborsite(row+1,col,rtot,ctot,clus,adj, scur, vis)
        neighborsite(row,col-1,rtot,ctot,clus,adj, scur, vis)
    return

# number of edges in adjacency matrix
def numedges(adj):
    tot = 0
    for row in adj:
        for ele in row:
            if ele == 1:
                tot += 1
    return tot/2

# returns adjacency matrix of cluster
def adjmatrix(clus, curorder2, r, c):
    indices = dict()
    count = 0
    curorder = curorder2 + 1 # since counting starts from 0
    for i in range(r):
        for j in range(c): 
            if clus[i][j]==1:
                indices[site(i,j,c)] = count
                count += 1
    adj = np.zeros((curorder,curorder))
    for node in indices:
        currow = row(node,c)
        curcol = col(node,c)

        # check left
        if curcol>0:
            if (node - 1) in indices:
                adj[indices[node]][indices[node - 1]] =1

        # check right
        if curcol < c - 1:
            if (node + 1) in indices:
                adj[indices[node]][indices[node + 1]] =1

        # check up
        if currow > 0:
            if (node - c) in indices:
                adj[indices[node]][indices[node - c]] =1

        # check down
        if currow < r - 1:
            if (node + c) in indices:
                adj[indices[node]][indices[node + c]] =1

    return adj

# returns key of cluster based on adjacency matrix for topological comparison
def key(adj):
    key = 0
    for i in range(np.shape(adj)[0]):
        for j in range(np.shape(adj)[1]):
            key += (adj[i][j])*(((i+1)*(j+1))**((i+1)+(j+1)))
    return key
            
# # # # # topologically distinct # # # # #
def onlyTopologicallyDistinct(symdistinct, linkedclusters, r, c, printoutput):
    topodistinct = [] # topo distinct

    # for each order
    for curorder in range(len(symdistinct)):
        first = True
        for matrix, weight in symdistinct[curorder].items():
            if first: 
                orderclus = {matrix: weight} # add first cluster in each order
                first = False
            else:
                newclus = True # is linkedclusters[j] a distinct cluster?
                adjclus = adjmatrix(linkedclusters[matrix], curorder, r, c)

                # for each matrix in list of topologically distinct matrices of current order
                for comp, val in orderclus.items():
                    if not newclus:
                        break
                    orig = linkedclusters[comp]
                    adjcomp = adjmatrix(orig, curorder, r, c)
                    permutations = generate_permutations(adjcomp, len(adjcomp))
                    
                    for perm in permutations:
                        if np.array_equal(perm,adjclus):
                            newclus = False
                            orderclus[comp] = weight + val
                            break
                if newclus:
                    orderclus[matrix] = weight
        if printoutput:
            print(orderclus)
            print(len(orderclus))
        topodistinct.append(orderclus)
    return topodistinct 

# swaps row i with row j and column i with column j
def swapIndices(arr1, i, j, defensiveCopy):
    if defensiveCopy:
        a1 = np.copy(arr1)
    else:
        a1 = arr1
    temp = np.copy(a1[:, i])
    a1[:, i] = a1[:, j]
    a1[:, j] = temp
    temp = np.copy(a1[i, :])
    a1[i, :] = a1[j, :]
    a1[j, :] = temp
    return a1

# this code was modified from: http://www.mafagafogigante.org/heaps-algorithm/
def generate_permutations(elements, n):
    # As by Robert Sedgewick in Permutation Generation Methods
    c = [0] * n
    yield elements
    i = 0
    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                swapIndices(elements, 0, i, False)
            else:
                swapIndices(elements, c[i], i, False)
            yield elements
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1

# calculates sum of L(c) for n sites
def latticeConstant(sites):
    L = [1,2,6,19,63,216,760,2725,9910,36446,135268,505861,1903890,7204874,27394666,104592937]
    return L[sites-1]

# generate spin matrices
def gensigma(new, site, total):
    sigma = np.identity(2**(site-1),dtype=complex)
    sigma = np.kron(sigma,new)
    if total > site:
        sigma = np.kron(sigma,np.identity(2**(total-site)))
    return sigma

# expectation value of Hamiltonian
def hexpect(adj, J, B):
    total = len(adj) # number of sites
    H = 0 # Hamiltonian
    for i in range(len(adj)):
        for j in range(i,len(adj[0])):
            if adj[i][j]: # for every pair of distinct adjacent sites
                H += 2*np.dot(gensigma(rsig, i+1, total),gensigma(lsig, j+1, total))
                H += 2*np.dot(gensigma(lsig, i+1, total),gensigma(rsig, j+1, total))
                H += np.dot(gensigma(zsig, i+1, total),gensigma(zsig, j+1, total))
    H = H*J/4
    eigvals, eigvec = np.linalg.eig(H) # eigenvalues of H
    prob = np.exp(-B*eigvals) # probability
    hexpect = np.sum(np.multiply(prob,eigvals))/np.sum(prob) # expectation value
    return hexpect

# returns L
'''
def genLvalues():
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
    return(Lvalues) '''
