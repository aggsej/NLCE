'''
This program loads data from topological cluster calculations and uses it to perform NLCE computations.
The first half looks at correlation as a function of temperature.
THe second half looks at correlation as a function of magnetic field.
'''

import pickle
with open('topologicalClusters', 'rb') as f:
    topodistinct = pickle.load(f)
with open('YMatrix', 'rb') as f:
    Y = pickle.load(f)
with open('linkedclusters', 'rb') as f:
    linkedclusters = pickle.load(f)      
import numpy as np
import matplotlib.pyplot as plt
import bisect
import math
import jp3fctn as jp
import itertools as itt

r = 4 # num rows
c = 4 # num col
N = r*c # number of sites
J = 1
orderexp = 8 # order of expansion

# constants
rsig = np.array([[0,1],[0,0]])
lsig = np.array([[0,0],[1,0]])
zsig = np.array([[1,0],[0,-1]])

sigx = np.array([[0,1],[1,0]],dtype=complex)
sigy = np.array([[0,-1j],[1j,0]],dtype=complex)
sigz = np.array([[1,0],[0,-1]],dtype=complex)

# construct sum<i,j> SzSz
def sumSzSz(adj):
    # construct sum<i,j> SzSz
    total = len(adj) # number of sites
    size = np.array([0,0])
    size[0] = np.size(np.dot(jp.gensigma(sigx, 1, total),jp.gensigma(sigx, 1, total)),0)
    size[1] = np.size(np.dot(jp.gensigma(sigx, 1, total),jp.gensigma(sigx, 1, total)),1)
    S = np.zeros(size,dtype=complex)
    for i in range(len(adj)):
        for j in range(i,len(adj[0])):
            if adj[i][j]: # for every pair of distinct adjacent sites
                S += np.dot(jp.gensigma(sigz, i+1, total),jp.gensigma(sigz, j+1, total))
    S = S/4 # pauli to spin matrix convesion
    return S

# construct sum<i,j> SxSx
def sumSxSx(adj,size):
    # construct sum<i,j> SxSx
    total = len(adj) # number of sites
    S = np.zeros(size,dtype=complex)
    for i in range(len(adj)):
        for j in range(i,len(adj[0])):
            if adj[i][j]: # for every pair of distinct adjacent sites
                S += np.dot(jp.gensigma(sigx, i+1, total),jp.gensigma(sigx, j+1, total))
    S = S/4 # pauli to spin matrix convesion
    return S

# construct sum<i,j> SxSx
def sumSySy(adj,size):
    # construct sum<i,j> SxSx
    total = len(adj) # number of sites
    S = np.zeros(size,dtype=complex)
    for i in range(len(adj)):
        for j in range(i,len(adj[0])):
            if adj[i][j]: # for every pair of distinct adjacent sites
                S += np.dot(jp.gensigma(sigy, i+1, total),jp.gensigma(sigy, j+1, total))
    S = S/4 # pauli to spin matrix convesion
    return S

# no dependence on J, but can be added in
def hamiltonian(adj, J,B,h):
    # construct hamiltonian
    total = len(adj) # number of sites
    size = np.array([0,0])
    size[0] = np.size(np.dot(jp.gensigma(sigx, 1, total),jp.gensigma(sigx, 1, total)),0)
    size[1] = np.size(np.dot(jp.gensigma(sigx, 1, total),jp.gensigma(sigx, 1, total)),1)
    H = sumSzSz(adj)+sumSySy(adj,size)+sumSxSx(adj,size)
    for i in range(len(adj)):
        H += (h/4)*jp.gensigma(zsig, i+1, total)
    eigvals, eigvec = np.linalg.eig(H) # eigenvalues of H
    prob = np.exp(-B*eigvals) # probability
    return H, eigvals, eigvec, prob

#return sz
def returnSz(adj):
    total = len(adj)
    ans = 0
    for i in range(len(adj)):
        ans += jp.gensigma(sigz, i+1, total)
    return ans

def returnSy(adj):
    total = len(adj)
    ans = 0
    for i in range(len(adj)):
        ans += jp.gensigma(sigy, i+1, total)
    return ans

def returnSx(adj):
    total = len(adj)
    ans = 0
    for i in range(len(adj)):
        ans += jp.gensigma(sigx, i+1, total)
    return ans

# given some set of values and some set of probabilities, outputs thermal expectation value
def thermalExpectation(prob,values):
    # find thermal expectaton value
    expect = np.sum(np.multiply(prob,values))/np.sum(prob) # expectation value
    return expect

# returns expectation value of <psi|matrix|psi>
def expectation(psi, matrix):
    result = np.zeros(len(psi),dtype=complex)
    for i in range(len(psi)):
        result[i] = np.dot(np.conj(psi),matrix[:,i])
    return np.dot(result,psi)

# wrapper for siteIndices function, returns indices of nonzero elements of cluster
def arrayIndices(r,c,cluster):
    # first find any nonzero entry
    found = False
    rnonzero = 0
    cnonzero = 0
    
    for i in range(r):
        if found:
            break;
        for j in range(c):
            if cluster[i][j] == 1:
                rnonzero = i
                cnonzero = j
                found = True
                break;
                    
    # find array of indices of nonzero 
    indnonzero = []
    vis = np.zeros((r,c),dtype=complex)
    indnonzero = jp.siteIndices(rnonzero, cnonzero, r, c, cluster, indnonzero, vis)
    return indnonzero

######

# energy/correlation at h=0
numtimes = 4 # how many orders to expand in
J = 1
Eoverall = [] 
Coverall = []

mstart = 3 # order of expansion to start with
for m in range(mstart,mstart+numtimes):
    # find number of clusters involved up to mth order
    numclus = 0
    for i in range(m):
        numclus += len(topodistinct[i])

    # independent variable
    T = np.arange(0.1,3,0.1)
    T = np.append(T,np.arange(3,8,0.5))
    #hfield = np.arange(0,5,0.2)
    # eigenvalues
    Evals = np.zeros(len(T),dtype=complex)
    Cvals = np.zeros(len(T),dtype=complex)
    kb = 1 #boltzmann constant
    B = 1/T #beta

    # for each value of temperture
    for hval in range(len(T)):
        # calculate property array
        prop = []
        correlation1 = []
        correlation2 = []

        # up to the order of expansion (i.e. m)
        for i in range(m):
            for clus,val in topodistinct[i].items():
                adj = jp.adjmatrix(linkedclusters[clus],i,r,c)
                if jp.numedges(adj) > 0:
                    H, eigvals, eigvec, prob = hamiltonian(adj, J,B[hval], 0)

                    propi = thermalExpectation(prob,eigvals)
                    #propi = expectation(eigvec[:,0],H)
                    prop.append(propi)
                    
                    Sz = returnSz(adj)
                    S = sumSzSz(adj)
                    values1 = np.zeros(len(eigvals),dtype=complex)
                    values2 = np.zeros(len(eigvals),dtype=complex)
                    for vec in range(len(eigvals)):
                        psi = eigvec[:,vec] # nth eigenvector
                        values2[vec] = expectation(psi, Sz) # append expectation value
                        values1[vec] = expectation(psi, S) # append expectation value
                    corri = thermalExpectation(prob,values1)
                    #corri = expectation(eigvec[:,0],S)
                    correlation1.append(corri)
                    #thing,eigvals,eigvec = returnSz(adj)
                    corri = thermalExpectation(prob,values2)
                    correlation2.append(corri)

                else:
                    prop.append(0)
                    correlation1.append(0) # not sure
                    correlation2.append(0) # not sure

        # calculate weight of cluster i 
        #print(numclus)
        for i in range(numclus):
            for j in range(i):
                prop[i] -= Y[i][j]*prop[j]
                correlation1[i] -= Y[i][j]*correlation1[j]
                correlation2[i] -= Y[i][j]*correlation2[j]
                
        #print(correlation)
        # the sum of the contributions to property P from all n-site topological clusters in the series (Cn)
        Sh = np.zeros(m,dtype=complex)
        Sc1 = np.zeros(m,dtype=complex)
        Sc2 = np.zeros(m,dtype=complex)
        counter = 0
        for i in range(m):
            # for clusters in ith order
            for clus, val in topodistinct[i].items():
                Lcn = topodistinct[i][clus]
                Sh[i] += Lcn*prop[counter]
                Sc1[i] += Lcn*correlation1[counter]
                Sc2[i] += Lcn*correlation2[counter]
                counter += 1
        # The property in the mth order of NLCE is the sum of Si's
        PmLh = 0
        PmLc1 = 0
        PmLc2 = 0
        
        #PmLc3 = 0
        for i in range(m):
            PmLh = PmLh + Sh[i]
            PmLc1 = PmLc1 + Sc1[i]
            PmLc2 = PmLc2 + Sc2[i]
        Evals[hval] = PmLh/6
        Cvals[hval] = PmLc1/32 #- PmLc2*PmLc2/16
    Eoverall.append(Evals)
    Coverall.append(Cvals)

fig, ax = plt.subplots()
ax.plot(T, Eoverall[0], 'r--', label='3rd order')
ax.plot(T, Eoverall[1], 'y--', label='4th order')
ax.plot(T, Eoverall[2], 'g--', label='5th order')
ax.plot(T, Eoverall[3], 'b--', label='6th order')
legend = ax.legend(loc='lower right', shadow=True, fontsize='large')
plt.xscale('log')
x1,x2,y1,y2 = plt.axis()
plt.axis((x1,x2,-0.2,0))

plt.title("Correlation per Site vs. Temperature for Canted Antiferromagnet")
plt.xlabel("T (log)")
plt.ylabel("Correlation")
plt.show()


'''
Now we look at varying magnetic field with temperature constant
'''


# energy/correlation at B = 1
numtimes = 3 # how many orders to expand in
J = 1
Eoverall = [] 
Coverall = []

mstart = 3 # order of expansion to start with
for m in range(mstart,mstart+numtimes):
    # find number of clusters involved up to mth order
    numclus = 0
    for i in range(m):
        numclus += len(topodistinct[i])

    # independent variable
    #T = np.arange(0.1,3,0.1)
    #T = np.append(T,np.arange(3,5,0.3))
    hfield = np.arange(0,5,0.2)
    # eigenvalues
    Evals = np.zeros(len(hfield),dtype=complex)
    Cvals = np.zeros(len(hfield),dtype=complex)
    kb = 1 #boltzmann constant
    B = 1 #beta

    # for each value of temperture
    for hval in range(len(hfield)):
        # calculate property array
        prop = []
        correlation1 = []
        correlation2 = []

        # up to the order of expansion (i.e. m)
        for i in range(m):
            for clus,val in topodistinct[i].items():
                adj = jp.adjmatrix(linkedclusters[clus],i,r,c)
                if jp.numedges(adj) > 0:
                    H, eigvals, eigvec, prob = hamiltonian(adj, J,B, hfield[hval])

                    propi = thermalExpectation(prob,eigvals)
                    #propi = expectation(eigvec[:,0],H)
                    prop.append(propi)
                    
                    S1 = returnSz(adj)
                    S2 = sumSzSz(adj)
                    
                    values1 = np.zeros(len(eigvals),dtype=complex)
                    values2 = np.zeros(len(eigvals),dtype=complex)
                    for vec in range(len(eigvals)):
                        psi = eigvec[:,vec] # nth eigenvector
                        values1[vec] = expectation(psi, S1) # append expectation value
                        values2[vec] = expectation(psi, S2) # append expectation value
                    corri = thermalExpectation(prob,values1)
                    #corri = expectation(eigvec[:,0],S)
                    correlation1.append(corri)
                    #thing,eigvals,eigvec = returnSz(adj)
                    corri = thermalExpectation(prob,values2)
                    correlation2.append(corri)

                else:
                    prop.append(0)
                    correlation1.append(0) # not sure
                    correlation2.append(0) # not sure

        # calculate weight of cluster i 
        #print(numclus)
        for i in range(numclus):
            for j in range(i):
                prop[i] -= Y[i][j]*prop[j]
                correlation1[i] -= Y[i][j]*correlation1[j]
                correlation2[i] -= Y[i][j]*correlation2[j]
                
        #print(correlation)
        # the sum of the contributions to property P from all n-site topological clusters in the series (Cn)
        Sh = np.zeros(m,dtype=complex)
        Sc1 = np.zeros(m,dtype=complex)
        Sc2 = np.zeros(m,dtype=complex)
        counter = 0
        for i in range(m):
            # for clusters in ith order
            for clus, val in topodistinct[i].items():
                Lcn = topodistinct[i][clus]
                Sh[i] += Lcn*prop[counter]
                Sc1[i] += Lcn*correlation1[counter]
                Sc2[i] += Lcn*correlation2[counter]
                counter += 1
        # The property in the mth order of NLCE is the sum of Si's
        PmLh = 0
        PmLc1 = 0
        PmLc2 = 0
        
        #PmLc3 = 0
        for i in range(m):
            PmLh = PmLh + Sh[i]
            PmLc1 = PmLc1 + Sc1[i]
            PmLc2 = PmLc2 + Sc2[i]
        Evals[hval] = PmLh
        Cvals[hval] = PmLc1/32 - PmLc2*PmLc2/(16*16)
    Eoverall.append(Evals)
    Coverall.append(Cvals)
    
fig, ax = plt.subplots()
ax.plot(T, Coverall[0], 'r--', label='3rd order')
ax.plot(T, Coverall[1], 'ys', label='4th order')
ax.plot(T, Coverall[2], 'g', label='5th order')
#ax.plot(T, Coverall[3], 'b', label='6th order')
legend = ax.legend(loc='lower right', shadow=True, fontsize='large')

#x1,x2,y1,y2 = plt.axis()
#plt.axis((x1,x2,-1,1))

plt.title("Correlation per Site vs. Temperature for Canted Antiferromagnet at h=0")
plt.xlabel("T")
plt.ylabel("Correlation")
plt.show()