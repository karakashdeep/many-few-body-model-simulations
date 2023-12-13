import numpy as np
import sympy as sp
from math import comb
import matplotlib.pyplot as plt
import networkx as nx
from scipy.sparse.linalg import eigsh, inv
from scipy.linalg import block_diag
from scipy.io import mmread, mmwrite
from itertools import product
from scipy.integrate import trapezoid as trapz
from scipy.linalg import bandwidth
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee as RCM
plt.rcParams["figure.figsize"] = (20,10)
#----------------------------------------------------------


#generates base10 representation of Permutation Ordered Binary (POB)
def pob2number(L):

    L.reverse()
    sum = 0
    for i in range(len(L)):
        temp = 0
        for j in range(i+1):
            temp += L[j]
        sum += L[i]*comb(i,temp)
    return sum
#----------------------------------------------------------


#generates POB representation of base10 number
def number2pob(n,r,v):

    s = ""
    k = r
    j = n

    while k >= 1:
        
        while True:
            j = j - 1
            p = comb(j, k)
            if v >= p:
                v = v - p
                bj = 1
                #print("k, j, p = ", str(k), str(j), str(p))
            else:
                bj = 0
            s += str(bj)
            if bj == 1:
                break
        k = k - 1


    if j >= 0:
        while j >= 0:
            j = j - 1
            bj = 0
            s += str(bj)


    s_list = [int(digit) for digit in s][:len(s)-1]
    return s_list

#----------------------------------------------------------


#takes a state and generates it's corresponding "hoppable" states.

def hopping_state_list(s_list):
    N = len(s_list)  # Number of lattice sites

    # Perform the hopping operation
    state_list = []  # Initialize a list to store states

    for i in range(N):
        if s_list[i] == 1:  # If there is an electron at site i
            # Hop to the next site (periodic boundary)
            sign = 1  # Sign for hopping to the next site
            new_state = s_list.copy()  # Create a copy of the current state
            new_state[(i + 1) % N] = 1
            new_state[i] = 0
            if new_state.count(1) == s_list.count(1):
                state_list.append((sign, new_state))

            # Hop to the previous site (periodic boundary)
            sign = 1  # Sign for hopping to the previous site
            new_state = s_list.copy()  # Create a copy of the current state
            new_state[(i - 1) % N] = 1
            new_state[i] = 0
            if new_state.count(1) == s_list.count(1):
                state_list.append((sign, new_state))

    return state_list
##########

#----------------------------------------------------------

#generates matrix representation of Hopping term of Hamiltonian

def hopping_term_matrix(n,r,t=1,justhalf=False):

    K = np.zeros((comb(n,r),comb(n,r)))
    for i in range(comb(n,r)):
        s_list_1 = hopping_state_list(number2pob(n,r,i))
        for j in s_list_1:
            temp = pob2number(j[1])
            K[i,temp] = j[0]
    if justhalf == False:
        return t*(np.kron(K,np.identity(comb(n,r)))+np.kron(np.identity(comb(n,r)),K))
    else:
        return t*K
#----------------------------------------------------------
#generates matrix representation of Coloumbic interaction term of Hamiltonian

def interaction_term_matrix(n,r,U=-1.5):
    l = comb(n,r)
    M = np.zeros((l,l))
    for i in range(l):
        s=number2pob(n,r,i)
        num_list = [i for i in range(comb(n,r))]
        mn_list = [(i,j) for i,j in product(num_list,num_list)]
        x_up = number2pob(n,r,mn_list[i][0])
        x_down = number2pob(n,r,mn_list[i][1])
        x_occupancy = np.dot(x_up,x_down)
        M[i][i] = x_occupancy
    return U*M
#----------------------------------------------------------

#Misc methods to diagonalize matrix
def diagonalize_matrix(matrix):
    M = sp.Matrix(matrix)
    P,D = M.diagonalize()
    return np.array(D).astype(np.complex64)

def diagonalize_matrix_approx(matrix, num_eigenvalues=36):
    # Diagonalize the sparse matrix using eigsh (sparse eigenvalue solver)
    eigenvalues, eigenvectors = eigsh(matrix, k=num_eigenvalues)

    # Sort eigenvalues and eigenvectors
    sorted_indices = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[sorted_indices]
    eigenvectors = eigenvectors[:, sorted_indices]

    # Create the diagonal matrix from eigenvalues
    diagonal_matrix = np.diag(eigenvalues)

    return diagonal_matrix
#----------------------------------------------------------
#regular spectrum plotter for hubbard hamiltonian
def spec_plotter(n,r,t=-1,u=0,mu=0, plot=False):

    K = hopping_term_matrix(n,r,t,mu)
    U = interaction_term_matrix(n,r,u)
    H = K + U
    sol = np.linalg.eigh(H)
    H_1 = diagonalize_matrix_approx(H,num_eigenvalues=H.shape[0])

    x=np.linspace(min(sol[0]-0.5),max(sol[0]+0.5),5000)
    print("min / max eigs:",min(sol[0]),"/",max(sol[0]))
    #x=np.linspace(-2.5,2.5,100000)
    y=np.array([0 for i in range(len(x))])
    for i in range(len(x)):
        temp = np.identity((comb(n,r)**2))*(x[i]+0.1j)
        y[i] = -(1/np.pi)*np.imag(np.trace(np.linalg.inv((temp - H_1))))
        
    if plot == True:
        plt.plot(x,y,"-")
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'A($\omega$)')
        #print(y)
        plt.show()
       
    return (x,y)
#----------------------------------------------------------
#spec plotter that only takes in a matrix
def spec_plotter_custom(M, plot=True):
    
    H_1 = M
    sol = np.linalg.eigh(H_1)
    l = len(M)
    #H_1 = diagonalize_matrix_approx(H,num_eigenvalues=H.shape[0])
    x=np.linspace(-4.5,4.5,500)
    #print("min / max eigs:",min(sol[0]),"/",max(sol[0]))
    #x=np.linspace(-2.5,2.5,100000)
    y=np.array([0 for i in range(len(x))])
    for i in range(len(x)):
        temp = np.identity(l)*(x[i]+0.1j)
        y[i] = -(1/np.pi)*np.imag(np.trace(np.linalg.inv((temp - H_1))))
    if plot == True:
        print("Trapezoidal=",trapz(y,x))
        plt.plot(x,y,"-")
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'A($\omega$)')
        #plt.legend(["Exact Diagonalization"])
        #print(y)
    sum=0
    return (x,y)

#----------------------------------------------------------



#graphs.....



#Plots the total fock space lattice
def plot_fock_lattice(n,r):
    K = hopping_term_matrix(n,r,t=1,justhalf=True)
    G = nx.Graph()
    l = [i for i in range(comb(n,r))]
    G.add_nodes_from(l)
    for i  in l:
        for j in l:
            if K[i,j] !=0:
                G.add_edge(i,j)
    #nx.draw(G,with_labels = True)
    return G
    #plt.show()

#----------------------------------------------------------
#Duplicate of above that converts adjecent matrix to graph
def plot_adj_2_graph(K):
    G = nx.Graph()
    l = [i for i in range(len(K))]
    G.add_nodes_from(l)
    for i  in l:
        for j in l:
            if K[i,j] !=0:
                G.add_edge(i,j)
    nx.draw(G,with_labels = True)
    return G

#----------------------------------------------------------
#only plots the maps from one sector to next i.e |n_l, n_r> ---> |n_l +/- 1, n_r -/+ 1>
def plot_fock_sector_lattice_maps(n,r):
    K = hopping_term_matrix(n,r,t=1,justhalf=True)
    G = nx.Graph()
    l = [i for i in range(comb(n,r))]
    G.add_nodes_from(l)
    for i  in l:
        for j in l:
            if K[i,j] !=0:
                state1 = number2pob(n,r,i)
                state2 = number2pob(n,r,j)
                ls1 = sum(state1[:r])
                rs1 = sum(state1[r:])
                ls2 = sum(state2[:r])
                rs2 = sum(state2[r:])
                if (ls1 == ls2-1 and rs1 == rs2+1) or (ls1 == ls2+1 and rs1 == rs2-1):
                    G.add_edge(i,j)
    nx.draw(G,with_labels = True)
    #plt.show()
#----------------------------------------------------------

#Plots each disconnected sector [N,0], [N-1,1], .... [0,N] i.e N+1 sectors for half filling case
def plot_fock_sector_lattice(n,r):
    K = hopping_term_matrix(n,r,t=1,justhalf=True)
    G = nx.Graph()
    l = [i for i in range(comb(n,r))]
    G.add_nodes_from(l)
    for i  in l:
        for j in l:
            if K[i,j] !=0:
                state1 = number2pob(n,r,i)
                state2 = number2pob(n,r,j)
                ls1 = sum(state1[:r])
                rs1 = sum(state1[r:])
                ls2 = sum(state2[:r])
                rs2 = sum(state2[r:])
                if (ls1 == ls2 and rs1 == rs2):
                    G.add_edge(i,j)
    nx.draw(G,with_labels = True)
    #plt.show()

def plot_degree_dist(G):
    degrees = [G.degree(n) for n in G.nodes()]
    nodes = [n for n in G.nodes()]
    plt.plot(nodes,degrees,"-")
    plt.show()
def plot_degree_dist_alt_1(G,m):
    degree_freq = nx.degree_histogram(G)
    degrees = range(len(degree_freq))
    plt.loglog(degrees[m:], degree_freq[m:],'go-') 
    plt.xlabel('Degree')
    plt.ylabel('Frequency')
    plt.show()
def plot_fock_lattice_degree_dist(n,r):
    K = hopping_term_matrix(n,r,t=1,justhalf=True)
    G = nx.Graph()
    l = [i for i in range(comb(n,r))]
    G.add_nodes_from(l)
    for i  in l:
        for j in l:
            if K[i,j] !=0:
                G.add_edge(i,j)
    plot_degree_dist(G)
    
    #plt.show()

#----------------------------------------------------------
#generates the permutation matrix given a permutation...manually
def perm_matrix(mat,perm):
    perm = list(perm)
    #perm.reverse()

    l = len(perm)

    orig = [i for i in range(l)]
    P = np.zeros((l,l))
    for i in orig:
        P[perm[i],i] = 1

    
    
    return np.matmul(P.T,np.matmul(mat,P))



def permutation_matrix(perm):
    perm = list(perm)
    #perm.reverse()

    l = len(perm)

    orig = [i for i in range(l)]
    P = np.zeros((l,l))
    for i in orig:
        P[perm[i],i] = 1

    
    
    return P
#Swaps...
def swap_operator(M,i,j):
    C = np.zeros(M.shape)
    C[i,j] = 1
    C[j,i] = 1
    for k in range(C.shape[0]):
        if k != i and k != j:
            C[k,k] = 1
    return np.matmul(C.T,np.matmul(M,C))
#takes in perm and matrix, to give permuted matrix, in house scipy ...

def perm_on_matrix(M,perm):
    G = csr_matrix(M)
    return csr_matrix.toarray(G[perm])

swap = lambda M,i,j : swap_operator(M,i,j)
submatrix = lambda M,i,j : M[tuple([slice(i,j)]*M.ndim)]


#-------------------------------------------------------------
def wolfram_rcm(mat):
    np.save('temp_matrix.npy', mat)
    # Bandwidth reduction-- in Wolfram Mathematica...
    from wolframclient.evaluation import WolframLanguageSession
    from wolframclient.language import wl, wlexpr
    session = WolframLanguageSession()

    session.evaluate(wlexpr('''Needs["GraphUtilities`"]
                            << NumPyArray`
                            M = ReadNumPyArray["temp_matrix.npy"]
                            {r,c} = MinimumBandwidthOrdering[M,Method->"RCM"]
                            Export["Mat.mtx", M[[r,c]], "MTX"]
                            Export["r_perm.csv", r, "CSV"]
                            Export["c_perm.csv", c, "CSV"]
                            ClearAll["Global`*"]

    '''))
    rcm_mat = mmread("Mat.mtx")
    r_perm = np.loadtxt("r_perm.csv")
    r_perm -= [1 for _ in range(len(r_perm))]
    c_perm = np.loadtxt("c_perm.csv")
    c_perm -= [1 for _ in range(len(c_perm))]
    print("r=",r_perm)
    print("c=",c_perm)
    ##plot_adj_2_graph(permutation_matrix(r_perm))
    ##rcm_mat_2 = perm_on_matrix(mat,perm)
    return rcm_mat

#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#Tridiagonalization functions

#Generates pattern representing block diagonal pivots. Eg. l=4-->[0, 1, 5, 11, 15, 16]
def pattern(l):
    pat = [1]
    m = int(l/2)
    for i in range(m-1):
        pat.append(4+(i*4))
    temp = reversed(pat)
    pat.append(4+(i*4)+2)
    #print("i=",i)
    #pat.append(temp)
    for i in temp:
        pat.append(i)
    #print(pat,"=",sum(pat))
    #block pattern
    block_pat = [0]
    count=0
    for i in pat:
        count += i
        block_pat.append(count)
    return block_pat

#-------------------------------------------

#The tridiagonal representation returned as a 3-tuple (A,B,C) where each element of the tuple is a list of block matrices of length n+1,n,n respectively
def block_partition(Mat):
    n = np.shape(Mat)[0]**0.5
    block_pat = pattern(n)
    A = []
    B = []
    C = []
    e_i = 0 #initial element
    e_m = 0 #middle element
    e_f = 0 #final element


    ###########################
    #COLLECTS DIAGONAL BLOCKS A
    for i in range(len(block_pat)-1):
        e_i = block_pat[i]
        e_f = block_pat[i+1]
        A.append(Mat[e_i:e_f,e_i:e_f])
        
    


    ###########################
    #COLLECTS upper DIAGONAL BLOCKS B
    for i in range(len(block_pat)-2):
        e_i = block_pat[i]
        e_m = block_pat[i+1]
        e_f = block_pat[i+2]
        B.append(Mat[e_i:e_m,e_m:e_f])
        


    ###########################
    #COLLECTS lower DIAGONAL BLOCKS C
    for i in range(len(block_pat)-2):
        e_i = block_pat[i]
        e_m = block_pat[i+1]
        e_f = block_pat[i+2]
        C.append(Mat[e_m:e_f,e_i:e_m])
    

    return (A,B,C)

#---------------------------------------------
#Forward and backward recursions...
#returns the green's function "G_F_N" as a diagonal blocked matrix, for a given frequency f...
def recursive_G_F(M_tri,f,flags=False):
    diag_list = M_tri[0]            #  i = 0,...,2*m 
    u_diag_list = M_tri[1]          #  i = 0,...,2*m-1
    l_diag_list = M_tri[2]          #  i = 0,...,2*m-1
    G_F_list = []

    n = len(u_diag_list)    #  =2*m
    
    def G_0(i):
        block = diag_list[i]
        l = np.shape(block)[0]
        return np.linalg.inv( np.identity(l)*(f  + 0.1j) - block )
    
    def u_tau(i):
        return u_diag_list[i]
    
    def l_tau(i):
        return l_diag_list[i]
    
    #first block
    G_F_list.append(G_0(0))
    
    
    for i in range(n):
        
        G_F = np.linalg.inv(np.linalg.inv(G_0(i+1)) - np.matmul(l_tau(i),np.matmul(    G_F_list[i]    ,u_tau(i))))
        if flags==True:
            print(np.shape(G_0(i+1)),"-",np.shape(l_tau(i)),"X",np.shape(G_F_list[i]),"X",np.shape(u_tau(i)),"=",np.shape(G_F))
        G_F_list.append(G_F)
    if flags==True:
        print("---------------")
    
    #print(np.imag(G_0(10)))
    #print(np.imag(G_F_list[len(G_F_list)-1]))


    G_N_list = [0 for i in range(n+1)]
    G_N_list[n] = G_F_list[len(G_F_list)-1]
    for j in range(n-1,-1,-1):
        temp1 = np.identity(np.shape(u_tau(j))[0])
        temp2 = np.matmul(u_tau(j), np.matmul(G_N_list[j+1],np.matmul(l_tau(j),G_F_list[j])))
        G_N = np.matmul(G_F_list[j],temp1+temp2)
        G_N_list[j] = G_N
        if flags==True:
            print(np.shape(G_F_list[j]),"X",np.shape(np.identity(np.shape(u_tau(j))[0])),"+",np.shape(u_tau(j)),"X",np.shape(G_N_list[j+1]),"X",np.shape(l_tau(j)),"X",np.shape(G_F_list[j]),"=",np.shape(G_N))
    if flags==True:
        print("---------------")

    total_G_N = 0
    for i in G_N_list:
        total_G_N = block_diag(total_G_N,i)



    return (total_G_N[1:,1:])

#----------------------------------------------
def rgf_2_spectra(H_matrix, plot=True):
    rcm_H_matrix = wolfram_rcm(H_matrix)
    M_tri = block_partition(rcm_H_matrix)
    
    #frequency space
    x = np.linspace(-4.5,4.5,500)
    y=[]
    for i in x:
        y.append(-(1/np.pi)*np.trace(np.imag(recursive_G_F(M_tri,i))))
    print("Trapezoidal=",trapz(y,x))
    if plot==True:
        plt.plot(x,y)
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'A($\omega$)')
        plt.legend(["Calculated via RCM-reduced F-RGF"])
        plt.show()
    return (x,y)

#----------------------------------------------
def two_particle_spec(n,plot=True):

    def rel_1(n,r,val):
        state = number2pob(n,r,val)
        if (state[0]==1 and state[n-1]==1):
            return True
        for i in range(n-1):
            if (state[i]==1 and state[i+1]==1):
                return True
        else:
            return False
    H = hopping_term_matrix(n,2,1,justhalf=True)+interaction_term_matrix(n,2,-3)


    sol = np.linalg.eigh(H)
    l = len(H)
    #H_1 = diagonalize_matrix_approx(H,num_eigenvalues=H.shape[0])
    H_1 = H
    x=np.linspace(min(sol[0]-0.5),max(sol[0]+0.5),1000)
    #print("min / max eigs:",min(sol[0]),"/",max(sol[0]))
    #x=np.linspace(-2.5,2.5,100000)
    y=np.array([0 for i in range(len(x))])
    for i in range(len(x)):
        temp = np.identity(l)*(x[i]+0.01j)
        mat_custom_trace = 0
        for j in range(H_1.shape[0]):
            
            if rel_1(n,2,j):
                mat_custom_trace +=np.linalg.inv((temp - H_1))[j,j]
                #print("j=",j)





        y[i] = -(1/np.pi)*np.imag(mat_custom_trace)
    if plot == True:
        print("Trapezoidal=",trapz(y,x))
        plt.plot(x,y,"-")
        plt.xlabel(r'$\omega$')
        plt.ylabel(r'A($\omega$)')
        #print(y)
    sum=0
    return (x,y)



