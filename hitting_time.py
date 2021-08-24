import numpy as np
import networkx as nx
from scipy.linalg import expm
from scipy.stats import expon

def giveMeHamiltonian(G, gamma, typeMatrix="laplacian"):
  # typeMatrix could be "adjacency" to use the adjacency matrix based Hamiltoninan 
  # or "laplacian" (default) to use the Laplacian matrix based Hamiltonian
  # gamma is the mutation rate 

  if typeMatrix == "adjacency":
    A = nx.adjacency_matrix(G)
  else:
    A = nx.laplacian_matrix(G)

  H = -gamma* A

  return H

# hitting time (mean first passage time) to F of a Markov chain 
def computeHittingTimes(transition_matrix, F):

  A = np.delete(transition_matrix, F, axis=0)
  A = np.delete(A, F, axis = 1)
  A = A - np.identity(A.shape[0])

  b = np.array([-1 for i in range(transition_matrix.shape[0] - len(F))])

  x = np.linalg.solve(A, b)

  for i in F:
    # F it is supposed to be ordered
    x = np.insert(x, i, 0)

  return x

def estimateTransitionMatrix(H, measurement_rate, N):
  # H is the Hamiltonian of the quantum walk
  # N is the number of realizations of the exponential random variable that determines the time at which we measure the quantum walk
  
  M = H.shape[0]

  times = expon.rvs(scale=measurement_rate, size=N)
  
  expected_matrix = np.zeros((M, M))

  for t in times:
    evolution_matrix = expm(-1j*t*H)
    transition_matrix = np.zeros((M, M))

    for i in range(M):
      for j in range(M):
        transition_matrix[i, j] = abs(evolution_matrix[i,j])**2

    expected_matrix += transition_matrix
  
  expected_matrix = expected_matrix/N

  return expected_matrix