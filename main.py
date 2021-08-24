import hitting_time
import networkx as nx
import json

def read_phenotypes(gspace_name):
  phenotypes_file = open("/data/phenotypes_"+gspace_name+".txt", "r")
  phenotypes = phenotypes_file.read().splitlines()
  phenotypes_file.close()

  return phenotypes

if __name__ == "__main__":
  gspace_name = "four_phen"

  gspace = nx.read_gml(gspace_name+".gml", label='id')
  phenotypes = read_phenotypes(gspace_name)

  genotype_networks_file = open("/data/gn_"+gspace_name+".json", "r")
  genotype_networks = json.load(genotype_networks_file)
  genotype_networks_file.close()

  final_phenotype = "Bbx"

  F = genotype_networks[final_phenotype]['nodes']

  gamma = 10**(-3)
  measurement_rate = 2000
  N = 500

  H = hitting_time.giveMeHamiltonian(gspace, gamma).toarray()
  T = hitting_time.estimateTransitionMatrix(H, measurement_rate, N)
  hitting_times = hitting_time.computeHittingTimes(T, F)

