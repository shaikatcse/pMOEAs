pAGE.F: pAGE offline (population size: 30, evaluations: 12000)
pAGE.O: pAGE online (population size: 30, evaluations: 12000)
AGE.1: AGE (population size: 30, evaluations: 12000)
AGE.2: AGE (population size: 100, evaluations: 12000)
AGE.3: AGE (population size: 100, evaluations: 24000)
pAGE.F.Max: pAGE offline (population size: 30, evaluations: 24000)
pAGE.O.Max: pAGE online (population size: 30, evaluations: 24000)


pNSGAII: pNSGAII (population size: 30, evaluations: 12000)
NSGAII.1: NSGAII (population size: 30, evaluations: 12000)
NSGAII.2: NSGAII (population size: 100, evaluations: 12000)
NSGAII.3: NSGAII (population size: 100, evaluations: 24000)
pNSGAII.Max: pNSGAII (population size: 30, evaluations: 24000)

pSPEA2: pSPEA2 (population size: 30, evaluations: 12000)
SPEA2.1: SPEA2 (population size: 30, evaluations: 12000)
SPEA2.2: SPEA2 (population size: 100, evaluations: 12000)
SPEA2.3: SPEA2 (population size: 100, evaluations: 24000)
pSPEA2.Max: SPEA2 (population size: 30, evaluations: 24000)


HV.M: hypervolume, mean
HV.Md: hypervolume, median
HV.St: hypervolum, standard deviation
HV.IQR: hypervolume, interquartile range

EPS.M: Epsilon, mean
EPS.Md: Epsilon, median
EPS.St: Epsilon, standard deviation
EPS.IQR: Epsilon, interquartile range

IGD.M: Inverted generational distance, mean
IGD.Md: Inverted generational distance, median
IGD.St: Inverted generational distance, standard deviation
IGD.IQR: Inverted generational distance, interquartile range

Regions:
R0: {0.8, 0.95});
R1: {0.40, 0.50});
R2: {0.15, 0.20});
    

Other parameters (resultsZDT, resultsDTLZ):
Crossover: SBX
Mutation: polynomial mutation
crossover probability: 0.9
Mutation probability: 1/number of variables
Crossover distribution index: 20
Mutation distribution index: 20
Number of tournaments: 5 (for pNSGAII and pSPEA2)
epsilon grid: 0 (pAGE offline)
epsilon grid: 0.01 (pAGE online)
epsilon grid: 0.01 (AGE)
Ranking choice strategy parameter (for pNSGAII and pSPEA2): 10 
Individual elimination parameter (for pAGE): 10
probability form selecting individual from different region: 0.2 (for pNSGAII and pSPEA2)


Other parameters (resultsZDT(previous version), resultsDTLZ(previous version)):
Ranking choice strategy parameter (for pNSGAII and pSPEA2): 5 


DTLZ

pAGE.F.1: pAGE offline (population size: 30, evaluations: 50000)
pAGE.F.2: pAGE offline (population size: 40, evaluations: 50000)
pAGE.O.1: pAGE online (population size: 30, evaluations: 50000)
pAGE.O.2: pAGE online (population size: 40, evaluations: 50000)
AGE.1: AGE  (population size: 100, evaluations: 50000)
AGE.2: AGE (population size: 150, evaluations: 49950)

pNSGAII.1: pNSGAII (population size: 30, evaluations: 50000, 10 solutions are associated with each region) 
pNSGAII.2: pNSGAII (population size: 40, evaluations: 50000, 15 solutions are associated with region 0 and 1, 10 solutions are associated with region 2)
NSGAII.1: NSGAII (population size: 100, evaluations: 50000)
NSGAII.2: NSGAII (population size: 150, evaluations: 49950)

pSPEA2.1: pSPEA2 (population size: 30, evaluations: 50000, 10 solutions are associated with each region) 
pSPEA2.2: pSPEA2 (population size: 40, evaluations: 50000, 15 solutions are associated with region 0 and 1, 10 solutions are associated with region 2)
SPEA2.1: SPEA2 (population size: 100, evaluations: 50000)
SPEA2.2: SPEA2 (population size: 150, evaluations: 49950)
