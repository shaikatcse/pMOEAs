//  NSGAII.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package ares.fbk.eu.pMOEA.jmetal.metaheuristics.nsgaII;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;

//import org.uma.jmetal.solution.Solution;
//import org.uma.jmetal.util.comparator.DominanceComparator;

import ares.fbk.eu.pMOEA.jmetal.util.ModifiedDominanceRanking;
import ares.fbk.eu.pMOEA.jmetal.core.*;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Epsilon;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.InvertedGenerationalDistance;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.pHypervolume;
import ares.fbk.eu.pMOEA.jmetal.util.Distance;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;
import ares.fbk.eu.pMOEA.jmetal.util.Ranking;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.CrowdingComparator;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator;

/**
 * Implementation of preference based NSGA-II.
 * 
 */

public class pNSGAII extends Algorithm {
	
	
	File fileHV[], fileIGD[], fileEpsilon[];
	FileWriter fwHV[],  fwIGD[],  fwEpsilon[];
	BufferedWriter bwHV[],  bwIGD[],  bwEpsilon[];
	
	boolean indicatorTracking = false;
	// declaration: regional true Pretofront
	List<double[]>[] rTrueFronts;
	
	/**
	 * Constructor
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public pNSGAII(Problem problem) {
		super(problem);
	} // NSGAII

	
	public pNSGAII(Problem problem, long seed, String folderName, int numberOfRegion) {
		super(problem);
		
		indicatorTracking = true;
		
		fileHV = new File[numberOfRegion];
		fileIGD = new File[numberOfRegion];
		fileEpsilon = new File[numberOfRegion];
		fwHV = new FileWriter[numberOfRegion];
		fwIGD = new FileWriter[numberOfRegion];
		fwEpsilon= new FileWriter[numberOfRegion];
		bwHV = new BufferedWriter[numberOfRegion];
		bwIGD = new BufferedWriter[numberOfRegion];
		bwEpsilon = new BufferedWriter[numberOfRegion];
		
		
				
		//repairSolution = new RepairSolution();
		if(!(new File(folderName+"\\HV").exists()))
			new File(folderName+"\\HV").mkdirs();
		
		if(!(new File(folderName+"\\IGD").exists()))
			new File(folderName+"\\IGD").mkdirs();
		
		if(!(new File(folderName+"\\Epsilon").exists()))
			new File(folderName+"\\Epsilon").mkdirs();
		
		
		for(int i=0;i<numberOfRegion;i++){
			
			fileHV[i] = new File(folderName + "\\HV\\trackHV_r"+i+"_" + seed);
		
			fileIGD[i] = new File(folderName + "\\IGD\\trackIGD_r"+i+"_"  + seed);
		
			fileEpsilon[i] = new File(folderName + "\\Epsilon\\trackEpsilon_r"+i+"_"  + seed);
		}
		
		// if file doesnt exists, then create it
		for(int i=0;i<numberOfRegion;i++){
		try {
			if (!fileHV[i].exists()) {
				fileHV[i].createNewFile();
			
				fileIGD[i].createNewFile();
				
				fileEpsilon[i].createNewFile();
				
			}

			fwHV[i] = new FileWriter(fileHV[i].getAbsoluteFile());
		
			fwIGD[i] = new FileWriter(fileIGD[i].getAbsoluteFile());
			
			fwEpsilon[i] = new FileWriter(fileEpsilon[i].getAbsoluteFile());
			

			bwHV[i] = new BufferedWriter(fwHV[i]);
			
			bwIGD[i] = new BufferedWriter(fwIGD[i]);
		
			bwEpsilon[i] = new BufferedWriter(fwEpsilon[i]);
			

		} catch (IOException e) {
			e.printStackTrace();
		}
		}
	}
	
	
	/**
	 * Runs the pNSGA-II algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */

	int populationSize;
	int maxEvaluations;
	int evaluations;

	Distance distance;

	List<Double[]> regions;
	int numberOfSolutionsPerRegion[], currentRegionNumber;

	public SolutionSet execute() throws JMException, ClassNotFoundException {

		QualityIndicator indicators; // QualityIndicator object
		int requiredEvaluations; // Use in the example of use of the
		// indicators object (see below)

		SolutionSet population;
		SolutionSet offspringPopulation;
		SolutionSet jointPopulation;

		Operator mutationOperator;
		Operator crossoverOperator;
		Operator selectionOperator;

		String paretoFrontFilePath;
		
		distance = new Distance();

		// Read the parameters
		populationSize = ((Integer) getInputParameter("populationSize"))
				.intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations"))
				.intValue();
		indicators = (QualityIndicator) getInputParameter("indicators");
		
				
		regions = ((List<Double[]>) getInputParameter("regions"));

		numberOfSolutionsPerRegion = ((int[]) getInputParameter("numberOfSolutionsPerRegion"));
		
			
		// Initialize the variables
		population = new SolutionSet(populationSize);
		evaluations = 0;

		// Read the operators
		mutationOperator = operators_.get("mutation");
		crossoverOperator = operators_.get("crossover");
		selectionOperator = operators_.get("selection");

		//for testing
		int acounterforPinting=0;
		
		// Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			evaluations++;
			population.add(newSolution);
		} // for

		//write population for each generation
	//	population.printObjectivesToFile("WithTrack\\pNSGAII\\FUN_"+ (int) evaluations / populationSize);
		
		doAssignRegion(population, regions);
		
		if(indicatorTracking){
			identifyRegionalTrueParetoFront();
			trackIndicators(population);
		}
		
		while(evaluations<maxEvaluations){

		// Generations
		offspringPopulation = new SolutionSet(populationSize);
		Solution[] parents = new Solution[2];

		for (int i = 0; i < regions.size(); i++) {
			selectionOperator.setParameter("currentRegion", i);
			currentRegionNumber = i;

			for (int j = 0; j < numberOfSolutionsPerRegion[i] / 2; j++) {
				
						// obtain parents
						parents[0] = (Solution) selectionOperator
								.execute(population);
						parents[1] = (Solution) selectionOperator
								.execute(population);
						Solution[] offSpring = (Solution[]) crossoverOperator
								.execute(parents);
						mutationOperator.execute(offSpring[0]);
						mutationOperator.execute(offSpring[1]);
						problem_.evaluate(offSpring[0]);
						problem_.evaluateConstraints(offSpring[0]);
						problem_.evaluate(offSpring[1]);
						problem_.evaluateConstraints(offSpring[1]);
						offspringPopulation.add(offSpring[0]);
						offspringPopulation.add(offSpring[1]);
						evaluations += 2;
					} // if
				} // for
		
			// Create the solutionSet union of solutionSet and offSpring
			jointPopulation = ((SolutionSet) population)
					.union(offspringPopulation);
			
			//jointPopulation.printObjectivesToFile("WithTrack\\pNSGAII\\"+ acounterforPinting++);

			doAssignRegion(jointPopulation, regions	 );
			
			/*
			 * temporary code
			 */
		/*	SolutionSet p1 = new SolutionSet(numberOfSolutionsPerRegion[0]*2);
			SolutionSet p2 = new SolutionSet(numberOfSolutionsPerRegion[1]*2);
			
			for(int zz =0; zz< jointPopulation.size();zz++){
				if(jointPopulation.get(zz).getAssignedRegionNumber() == 0)
					p1.add(jointPopulation.get(zz));
				else
					p2.add(jointPopulation.get(zz));
			}
			
			p1.printObjectivesToFile("assignRegion\\Region1_"+ (int) evaluations / populationSize);
			p2.printObjectivesToFile("assignRegion\\Region2_"+ (int) evaluations / populationSize);
			*/
			
			/*
			 * temporary code End
			 */
			

			population.clear();

			
			for (int i = 0; i < regions.size(); i++) {
				currentRegionNumber = i;

				SolutionSet regionalSolutions = unionSolutionsOfSameRegion(
						jointPopulation);
				SolutionSet globalPopulation = new SolutionSet(jointPopulation);
				globalPopulation.removeSolutionSetSolutions(regionalSolutions);

				ModifiedDominanceRanking ranking = computeRanking(
						regionalSolutions,globalPopulation, evaluations, maxEvaluations);

								
				/*SolutionSet tempPopulation = crowdingDistanceSelection(ranking);
				Comparator DOMINANCE_COMPARATOR = new DominanceComparator();
				int dominateMe[] = new int[tempPopulation.size()];
				for (int j = 0; j < tempPopulation.size() - 1; j++) {
					for (int k = j + 1; j < tempPopulation.size(); j++) {
						int flagDominate = DOMINANCE_COMPARATOR.compare(
								tempPopulation.get(j), tempPopulation.get(k));
						if (flagDominate == -1) {
							dominateMe[k]++;
						} else if (flagDominate == 1) {
							dominateMe[j]++;
						}

					}
				}

				for (int j = 0; j < tempPopulation.size(); j++)
					System.out.println(j + " " + dominateMe[j]);*/

				population.addAll(crowdingDistanceSelection(ranking));

			}
			
			//write population for each generation
			//population.printObjectivesToFile("WithTrack\\pNSGAII\\FUN_"+ (int) evaluations / populationSize);
			//population.printObjectivesToFile("WithTrack\\pNSGAII\\"+ acounterforPinting++);
					
			if(indicatorTracking){
				trackIndicators(population);
			}
		}
		

		// Return the first non-dominated front
		//ModifiedDominanceRanking ranking1 = new ModifiedDominanceRanking(population);
		//ranking1.getSubfront(0).printFeasibleFUN("FUN_NSGAII");
		
		if(indicatorTracking){
			closeAllFiles();
		}
		
		//testing
		//return population;

		// check: return only feasible solutions
		return getResult(population);
	}// execute

	public SolutionSet getResult(SolutionSet population) {
		doAssignRegion(population, regions);

		List<Solution> temp = new ArrayList<Solution>();

		for (int i = 0; i < population.size(); i++) {
			Solution sol = population.get(i);
			if (sol.getDistanceFromCloestRegion() == 0.0) {
			//if (sol.getDistanceFromCloestRegion() == 0.0 && sol.getDominateMeByGP() == 0 ) {
				temp.add(sol);
			}
		}

		SolutionSet finalFieasiblepopulation = new SolutionSet(temp.size());
		for (int i = 0; i < temp.size(); i++) {
			finalFieasiblepopulation.add(temp.get(i));
		}

		/*System.out.println("final");
		Comparator DOMINANCE_COMPARATOR = new ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator();
		int dominateMe[] = new int[finalFieasiblepopulation.size()];
		for (int i = 0; i < finalFieasiblepopulation.size() - 1; i++) {
			for (int j = i + 1; j < finalFieasiblepopulation.size(); j++) {
				int flagDominate = DOMINANCE_COMPARATOR.compare(
						finalFieasiblepopulation.get(i),
						finalFieasiblepopulation.get(j));
				if (flagDominate == -1) {
					dominateMe[j]++;
				} else if (flagDominate == 1) {
					dominateMe[i]++;
				}

			}
		}

		int numberOfDominantIndv = 0;
		for (int i = 0; i < finalFieasiblepopulation.size(); i++) {
			if (dominateMe[i] > 0)
				numberOfDominantIndv++;
		}

		System.out.println("Number of dominant individual: "
				+ numberOfDominantIndv);*/

		return finalFieasiblepopulation;

	}

	/*protected Ranking computeRanking(SolutionSet solutionList) {
		Ranking ranking = new Ranking(solutionList);
		return ranking;
	}*/

	protected ModifiedDominanceRanking computeRanking(
			SolutionSet localPolulation, int evaluations, int maxEvaluations) {
		ModifiedDominanceRanking ranking = new ModifiedDominanceRanking(
				localPolulation, evaluations, maxEvaluations);
		// ranking.setCurrentGeneration(iterations/populationSize);
		// ranking.computeRanking(localPolulation, globalPolulation);
		// ranking.computeRanking(localPolulation);

		return ranking;
	}
	
	protected ModifiedDominanceRanking computeRanking(
			SolutionSet localPolulation, SolutionSet globalPopulation,  int evaluations, int maxEvaluations) {
		ModifiedDominanceRanking ranking = new ModifiedDominanceRanking(
				localPolulation, globalPopulation, evaluations, maxEvaluations);
		// ranking.setCurrentGeneration(iterations/populationSize);
		// ranking.computeRanking(localPolulation, globalPolulation);
		// ranking.computeRanking(localPolulation);

		return ranking;
	}
	
	protected ModifiedDominanceRanking computeRanking(
			SolutionSet population) {
		ModifiedDominanceRanking ranking = new ModifiedDominanceRanking(
				population);
		// ranking.setCurrentGeneration(iterations/populationSize);
		// ranking.computeRanking(localPolulation, globalPolulation);
		// ranking.computeRanking(localPolulation);

		return ranking;
	}

	 protected SolutionSet crowdingDistanceSelection(
			ModifiedDominanceRanking ranking) {

		SolutionSet population = new SolutionSet(numberOfSolutionsPerRegion[currentRegionNumber]);
		int rankingIndex = 0;
		while (populationIsNotFull(population)) {
			if (subfrontFillsIntoThePopulation(ranking, rankingIndex,
					population)) {
				addRankedSolutionsToPopulation(ranking, rankingIndex,
						population);
				rankingIndex++;
			} else {
				// crowdingDistance.computeDensityEstimator(ranking.getSubfront(rankingIndex));
				distance.crowdingDistanceAssignment(
						ranking.getSubfront(rankingIndex),
						problem_.getNumberOfObjectives());
				addLastRankedSolutionsToPopulation(ranking, rankingIndex,
						population);
			}
		}

		return population;
	}
	
	
	
	protected void addLastRankedSolutionsToPopulation(
			ModifiedDominanceRanking ranking, int rank, SolutionSet population) {
		SolutionSet currentRankedFront = ranking.getSubfront(rank);

		currentRankedFront.sort(new CrowdingComparator());

		int i = 0;
		while (population.size() < numberOfSolutionsPerRegion[currentRegionNumber]) {
			population.add(currentRankedFront.get(i));
			i++;
		}
	}
	
	
	protected boolean populationIsNotFull(SolutionSet population) {
		return population.size() < numberOfSolutionsPerRegion[currentRegionNumber];
	}
	

	protected boolean subfrontFillsIntoThePopulation(
			ModifiedDominanceRanking ranking, int rank, SolutionSet population) {
		return ranking.getSubfront(rank).size() < (numberOfSolutionsPerRegion[currentRegionNumber] - population
				.size());
	}
	
	
	protected void addRankedSolutionsToPopulation(
			ModifiedDominanceRanking ranking, int rank, SolutionSet population) {
		SolutionSet front;

		front = ranking.getSubfront(rank);

		for (int i = 0; i < front.size(); i++) {
			Solution solution = front.get(i);
			population.add(solution);
		}
	}
	

	public void doAssignRegion(SolutionSet solutionSet, List<Double[]> regions	) {

		int numberOfSolutionsInEachResion [] = new int[regions.size()];
		if(solutionSet.size() > populationSize ){
			// when joint population is considered
			for(int i=0;i<numberOfSolutionsInEachResion.length;i++) {
				numberOfSolutionsInEachResion[i]= this.numberOfSolutionsPerRegion[i]* 2;
			}
		}else{
			for(int i=0;i<numberOfSolutionsInEachResion.length;i++) {
				numberOfSolutionsInEachResion[i]= this.numberOfSolutionsPerRegion[i];
			}
		}
		
		int numbreOfRegions = regions.size();
		int assignedNumberOfSolutions[] = new int[numbreOfRegions];

		for (int j = 0; j < solutionSet.size(); j++) {
			Solution sol = solutionSet.get(j);
			Boolean isAssigned = false;
			for (int i = 0; i < numbreOfRegions; i++) {
				if (sol.getObjective(0) >= regions.get(i)[0]
						&& sol.getObjective(0) <= regions.get(i)[1]) {

					if (assignedNumberOfSolutions[i] < numberOfSolutionsInEachResion[i]) {
						sol.setAssignedRegionNumber(i);
						sol.setDistanceFromCloestRegion(0.0);
						isAssigned = true;
						assignedNumberOfSolutions[i]++;
						break;
					}

				}
			}

			if (!isAssigned) {
				double distance = Double.POSITIVE_INFINITY;
				double dis = Double.POSITIVE_INFINITY;
				int assignedRegion = -1;

				for (int i = 0; i < numbreOfRegions; i++) {

					if (assignedNumberOfSolutions[i] < numberOfSolutionsInEachResion[i])
						dis = calculateDistanceFromARegion(sol, i, regions);
					if (dis < distance
							&& assignedNumberOfSolutions[i] < numberOfSolutionsInEachResion[i]) {
						distance = dis;
						assignedRegion = i;
					}
				}

				sol.setAssignedRegionNumber(assignedRegion);
				sol.setDistanceFromCloestRegion(distance);
				assignedNumberOfSolutions[assignedRegion]++;
			}
		}
	}

	double calculateDistanceFromARegion(Solution sol, int regionNumber,
			List<Double[]> regions) {
		double distance = Double.POSITIVE_INFINITY;

		for (int i = 0; i < 2; i++) {
			double dis = Math.abs(regions.get(regionNumber)[i]
					- sol.getObjective(0));
			if (dis < distance) {
				distance = dis;
			}
		}

		return distance;
	}

	protected SolutionSet unionSolutionsOfSameRegion(SolutionSet union) {
		SolutionSet resultSolutionSet = new SolutionSet(
				numberOfSolutionsPerRegion[currentRegionNumber] * 2);
		for (int i = 0; i < union.size(); i++) {
			Solution s = union.get(i);
			if (s.getAssignedRegionNumber() == currentRegionNumber) {
				resultSolutionSet.add(s);
			}
		}
		return resultSolutionSet;
	}
	
	void identifyRegionalTrueParetoFront(){
		String paretoFrontFilePath = (String) getInputParameter("paretoFrontFilePath");
		 rTrueFronts= new List[regions
					.size()];
		
		if(indicatorTracking && paretoFrontFilePath != null){
			double[][] trueFront = new Hypervolume().utils_
					.readFront(paretoFrontFilePath);

			// declaration: regional solution Pretofront
			List<double[]>[] rFronts = new List[regions.size()];

			for (int i = 0; i < rTrueFronts.length; i++) {
				rTrueFronts[i] = new ArrayList<double[]>();
				rFronts[i] = new ArrayList<double[]>();
			}

			// identify regional true Paretofront
			for (int i = 0; i < trueFront.length; i++) {
				for (int regionNumber = 0; regionNumber < regions
						.size(); regionNumber++) {
					if (trueFront[i][0] >= regions
							.get(regionNumber)[0] /* left side */
							&& trueFront[i][0] <= regions
									.get(regionNumber)[1]) /*
															 * right
															 * side
															 */{
						double[] sol = trueFront[i];
						rTrueFronts[regionNumber].add(sol);

					}
				}
			}
		}
	}
	
	void trackIndicators(SolutionSet population){

		if (indicatorTracking) {
			int genNo = (int) evaluations / populationSize;
			double HV=-1.0, IGD=-1.0, EPS=-1.0 ;
			
			// identify regional solution front
			// declaration: regional true Pretofront
			List<double[]>[] rFronts = new List[regions
								.size()];
			// declaration: regional true Pretofront
			List <Solution>[] rFrontsSol = new ArrayList[regions.size()];
			
			for(int i=0;i<regions.size();i++){
				rFrontsSol[i]=new ArrayList<Solution>();
			}
			
			for (int i = 0; i < population.size(); i++) {
				for (int regionNumber = 0; regionNumber < regions
						.size(); regionNumber++) {
					if (population.get(i).getObjective(0) >= regions
							.get(regionNumber)[0]
							&& population.get(i).getObjective(0) <= regions
									.get(regionNumber)[1]) {
						
						rFrontsSol[regionNumber].add(population.get(i));
							

					}
				}
			}
			
			
			for (int i = 0; i < regions.size(); i++) {
				
				//temporary solutionSet
				SolutionSet tempSolSet = new SolutionSet(rFrontsSol[i].size());
				for(int jj = 0 ; jj<rFrontsSol[i].size();jj++){
					tempSolSet.add(rFrontsSol[i].get(jj));
				}
				Ranking ranking = new Ranking(tempSolSet);
				
				
				if(rTrueFronts[i].size() != 0 && rFrontsSol[i].size() != 0){

					double [][] trueFront = new double[rTrueFronts[i].size()][];
					double [][] solutionFront = new double[ranking.getSubfront(0).size()][];

					for (int j = 0; j < rTrueFronts[i].size(); j++) {
						trueFront[j] = rTrueFronts[i].get(j);
					}

					for (int j = 0; j < ranking.getSubfront(0).size(); j++) {
						double [] sol = new double[ranking.getSubfront(0).get(j).getNumberOfObjectives()] ;
						for(int z=0;z<ranking.getSubfront(0).get(j).getNumberOfObjectives();z++){
							sol[z] = ranking.getSubfront(0).get(j).getObjective(z);
						}
						solutionFront[j] = sol;
					}
					
			
					if (solutionFront.length != 0) {

						pHypervolume HV_ind = new pHypervolume();
						// Hypervolume indicators = new
						// Hypervolume();
						HV = HV_ind.hypervolume(
								solutionFront, trueFront,
								trueFront[0].length,
								regions.get(i)[1]);
					
						InvertedGenerationalDistance IGD_ind = new InvertedGenerationalDistance();
						// double[][] solutionFront =
						// indicators.utils_
						// .readFront(solutionFrontFile);
						// double[][] trueFront =
						// indicators.utils_.readFront(paretoFrontPath);
						IGD = IGD_ind
								.invertedGenerationalDistance(
										solutionFront,
										trueFront,
										trueFront[0].length);
					
						
						Epsilon EPS_ind = new Epsilon();
						// double[][] solutionFront =
						// indicators.utils_
						// .readFront(solutionFrontFile);
						// double[][] trueFront =
						// indicators.utils_.readFront(paretoFrontPath);
						EPS = EPS_ind.epsilon(
								solutionFront, trueFront,
								trueFront[0].length);
					}
					
					
					
				}
				try {
					bwHV[i].write(genNo + " " + HV + "\n");
					bwIGD[i].write(genNo + " " + IGD + "\n");
					bwEpsilon[i].write(genNo + " " + EPS + "\n");
					
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			

			
		}
	}

	void closeAllFiles(){
		for(int i=0;i<regions.size();i++){
			try {
				bwHV[i].close();
		
			bwIGD[i].close();

			bwEpsilon[i].close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
}
	}
	
	

} // NSGA-II
