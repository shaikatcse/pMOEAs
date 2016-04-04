//  SPEA2.java
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

package ares.fbk.eu.pMOEA.jmetal.metaheuristics.spea2;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import ares.fbk.eu.pMOEA.jmetal.core.*;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Epsilon;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.InvertedGenerationalDistance;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.pHypervolume;
import ares.fbk.eu.pMOEA.jmetal.util.Distance;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;
import ares.fbk.eu.pMOEA.jmetal.util.ModifiedSpea2Fitness;
import ares.fbk.eu.pMOEA.jmetal.util.Ranking;
import ares.fbk.eu.pMOEA.jmetal.util.Spea2Fitness;

/**
 * This class representing the SPEA2 algorithm
 */
public class pSPEA2 extends Algorithm {
	
	File fileHV[], fileIGD[], fileEpsilon[];
	FileWriter fwHV[],  fwIGD[],  fwEpsilon[];
	BufferedWriter bwHV[],  bwIGD[],  bwEpsilon[];
	
	boolean indicatorTracking = false;
	// declaration: regional true Pretofront
	List<double[]>[] rTrueFronts;

	/**
	 * Defines the number of tournaments for creating the mating pool
	 */
	public static final int TOURNAMENTS_ROUNDS = 1;

	/**
	 * Constructor. Create a new SPEA2 instance
	 * 
	 * @param problem
	 *            Problem to solve
	 */
	public pSPEA2(Problem problem) {
		super(problem);
	} // Spea2
	
	public pSPEA2(Problem problem, long seed, String folderName, int numberOfRegion) {
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

	int populationSize, archiveSize;
	int maxEvaluations;
	int evaluations;

	List<Double[]> regions;
	int numberOfSolutionsPerRegion[],currentRegionNumber;

	/**
	 * Runs of the Spea2 algorithm.
	 * 
	 * @return a <code>SolutionSet</code> that is a set of non dominated
	 *         solutions as a result of the algorithm execution
	 * @throws JMException
	 */
	public SolutionSet execute() throws JMException, ClassNotFoundException {
		Operator crossoverOperator, mutationOperator, selectionOperator;
		SolutionSet population, archive, offSpringSolutionSet;

		// Read the params
		populationSize = ((Integer) getInputParameter("populationSize"))
				.intValue();
		archiveSize = ((Integer) getInputParameter("archiveSize")).intValue();
		maxEvaluations = ((Integer) getInputParameter("maxEvaluations"))
				.intValue();

		regions = ((List<Double[]>) getInputParameter("regions"));

		numberOfSolutionsPerRegion = ((int []) getInputParameter("numberOfSolutionsPerRegion"));

		// Read the operators
		crossoverOperator = operators_.get("crossover");
		mutationOperator = operators_.get("mutation");
		selectionOperator = operators_.get("selection");

		// Initialize the variables
		population = new SolutionSet(populationSize);
		archive = new SolutionSet(archiveSize);
		evaluations = 0;

		// -> Create the initial solutionSet
		Solution newSolution;
		for (int i = 0; i < populationSize; i++) {
			newSolution = new Solution(problem_);
			problem_.evaluate(newSolution);
			problem_.evaluateConstraints(newSolution);
			evaluations++;
			population.add(newSolution);
		}

		doAssignRegion(population, regions);
	

		while (evaluations < maxEvaluations) {
			
			
			SolutionSet union = ((SolutionSet) population).union(archive);
		
			doAssignRegion(union, regions);
			
			
			archive.clear();

			for (int i = 0; i < regions.size(); i++) {
				currentRegionNumber = i;

				SolutionSet regionalSolutions = unionSolutionsOfSameRegion(union);
				SolutionSet globalPopulation = new SolutionSet(union);
				globalPopulation.removeSolutionSetSolutions(regionalSolutions);

				ModifiedSpea2Fitness spea2Fitness = new ModifiedSpea2Fitness(
						regionalSolutions);
				spea2Fitness.fitnessAssign(regionalSolutions, globalPopulation, evaluations, maxEvaluations);
				SolutionSet tempArchiev = spea2Fitness
						.environmentalSelection(numberOfSolutionsPerRegion[currentRegionNumber]);

				// strenghtRawFitness.computeDensityEstimator(regionalSolutions,
				// globalPopulation);
				// tempArchive =
				// environmentalSelection.execute(regionalSolutions);

				archive.addAll(tempArchiev);

			}
			
			if(indicatorTracking){
				trackIndicators(population);
			}

			// Spea2Fitness spea = new Spea2Fitness(union);
			// spea.fitnessAssign();
			// archive = spea.environmentalSelection(archiveSize);
			// Create a new offspringPopulation
			offSpringSolutionSet = new SolutionSet(populationSize);
			Solution[] parents = new Solution[2];

			for (int i = 0; i < regions.size(); i++) {
				selectionOperator.setParameter("currentRegion", i);
				currentRegionNumber = i;

				for (int z = 0; z < numberOfSolutionsPerRegion[i] / 2 *2; z++) {
			
				int j = 0;
				do {
					j++;
					parents[0] = (Solution) selectionOperator.execute(archive);
				} while (j < pSPEA2.TOURNAMENTS_ROUNDS); // do-while
				int k = 0;
				do {
					k++;
					parents[1] = (Solution) selectionOperator.execute(archive);
				} while (k < pSPEA2.TOURNAMENTS_ROUNDS); // do-while

				// make the crossover
				Solution[] offSpring = (Solution[]) crossoverOperator
						.execute(parents);
				mutationOperator.execute(offSpring[0]);
				problem_.evaluate(offSpring[0]);
				problem_.evaluateConstraints(offSpring[0]);
				offSpringSolutionSet.add(offSpring[0]);
				evaluations++;
			}
			}
			population = offSpringSolutionSet;
		
		
		} // while

		// Ranking ranking = new Ranking(archive);
		// return ranking.getSubfront(0);
		
		if(indicatorTracking){
			closeAllFiles();
		}
		
		return getResult(archive);

	} // execute

	public SolutionSet getResult(SolutionSet population) {
		doAssignRegion(population, regions);

		List<Solution> temp = new ArrayList<Solution>();

		for (int i = 0; i < population.size(); i++) {
			Solution sol = population.get(i);
			if (sol.getDistanceFromCloestRegion() == 0.0) {
				temp.add(sol);
			}
		}

		SolutionSet finalFieasiblepopulation = new SolutionSet(temp.size());
		for (int i = 0; i < temp.size(); i++) {
			finalFieasiblepopulation.add(temp.get(i));
		}

	/*	System.out.println("final");
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
	
} // SPEA2
