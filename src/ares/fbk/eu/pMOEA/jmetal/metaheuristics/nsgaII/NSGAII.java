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

package ares.fbk.eu.pMOEA.jmetal
.metaheuristics.nsgaII;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ares.fbk.eu.pMOEA.jmetal
.core.*;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Epsilon;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.InvertedGenerationalDistance;
import ares.fbk.eu.pMOEA.jmetal
.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.pHypervolume;
import ares.fbk.eu.pMOEA.jmetal
.util.Distance;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;
import ares.fbk.eu.pMOEA.jmetal
.util.Ranking;
import ares.fbk.eu.pMOEA.jmetal
.util.comparators.CrowdingComparator;

/** 
 *  Implementation of NSGA-II.
 *  This implementation of NSGA-II makes use of a QualityIndicator object
 *  to obtained the convergence speed of the algorithm. This version is used
 *  in the paper:
 *     A.J. Nebro, J.J. Durillo, C.A. Coello Coello, F. Luna, E. Alba 
 *     "A Study of Convergence Speed in Multi-Objective Metaheuristics." 
 *     To be presented in: PPSN'08. Dortmund. September 2008.
 */

public class NSGAII extends Algorithm {
	
	File fileHV[], fileIGD[], fileEpsilon[];
	FileWriter fwHV[],  fwIGD[],  fwEpsilon[];
	BufferedWriter bwHV[],  bwIGD[],  bwEpsilon[];
	
	boolean indicatorTracking = false;
	// declaration: regional true Pretofront
	List<double[]>[] rTrueFronts;
	
  /**
   * Constructor
   * @param problem Problem to solve
   */
  public NSGAII(Problem problem) {
    super (problem) ;
  } // NSGAII
  
  
  public NSGAII(Problem problem, long seed, String folderName, int numberOfRegion) {
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
   * Runs the NSGA-II algorithm.
   * @return a <code>SolutionSet</code> that is a set of non dominated solutions
   * as a result of the algorithm execution
   * @throws JMException 
   */
  
  int populationSize;
  int maxEvaluations;
  int evaluations;
  
  List<Double[]> regions;
  
  public SolutionSet execute() throws JMException, ClassNotFoundException {


    QualityIndicator indicators; // QualityIndicator object
    int requiredEvaluations; // Use in the example of use of the
    // indicators object (see below)

    SolutionSet population;
    SolutionSet offspringPopulation;
    SolutionSet union;

    Operator mutationOperator;
    Operator crossoverOperator;
    Operator selectionOperator;

    Distance distance = new Distance();

    //Read the parameters
    populationSize = ((Integer) getInputParameter("populationSize")).intValue();
    maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();
    indicators = (QualityIndicator) getInputParameter("indicators");

    //Initialize the variables
    population = new SolutionSet(populationSize);
    evaluations = 0;

    requiredEvaluations = 0;

    //Read the operators
    mutationOperator = operators_.get("mutation");
    crossoverOperator = operators_.get("crossover");
    selectionOperator = operators_.get("selection");

    // Create the initial solutionSet
    Solution newSolution;
    for (int i = 0; i < populationSize; i++) {
      newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      evaluations++;
      population.add(newSolution);
    } //for   
    
  //write population for each generation
  	population.printObjectivesToFile("WithTrack\\NSGAII\\FUN_"+ (int) evaluations / populationSize);
  		
  	//for testing
  	int acounterforPinting=0;
    
    if(indicatorTracking){
    	regions = ((List<Double[]>) getInputParameter("regions"));
		identifyRegionalTrueParetoFront();
		trackIndicators(population);
	}

    // Generations 
    while (evaluations < maxEvaluations) {

      // Create the offSpring solutionSet      
      offspringPopulation = new SolutionSet(populationSize);
      Solution[] parents = new Solution[2];
      for (int i = 0; i < (populationSize / 2); i++) {
        if (evaluations < maxEvaluations) {
          //obtain parents
          parents[0] = (Solution) selectionOperator.execute(population);
          parents[1] = (Solution) selectionOperator.execute(population);
          Solution[] offSpring = (Solution[]) crossoverOperator.execute(parents);
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
      union = ((SolutionSet) population).union(offspringPopulation);
      
      //testing:
      //union.printObjectivesToFile("WithTrack\\NSGAII\\"+ acounterforPinting++);

      // Ranking the union
      Ranking ranking = new Ranking(union);

      int remain = populationSize;
      int index = 0;
      SolutionSet front = null;
      population.clear();

      // Obtain the next front
      front = ranking.getSubfront(index);

      while ((remain > 0) && (remain >= front.size())) {
        //Assign crowding distance to individuals
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        //Add the individuals of this front
        for (int k = 0; k < front.size(); k++) {
          population.add(front.get(k));
        } // for

        //Decrement remain
        remain = remain - front.size();

        //Obtain the next front
        index++;
        if (remain > 0) {
          front = ranking.getSubfront(index);
        } // if        
      } // while

      // Remain is less than front(index).size, insert only the best one
      if (remain > 0) {  // front contains individuals to insert                        
        distance.crowdingDistanceAssignment(front, problem_.getNumberOfObjectives());
        front.sort(new CrowdingComparator());
        for (int k = 0; k < remain; k++) {
          population.add(front.get(k));
        } // for

        remain = 0;
      } // if                               

      // This piece of code shows how to use the indicator object into the code
      // of NSGA-II. In particular, it finds the number of evaluations required
      // by the algorithm to obtain a Pareto front with a hypervolume higher
      // than the hypervolume of the true Pareto front.
      if ((indicators != null) &&
          (requiredEvaluations == 0)) {
        double HV = indicators.getHypervolume(population);
        if (HV >= (0.98 * indicators.getTrueParetoFrontHypervolume())) {
          requiredEvaluations = evaluations;
        } // if
      } // if
    
      
  		//write population for each generation
		population.printObjectivesToFile("WithTrack\\NSGAII\\FUN_"+ (int) evaluations / populationSize);
		//population.printObjectivesToFile("WithTrack\\NSGAII\\"+ acounterforPinting++);
		
		
      if(indicatorTracking){
			trackIndicators(population);
		}
    } // while

    // Return as output parameter the required evaluations
    setOutputParameter("evaluations", requiredEvaluations);

    // Return the first non-dominated front
    Ranking ranking = new Ranking(population);
    ranking.getSubfront(0).printFeasibleFUN("FUN_NSGAII") ;
    
    if(indicatorTracking){
		closeAllFiles();
	}

    return ranking.getSubfront(0);
  } // execute
  
  
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
