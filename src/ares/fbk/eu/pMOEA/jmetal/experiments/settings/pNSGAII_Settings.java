//  NSGAII_Settings.java 
//
//  Authors:
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
.experiments.settings;

import ares.fbk.eu.pMOEA.jmetal
.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal
.experiments.Settings;
import ares.fbk.eu.pMOEA.jmetal
.metaheuristics.nsgaII.NSGAII;
import ares.fbk.eu.pMOEA.jmetal.metaheuristics.nsgaII.pNSGAII;
import ares.fbk.eu.pMOEA.jmetal
.operators.crossover.Crossover;
import ares.fbk.eu.pMOEA.jmetal
.operators.crossover.CrossoverFactory;
import ares.fbk.eu.pMOEA.jmetal
.operators.mutation.Mutation;
import ares.fbk.eu.pMOEA.jmetal
.operators.mutation.MutationFactory;
import ares.fbk.eu.pMOEA.jmetal
.operators.selection.Selection;
import ares.fbk.eu.pMOEA.jmetal
.operators.selection.SelectionFactory;
import ares.fbk.eu.pMOEA.jmetal
.problems.ProblemFactory;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

/**
 * Settings class of algorithm NSGA-II (real encoding)
 */
public class pNSGAII_Settings extends Settings {
  public int populationSize_                 ;
  public int maxEvaluations_                 ;
  public double mutationProbability_         ;
  public double crossoverProbability_        ;
  public double mutationDistributionIndex_   ;
  public double crossoverDistributionIndex_  ;
  
  public int numberOfTournaments_;
  public int numberOfSolutionPerRegion_[];
  
  public List<Double[]> regions= new ArrayList<Double[]>(3);

  /**
   * Constructor
   */
  public pNSGAII_Settings(String problem) {
    super(problem) ;

    Object [] problemParams = {"Real"};
    try {
	    problem_ = (new ProblemFactory()).getProblem(problemName_, problemParams);
    } catch (JMException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
    }
    
    numberOfSolutionPerRegion_ = new int [] {10, 10, 10};
    
    //DTLZ
    //numberOfSolutionPerRegion_ = new int [] {15, 15, 10};
    
    // Default experiments.settings
    populationSize_              = 0   ;
    for(int i=0;i<numberOfSolutionPerRegion_.length; i++){
		populationSize_+=numberOfSolutionPerRegion_[i];
	}
    maxEvaluations_              = 12000 ;
    mutationProbability_         = 1.0/problem_.getNumberOfVariables() ;
    crossoverProbability_        = 0.9   ;
    mutationDistributionIndex_   = 20.0  ;
    crossoverDistributionIndex_  = 20.0  ;
    
    numberOfTournaments_ = 5;
    
    regions.add(new Double[] {0.8, 0.95});
    regions.add(new Double[] {0.40, 0.50});
    regions.add(new Double[] {0.15, 0.20});
    
    /*regions.add(new Double[] {0.8, 0.85});
    regions.add(new Double[] {0.35, 0.50});
    regions.add(new Double[] {0.15, 0.20});*/
    
    
  } // NSGAII_Settings


  /**
   * Configure NSGAII with default parameter experiments.settings
   * @return A NSGAII algorithm object
   * @throws ares.fbk.eu.pMOEA.jmetal.util.JMException
   */
  public Algorithm configure() throws JMException {
    Algorithm algorithm ;
    Selection  selection ;
    Crossover  crossover ;
    Mutation   mutation  ;

    HashMap  parameters ; // Operator parameters

	    
    
	algorithm = new pNSGAII(problem_);
	// algorithm = new ssNSGAII(problem);

	
	// Algorithm parameters
	
	algorithm.setInputParameter("populationSize", populationSize_);
	algorithm.setInputParameter("maxEvaluations", maxEvaluations_);
	algorithm.setInputParameter("numberOfSolutionsPerRegion", numberOfSolutionPerRegion_);
	algorithm.setInputParameter("regions", regions);
	
	// Mutation and Crossover for Real codification
	parameters = new HashMap();
	parameters.put("probability", crossoverProbability_);
	parameters.put("distributionIndex", crossoverDistributionIndex_);
	crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
			parameters);

	parameters = new HashMap();
	parameters.put("probability", mutationProbability_);
	parameters.put("distributionIndex", mutationDistributionIndex_);
	mutation = MutationFactory.getMutationOperator("PolynomialMutation",
			parameters);

	// Selection Operator
	parameters = new HashMap();
	parameters.put("numberOfTournaments", numberOfTournaments_);
	parameters.put("numberOfSolutionsPerRegion", numberOfSolutionPerRegion_);
	parameters.put("numberOfRegions", numberOfSolutionPerRegion_.length);
	selection = SelectionFactory.getSelectionOperator("PreferredTournamentSelection",
			parameters);
	

	// Add the operators to the algorithm
	algorithm.addOperator("crossover", crossover);
	algorithm.addOperator("mutation", mutation);
	algorithm.addOperator("selection", selection);


    return algorithm ;
  } // configure

 /**
  * Configure NSGAII with user-defined parameter experiments.settings
  * @return A NSGAII algorithm object
  */
  @Override
  public Algorithm configure(Properties configuration) throws JMException {
    Algorithm algorithm ;
    Selection  selection ;
    Crossover  crossover ;
    Mutation   mutation  ;

    HashMap  parameters ; // Operator parameters

    // Creating the algorithm. There are two choices: NSGAII and its steady-
    // state variant ssNSGAII
    algorithm = new NSGAII(problem_) ;
    //algorithm = new ssNSGAII(problem_) ;

    // Algorithm parameters
    populationSize_ = Integer.parseInt(configuration.getProperty("populationSize",String.valueOf(populationSize_)));
    maxEvaluations_  = Integer.parseInt(configuration.getProperty("maxEvaluations",String.valueOf(maxEvaluations_)));
    algorithm.setInputParameter("populationSize",populationSize_);
    algorithm.setInputParameter("maxEvaluations",maxEvaluations_);

    // Mutation and Crossover for Real codification
    crossoverProbability_ = Double.parseDouble(configuration.getProperty("crossoverProbability",String.valueOf(crossoverProbability_)));
    crossoverDistributionIndex_ = Double.parseDouble(configuration.getProperty("crossoverDistributionIndex",String.valueOf(crossoverDistributionIndex_)));
    parameters = new HashMap() ;
    parameters.put("probability", crossoverProbability_) ;
    parameters.put("distributionIndex", crossoverDistributionIndex_) ;
    crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);

    mutationProbability_ = Double.parseDouble(configuration.getProperty("mutationProbability",String.valueOf(mutationProbability_)));
    mutationDistributionIndex_ = Double.parseDouble(configuration.getProperty("mutationDistributionIndex",String.valueOf(mutationDistributionIndex_)));
    parameters = new HashMap() ;
    parameters.put("probability", mutationProbability_) ;
    parameters.put("distributionIndex", mutationDistributionIndex_) ;
    mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

    // Selection Operator
    parameters = null ;
    numberOfTournaments_ = Integer.parseInt(configuration.getProperty("numberOfTournaments",String.valueOf(numberOfTournaments_)));
   String numberOfSolutionPerRegion = configuration.getProperty("numberOfSolutionPerRegion",String.valueOf(numberOfSolutionPerRegion_));
    selection = SelectionFactory.getSelectionOperator("PreferredTournamentSelection", parameters) ;


    // Add the operators to the algorithm
    algorithm.addOperator("crossover",crossover);
    algorithm.addOperator("mutation",mutation);
    algorithm.addOperator("selection",selection);

    return algorithm ;
  }
} // NSGAII_Settings
