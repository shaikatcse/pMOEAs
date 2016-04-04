//  NSGAIIRandom_Settings.java
//
//  Authors:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//
//  Copyright (c) 2013 Antonio J. Nebro
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
.metaheuristics.nsgaII.NSGAIIRandom;
import ares.fbk.eu.pMOEA.jmetal
.operators.selection.Selection;
import ares.fbk.eu.pMOEA.jmetal
.operators.selection.SelectionFactory;
import ares.fbk.eu.pMOEA.jmetal
.problems.ProblemFactory;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;
import ares.fbk.eu.pMOEA.jmetal
.util.offspring.DifferentialEvolutionOffspring;
import ares.fbk.eu.pMOEA.jmetal
.util.offspring.Offspring;
import ares.fbk.eu.pMOEA.jmetal
.util.offspring.PolynomialMutationOffspring;
import ares.fbk.eu.pMOEA.jmetal
.util.offspring.SBXCrossoverOffspring;

import java.util.HashMap;
import java.util.Properties;

/**
 * Settings class of algorithm NSGAIIRandom
 * Reference: Antonio J. Nebro, Juan José Durillo, Mirialys Machin Navas, Carlos A. Coello Coello, Bernabé Dorronsoro:
 * A Study of the Combination of Variation Operators in the NSGA-II Algorithm.
 * CAEPIA 2013: 269-278
 * DOI: http://dx.doi.org/10.1007/978-3-642-40643-0_28
 */
public class NSGAIIRandom_Settings extends Settings {
  public int populationSize_                 ;
  public int maxEvaluations_                 ;
  public double mutationProbability_         ;
  public double crossoverProbability_        ;
  public double mutationDistributionIndex_   ;
  public double crossoverDistributionIndex_  ;
  public double CR_                          ;
  public double F_                           ;

  /**
   * Constructor
   * @throws ares.fbk.eu.pMOEA.jmetal.util.JMException
   */
  public NSGAIIRandom_Settings(String problem) throws JMException {
    super(problem) ;
    
    Object [] problemParams = {"Real"};
    try {
	    problem_ = (new ProblemFactory()).getProblem(problemName_, problemParams);
    } catch (JMException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
    }  
    // Default settings
    populationSize_              = 100   ;
    maxEvaluations_              = 150000 ;
    mutationProbability_         = 1.0/problem_.getNumberOfVariables() ;
    crossoverProbability_        = 0.9 ;
    mutationDistributionIndex_   = 20 ;
    crossoverDistributionIndex_  = 20 ;
    CR_                          = 1.0 ;
    F_                           = 0.5 ;
  } // NSGAII_Settings

  
  /**
   * Configure NSGAII with user-defined parameter settings
   * @return A NSGAII algorithm object
   * @throws ares.fbk.eu.pMOEA.jmetal.util.JMException
   */
  public Algorithm configure() throws JMException {
    Algorithm algorithm ;
    Selection  selection ;

    HashMap  parameters ; // Operator parameters
    
    algorithm = new NSGAIIRandom(problem_) ;
    
    // Algorithm parameters
    algorithm.setInputParameter("populationSize",populationSize_);
    algorithm.setInputParameter("maxEvaluations",maxEvaluations_);

    Offspring[] getOffspring = new Offspring[3];
    double CR, F;
    getOffspring[0] = new DifferentialEvolutionOffspring(CR_, F_);

    getOffspring[1] = new SBXCrossoverOffspring(crossoverProbability_, crossoverDistributionIndex_);

    getOffspring[2] = new PolynomialMutationOffspring(mutationProbability_, mutationDistributionIndex_);

    algorithm.setInputParameter("offspringsCreators", getOffspring);

    // Selection Operator 
    parameters = null ;
    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;     

    // Add the operators to the algorithm
    algorithm.addOperator("selection",selection);
    
    return algorithm ;
  } // configure

  /**
   * Configure NSGAIIRandom with user-defined parameter experiments.settings
   * @return A NSGAIIRandom algorithm object
   */
  @Override
  public Algorithm configure(Properties configuration) throws JMException {
    Algorithm algorithm ;
    Selection  selection ;

    HashMap  parameters ; // Operator parameters

    // Creating the algorithm.
    algorithm = new NSGAIIRandom(problem_) ;

    // Algorithm parameters
    populationSize_ = Integer.parseInt(configuration.getProperty("populationSize",String.valueOf(populationSize_)));
    maxEvaluations_  = Integer.parseInt(configuration.getProperty("maxEvaluations",String.valueOf(maxEvaluations_)));
    algorithm.setInputParameter("populationSize",populationSize_);
    algorithm.setInputParameter("maxEvaluations",maxEvaluations_);

    // Mutation and Crossover for Real codification
    crossoverProbability_ = Double.parseDouble(configuration.getProperty("crossoverProbability",String.valueOf(crossoverProbability_)));
    crossoverDistributionIndex_ = Double.parseDouble(configuration.getProperty("crossoverDistributionIndex",String.valueOf(crossoverDistributionIndex_)));
    mutationProbability_ = Double.parseDouble(configuration.getProperty("mutationProbability",String.valueOf(mutationProbability_)));
    mutationDistributionIndex_ = Double.parseDouble(configuration.getProperty("mutationDistributionIndex",String.valueOf(mutationDistributionIndex_)));
    CR_ = Double.parseDouble(configuration.getProperty("CR",String.valueOf(CR_)));
    F_ = Double.parseDouble(configuration.getProperty("F",String.valueOf(F_)));

    Offspring[] getOffspring = new Offspring[3];
    getOffspring[0] = new DifferentialEvolutionOffspring(CR_, F_);

    getOffspring[1] = new SBXCrossoverOffspring(crossoverProbability_, crossoverDistributionIndex_);

    getOffspring[2] = new PolynomialMutationOffspring(mutationProbability_, mutationDistributionIndex_);

    algorithm.setInputParameter("offspringsCreators", getOffspring);

    // Selection Operator
    parameters = null ;
    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;

    // Add the operators to the algorithm
    algorithm.addOperator("selection",selection);

    return algorithm ;
  }
} // NSGAIIRandom_Settings
