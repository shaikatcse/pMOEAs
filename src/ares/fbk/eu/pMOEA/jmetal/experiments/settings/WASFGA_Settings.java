package ares.fbk.eu.pMOEA.jmetal
.experiments.settings;

import java.util.HashMap;

import ares.fbk.eu.pMOEA.jmetal
.metaheuristics.wasfga.WASFGA;
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
.core.*;
import ares.fbk.eu.pMOEA.jmetal
.experiments.Settings;
import ares.fbk.eu.pMOEA.jmetal
.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal
.util.AchievementScalarizingFunction;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;
import ares.fbk.eu.pMOEA.jmetal
.util.ReferencePoint;

/**
 * @author Rubén Saborido Infantes
 * Settings class of algorithm WASF-GA (real encoding)
 */
public class WASFGA_Settings extends Settings {
  public int populationSize_                 ; 
  public int maxEvaluations_                 ;
  public double mutationProbability_         ;
  public double crossoverProbability_        ;
  public double mutationDistributionIndex_   ;
  public double crossoverDistributionIndex_  ;
  QualityIndicator convergenceIndicator_     ;   
  public boolean normalization_              ;
  public String weightsDirectory_            ;  
  public AchievementScalarizingFunction asf_ ;
  public ReferencePoint referencePoint_		 ;
  public boolean estimatePoints_             ; 
  public String folderForOutputFiles_		 ;
  
  /**
   * Constructor
   * @throws JMException 
   */
  public WASFGA_Settings(String problem) throws JMException {
    super(problem) ;
    String paretoFrontDirectory = new String("data/paretoFronts/");
    
    Object [] problemParams = {"Real"};
    problem_ = (new ProblemFactory()).getProblem(problemName_, problemParams);        
    
    // Default settings
    populationSize_              = 300   ;
    maxEvaluations_              = 150000 ;
    mutationProbability_         = 1.0/problem_.getNumberOfVariables() ;
    crossoverProbability_        = 0.9   ;
    mutationDistributionIndex_   = 20.0  ;
    crossoverDistributionIndex_  = 20.0  ;        
    weightsDirectory_            = new String("data/weights");
    normalization_               = true;
    estimatePoints_              = true;        
    
    if (problem_.getName().substring(0, problem_.getName().length()-1).equals("ZDT") || problem_.getName().substring(0, problem_.getName().length()-1).equals("LZ09_F"))
        convergenceIndicator_        = new QualityIndicator(problem_, paretoFrontDirectory + problem_.getName() + ".pf");
    else
        convergenceIndicator_        = new QualityIndicator(problem_, paretoFrontDirectory + problem_.getName() + "." + problem_.getNumberOfObjectives() + "D.pf");       
  } // WASFGA_Settings

  
  /**
   * Configure WASFGA with user-defined parameter settings
   * @return A WASFGA algorithm object
   * @throws ares.fbk.eu.pMOEA.jmetal.util.JMException
   */
  public Algorithm configure() throws JMException {
    Algorithm algorithm ;
    Selection  selection ;
    Crossover  crossover ;
    Mutation   mutation  ;

    HashMap  parameters ; // Operator parameters  
            
    // Creating the algorithm.
    algorithm = new WASFGA(problem_) ;    
    
    // Algorithm parameters
    algorithm.setInputParameter("populationSize",populationSize_);
    algorithm.setInputParameter("maxEvaluations",maxEvaluations_);
    algorithm.setInputParameter("weightsDirectory",weightsDirectory_);
    algorithm.setInputParameter("indicators",convergenceIndicator_);        
    algorithm.setInputParameter("normalization",normalization_);
    algorithm.setInputParameter("asf",asf_);
    algorithm.setInputParameter("estimatePoints",estimatePoints_);
    algorithm.setInputParameter("referencePoint",referencePoint_);
    algorithm.setInputParameter("folderForOutputFiles", folderForOutputFiles_);    

    // Mutation and Crossover for Real codification
    parameters = new HashMap() ;
    parameters.put("probability", crossoverProbability_) ;
    parameters.put("distributionIndex", crossoverDistributionIndex_) ;
    crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);                   

    parameters = new HashMap() ;
    parameters.put("probability", mutationProbability_) ;
    parameters.put("distributionIndex", mutationDistributionIndex_) ;
    mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                        

    // Selection Operator 
    parameters = null ;
    selection = SelectionFactory.getSelectionOperator("BinaryTournament2", parameters) ;     

    // Add the operators to the algorithm
    algorithm.addOperator("crossover",crossover);
    algorithm.addOperator("mutation",mutation);
    algorithm.addOperator("selection",selection);
   
    return algorithm ;
  } // configure
} // WASFGA_Settings