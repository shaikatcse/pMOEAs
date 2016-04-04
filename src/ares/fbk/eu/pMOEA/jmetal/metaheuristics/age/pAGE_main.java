// TO UPDATE: THIS IS NOW pAGE (just 1 line modified)




/**
 * AGE_main.java
 * This class represents the AGE algorithm (GECCO 2013 version)
 * 
 * "A Fast Approximation-Guided Evolutionary Multi-Objective Algorithm"
 * Markus Wagner and Frank Neumann
 * Genetic and Evolutionary Computation Conference 2013
 * http://cs.adelaide.edu.au/~markus/pub/2013gecco-age2.pdf
 * 
 * Creator: Markus Wagner (wagner@acrocon.com)
 * Feel free to contact me. There are many gems hidden in this code.
 *
 * @author Markus Wagner
 * @version 1.1 (GECCO 2013 version)
 */
package ares.fbk.eu.pMOEA.jmetal.metaheuristics.age;

import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;

import ares.fbk.eu.pMOEA.jmetal.core.*;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.FitnessComparator;
import ares.fbk.eu.pMOEA.jmetal.operators.crossover.*   ;
import ares.fbk.eu.pMOEA.jmetal.operators.mutation.*    ;
import ares.fbk.eu.pMOEA.jmetal.operators.selection.*   ;
import ares.fbk.eu.pMOEA.jmetal.problems.*                  ;
import ares.fbk.eu.pMOEA.jmetal.problems.DTLZ.*;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.*;
import ares.fbk.eu.pMOEA.jmetal.problems.WFG.*;
import ares.fbk.eu.pMOEA.jmetal.problems.LZ09.* ;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;

public class pAGE_main {

  /**
   * Usage example
   * 
   * VIA COMMANDLINE
   * ===============
   * 
   * /home/mwagner/scratch/age/age.sh ZDT1 ZDT1.pf 100 100000 BinaryTournament 0.01 zdtmu100
   * 
   * -->
   * 
   * function:                                      ZDT1
   * Pareto front file:                             ZDT1.pf
   * population size mu:                            100
   * evaluations:                                   100000
   * parent selection:                              BinaryTournament
   * epsilon grid for the approximating archive:    0.01            (note: try out different orders of magnitude)
   * output folder for the final population:        zdtmu100
   * 
   * VIA METHOD CALLS
   * ================
   * 
   * you need to go via AGE.java: 
   * - either execute a mini-study via AGE_main.java (see examples below)
   * - or execute a single run (see examples at the very bottom of AGE.java)
   * 
   * DTLZ 1/2/3/4
   * ============
   * - Several DTLZ fronts are available in originalParetoFronts, due to space
   *   reasons (several GBs), we do not provide fronts with >1000 points. However, 
   *   these can be easily generated using the uniform sampler in /ParetoFrontGenerator.java
   *   For our publications, we sample the front 1.000.000 times.
   * - DTLZ 1 has the fronts DTLZ1.*.pf
   * - DTLZ 2/3/4 have the same fronts DTLZ2.*.pf
   * 
   * Note:
   * the seventh parameter to AGE_main.main is for a powerful printer: 
   * - yes: deactivates the regular print of the hypervolume + epsilon approximation of the current population
   * - no: regular prints
   * Important: the regular prints invoke costly computations, and thus should not be run for time-restricted studies!
   * 
   */
  
 
    
  public static double median(double[] d) {
        Arrays.sort(d);
        int position = d.length/2;
        return d[position];
    }
    
  public static void main(String [] args) throws JMException, IOException, ClassNotFoundException {
      
      /*
       *  NOTE: the following switch can be used to run "mini studies" on a particular function
       */
      
      if (!true) {
          pAGE_main.mainORIG(args);
          
      } else {
          double epsilonSum = 0;
          int repetitions = 10;
          String[] setup = null;
          if (args.length==0) {
//              setup = new String[]{"LZ09_F5", "originalParetoFronts\\LZ09_F5.pf", "100", "true", "true", "0", "false","300000","BinaryTournament2","foo"};
                setup = new String[]{"ZDT1", "originalParetoFronts\\WFG4.3D.pf", "30", "true", "true", "0.005", "true","10000","BinaryTournament","foo"}; // BT3: dominance, crowding
//                setup = new String[]{"DTLZ1_2D", "originalParetoFronts\\DTLZ1.2D.1000.pf", "100", "true", "true", "0", "true","100000","BinaryTournament2","foo"}; // BT3: dominance, crowding
//                setup = new String[]{"DTLZ4_2D", "originalParetoFronts\\DTLZ2.2D.1000.pf", "100", "true", "true", "0", "true","15000","RandomSelection","foo"}; // BT3: dominance, crowding
          } else {
              setup = args;
          }
          for (int i = 1; i<=repetitions; i++) {
              System.out.println("miniStudy: repetition "+i);
              double temp = pAGE_main.mainORIG(setup);
            epsilonSum += temp;   //problem, Pareto front, pop size, doCrossover, doMutation
            System.out.println("miniStudy: repetition "+(i) + "/"+repetitions+" result epsilon="+temp + " (avg so far: "+ epsilonSum/(i+0d) +")");
            System.gc();
          }
          double epsilonAvg = epsilonSum/repetitions;
          System.out.println("\n\nminiStudy\n"+
                  Arrays.toString(setup) + 
                  " --> avg of "+repetitions+" runs is "+ epsilonAvg);
      }
  }
  
  public static double mainORIG(String [] args) throws JMException, IOException, ClassNotFoundException {
    Problem   problem   ;         // The problem to solve
    Algorithm algorithm ;         // The algorithm to use
    Operator  crossover ;         // Crossover operator
    Operator  mutation  ;         // Mutation operator
    Operator  selection ;         // Selection operator

    QualityIndicator indicators ; // Object to get quality indicators
    HashMap parameters; // Operator parameters

    indicators = null ;
    if (args.length == 1) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
    } // if
    else if (args.length == 2 || args.length >= 6) {
      Object [] params = {"Real"};
      problem = (new ProblemFactory()).getProblem(args[0],params);
      //indicators = new QualityIndicator(problem, args[1]) ;
    } // if
    else { // Default problem
//      problem = new Kursawe("Real", 3);
//      problem = new LinearFunction("Real");
      //problem = new Kursawe("BinaryReal", 3);
      //problem = new Water("Real");
      //problem = new ZDT1("ArrayReal", 100);
      //problem = new ConstrEx("Real");
//      problem = new Griewank("Real",1);
      problem = new DTLZ1("Real");
//      indicators = new QualityIndicator(problem, "originalParetoFronts\\DTLZ1.2D.pf") ;
    } // else

    algorithm = new pAGE(problem);

    // Algorithm parameters
    algorithm.setInputParameter("maxEvaluations",25000);

    if (args.length >= 5) {
        algorithm.setInputParameter("populationSize",Integer.parseInt(args[2]));
        algorithm.setInputParameter("doCrossover", Boolean.parseBoolean(args[3]));
        algorithm.setInputParameter("doMutation", Boolean.parseBoolean(args[4]));
        algorithm.setInputParameter("epsilonGridWidth", Double.parseDouble(args[5]));
        algorithm.setInputParameter("doOnMPICluster", Boolean.parseBoolean(args[6]));
        algorithm.setInputParameter("maxEvaluations",Integer.parseInt(args[7]));

    } else {
        algorithm.setInputParameter("populationSize",100);
        algorithm.setInputParameter("doMutation", true);
        algorithm.setInputParameter("doCrossover", true);
    }
    algorithm.setInputParameter("infoPrinterHowOften", 100);
    algorithm.setInputParameter("infoPrinterSubDir", args[9]);

    // Mutation and Crossover for Real codification
    parameters = new HashMap();
	parameters.put("probability", 0.9);
	parameters.put("distributionIndex", 20.0);
	crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
			parameters);
	
	parameters = new HashMap();
	parameters.put("probability", 1.0 / problem.getNumberOfVariables());
	parameters.put("distributionIndex", 20.0);
	mutation = MutationFactory.getMutationOperator("PolynomialMutation",
			parameters);
    
    /* Mutation and Crossover Binary codification */

    /* Selection Operator */
//    selection = new RandomSelection();
  //  selection = new BinaryTournament(new FitnessComparator());                  // minimizer!
    
    parameters = new HashMap();
    
    parameters.put("comparator" ,  new FitnessComparator());
	selection = SelectionFactory.getSelectionOperator("BinaryTournament",
			parameters);
//    selection = new BinaryTournament(new FitnessAndCrowdingDistanceComparator());
    // Add the operators to the algorithm
    algorithm.addOperator("crossover",crossover);
    algorithm.addOperator("mutation",mutation);
    algorithm.addOperator("selection",selection);
    algorithm.setInputParameter("indicators", indicators) ;
    algorithm.addOperator("selection", (Selection) SelectionFactory.getSelectionOperator(args[8], parameters));

    // Execute the Algorithm
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    long estimatedTime = System.currentTimeMillis() - initTime;

    
	System.out.println("Total execution time: " + estimatedTime + "ms");
	System.out.println("Variables values have been writen to file VAR");
	population.printVariablesToFile("VAR");
	System.out.println("Objectives values have been writen to file FUN");
	population.printObjectivesToFile("FUN");
    
      
    double epsilon = 0;
    if (indicators != null && Boolean.parseBoolean(args[6])==false) {
      System.out.println("Quality indicators") ;
      System.out.println("Hypervolume: " + indicators.getHypervolume(population)) ;
      System.out.println("GD         : " + indicators.getGD(population)) ;
      System.out.println("IGD        : " + indicators.getIGD(population)) ;
      System.out.println("Spread     : " + indicators.getSpread(population)) ;
      epsilon = indicators.getEpsilon(population);
      System.out.println("Epsilon    : " + epsilon) ;  
    } // if
    //return indicators.getEpsilon(population);
    return 0.0;
  } //main
} 