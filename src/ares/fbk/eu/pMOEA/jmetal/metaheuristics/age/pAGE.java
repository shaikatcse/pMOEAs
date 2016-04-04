/**
 * AGE.java 
 *
 * @author Markus Wagner
 * @version 1.1p (GECCO 2013 version + preferred regions)
 */
package ares.fbk.eu.pMOEA.jmetal.metaheuristics.age;

import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import ares.fbk.eu.pMOEA.jmetal.core.*;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator;
import ares.fbk.eu.pMOEA.jmetal.metaheuristics.age.util.DeepCopy;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal.util.*;

/**
 * This class represents the pAGE algorithm (GECCO 2013 AGE version).
 * The p stands for the variant that considers preferred regions.
 */
public class pAGE extends Algorithm {

    // pAGE START
    public static double[] preferredRanges = new double[]{0.15,0.2,  0.4,0.5,  0.8,0.95};
    public static int preferredObjective = 0; // 0 is 1, 1 is 2, ...
//    public static double preferredMin = preferredRanges[0];
    // pAGE END
    
    public static boolean isSolutionInPreferred(Solution x) {
        boolean debugPrint = true;
        boolean result = false;
        
        //double[] objectives = x.getObjectives();  // currently this code works for only one preferred objective
        double[] objectives = new double[x.getNumberOfObjectives()];
        for(int i=0;i<objectives.length;i++){
        	objectives[i] = x.getObjective(i);
        }
        
        if (debugPrint) System.out.println(Arrays.toString(objectives));
        double objective = objectives[preferredObjective];
        for (int i=0; i<preferredRanges.length; i=i+2) {
            if (preferredRanges[i]<=objective && objective<=preferredRanges[i+1]) {
                if (debugPrint) System.out.println("true");
                return true; // we can shortcut, since it is sufficient for a solution to be in one region
            }
        }
        if (debugPrint) System.out.println("false objective="+objective);
        return result;
    }
    
    static boolean debugPrintGlobal = !true;
    boolean debugPrintAGE = !true;
    /**
     * Defines the number of tournaments for creating the mating pool
     */
    public static final int TOURNAMENTS_ROUNDS = 1;
    /**
     * Stores the problem to solve
     */
   // private Problem problem_;

    /**
     * Constructor.
     * Create a new AGE instance
     * @param problem Problem to solve
     */
    public pAGE(Problem problem) {
        super (problem);
    }

    
    /**
     * Runs the AGE algorithm.
     * @return a <code>SolutionSet</code> that is a set of non dominated solutions
     * as a result of the algorithm execution
     * @throws JMException
     */
    public SolutionSet execute() throws JMException, ClassNotFoundException {

        QualityIndicator indicators = (QualityIndicator) getInputParameter("indicators"); // this line had to be added (mw)
        boolean doMutation = ((Boolean)getInputParameter("doMutation")).booleanValue();
        boolean doCrossover = ((Boolean)getInputParameter("doCrossover")).booleanValue();
       
        /*int infoPrinterHowOften;
        if (getInputParameter("infoPrinterHowOften")==null) infoPrinterHowOften=1000;
            else infoPrinterHowOften = ((Integer)getInputParameter("infoPrinterHowOften")).intValue();
        String infoPrinterSubDir = (String)getInputParameter("infoPrinterSubDir");*/
        double epsilonGridWidth = ((Double) getInputParameter("epsilonGridWidth")).doubleValue();
       /* boolean doOnMPICluster;
        if (getInputParameter("doOnMPICluster")==null) doOnMPICluster = false;
            else doOnMPICluster = ((Boolean) getInputParameter("doOnMPICluster")).booleanValue();*/

        int populationSize, maxEvaluations, evaluations;
        Operator crossoverOperator, mutationOperator, selectionOperator;
        SolutionSet solutionSet, archive, offSpringSolutionSet;

        //Read the params
        populationSize = ((Integer) getInputParameter("populationSize")).intValue();
        maxEvaluations = ((Integer) getInputParameter("maxEvaluations")).intValue();

        //Read the operators
        crossoverOperator = operators_.get("crossover");
        mutationOperator = operators_.get("mutation");
        selectionOperator = operators_.get("selection");
        System.out.println("selector:"+selectionOperator.toString());

        //Initialize the variables
        solutionSet = new SolutionSet(populationSize);
        evaluations = 0;

        //-> Create the initial solutionSet
        Solution newSolution;
        for (int i = 0; i < populationSize; i++) {
            newSolution = new Solution(problem_);
            problem_.evaluate(newSolution);
            problem_.evaluateConstraints(newSolution);
            evaluations++;
            solutionSet.add(newSolution);
        }
        
        /* Initialize the archive with the solutionSet. In subsequent iterations, 
         * the newly constructed point will first be added to the archive, and then
         * the best mu points out of mu+1 are selected that approximate the new 
         * archive best.
         */
        archive = new SolutionSet(solutionSet.size());
        for(int i=0; i<solutionSet.size(); i++) {
            archive.add( new Solution(solutionSet.get(i)) );
        }
        
        /* or: initialise with epsilonboxes, where a point is in the box */
        boolean useEpsilonBoxesArchive;
        
        if (epsilonGridWidth==0) useEpsilonBoxesArchive = false;
            else useEpsilonBoxesArchive = true;
        
        System.out.println("useEpsilonBoxesArchive="+useEpsilonBoxesArchive+" epsilonGridWidth="+epsilonGridWidth);
        
        if (useEpsilonBoxesArchive) {
            archive = new SolutionSet(populationSize);
            for (int i = 0; i<populationSize; i++) {
                Solution converted = convertSolutionToEpsilonGridVectorFLOOR(solutionSet.get(i),epsilonGridWidth);
                archive.add(converted); 
           }
        }
        
        System.out.println("initial: population.size()="+solutionSet.size()+" archive.size()="+archive.size());

        Comparator cNormal = new DominanceComparator();

        
//        int newPointIsDominatedByOldArchiveTakeNeverthelessCounter = 0;
//        int newPointIsDominatedByOldArchiveCounter = 0;
        
        // main loop starts...
        while (evaluations <= maxEvaluations) {
            
            /* START debug printouts */
         /*   if (infoPrinter==null) 
                if (doOnMPICluster) {
                    infoPrinter = new InfoPrinter(this, problem_, infoPrinterSubDir); 
                } else {
                    infoPrinter = new InfoPrinter(this, problem_, infoPrinterSubDir, infoPrinterHowOften, infoPrinterHowOften);
                }
            
//            if (infoPrinter==null) infoPrinter = new InfoPrinter(this, problem_, infoPrinterSubDir);

            if (evaluations%1000==0)
                if (doOnMPICluster) {
                    infoPrinter.printLotsOfValuesToFile(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, false);
                    infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, false);
                } else {
                    infoPrinter.printLotsOfValuesToFile(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, false, false);
                    infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, false, false);
                }
            if (evaluations>=maxEvaluations) {
                if (doOnMPICluster) {
                    infoPrinter.printLotsOfValuesToFile(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, true);
                    infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, true); //correct line
//                    infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, false, false);
                }
                break;
            } */
            /* END debug printouts */

            /* AGE: Create a new offspringPopulation */

            /* START AGE block */


            /* AGE 1. step: generate one new solution */
            offSpringSolutionSet = new SolutionSet(populationSize);                          // generate mu solutions
//            SolutionSet offSpringSolutionSetForArchive = new SolutionSet(populationSize);    // generate mu solutions
            Solution[] parents = new Solution[2];
            Solution[] offSpring = null;
            
            
            
            /* from CEC 2012 paper: thinning and crowding */
          
            // 1. take first front(s) and extreme points
            boolean debugSelection = !true;

            if (debugSelection) System.out.print("orig:"+solutionSet.size()+ " ");

            Distance distance = new Distance();

            Ranking tempr = new Ranking(solutionSet);
            int crowdingDistanceSwitch = -1;                                // for -1 for AGE-specific variants

            SolutionSet t = new SolutionSet(populationSize);

            if (debugSelection) System.out.print(" fronts:"+tempr.getNumberOfSubfronts()+ " ");
            t = t.union(tempr.getSubfront(0));                              // definitely take first front
            int included = 1;
            if (debugSelection) System.out.print(" 1st:"+t.size()+ " ");

            for (int i = included; i<tempr.getNumberOfSubfronts(); i++) {   // from next fronts take probabilistically
                SolutionSet x = tempr.getSubfront(i);
                distance.crowdingDistanceAssignment(x, problem_.getNumberOfObjectives(), crowdingDistanceSwitch, archive);

                for (int j = 0; j<x.size(); j++) {
                    if (Math.random()< (1d/(i+1)) ) {                 // decreasing probability depending on the rank number
                        t.add(x.get(j));
                        if (debugSelection) System.out.print("+");
                    } else { if (debugSelection) System.out.print("/");}

                }
                if (debugSelection) System.out.print(",");
//                    System.out.println(" new:"+t.size());
            }
            solutionSet = t;
            if (debugSelection) System.out.println(" new:"+solutionSet.size());
            
            // 2. crowding distances for each front
//            Distance distance = new Distance();
            tempr = new Ranking(solutionSet);
            for (int i = 0; i<tempr.getNumberOfSubfronts(); i++) {
                distance.crowdingDistanceAssignment(tempr.getSubfront(i), problem_.getNumberOfObjectives(), 0, archive);    // careful: use correct switch!
            }

            
            
            
            
//  /*mu+1*/      for (int kk = 0; kk<1 ; kk++){
   /*mu+mu*/    for (int kk = 0; kk<populationSize; kk++){                     // loop condition: generate lambda inividuals
                        parents[0] = (Solution) selectionOperator.execute(solutionSet); // carefull! the operator may work on the fitness values (which we did not really have in the beginning)
                        parents[1] = (Solution) selectionOperator.execute(solutionSet);

                    //make the crossover and generate a single child
                    if (doCrossover) offSpring = (Solution [])crossoverOperator.execute(parents);    // 2 parents are XOed
                        else offSpring = parents;                                                    // no XO
                    if (doMutation) mutationOperator.execute(offSpring[0]);                          // mutation

                    // FITNESS EVALUATION - note: this does not set fitness, just runs the problem functions
                    problem_.evaluate(offSpring[0]);
                    problem_.evaluateConstraints(offSpring[0]);
                    evaluations++;

                    

                    
                    /* START check if new offSpring is not (epsilon) dominated by an archive point */
                    boolean newPointIsDominatedByOldArchive = false;
                    boolean newPointIsDominatedByOldArchiveTakeNevertheless = false;
                    
                    boolean debugPrintEpsilon = false;
                    
                    for (int i = 0; i<archive.size(); i++) {
                        Solution archiveSolution = archive.get(i);
                        /*
                         * COMPARE: return -1, or 0, or 1 if solution1 dominates solution2, both are
                         *          non-dominated, or solution1 is dominated by solution2, respectively.
                         */                        
                        int result = cNormal.compare(archiveSolution, offSpring[0]);
                        
                        
                        if (debugPrintEpsilon) System.out.println(archiveSolution.toString() + " " + offSpring[0].toString() + " -> "+result);
                        
                            
                        int offspringVsArchiveByAtMostEpsilon = Integer.MAX_VALUE;
                        if (useEpsilonBoxesArchive) {
                            // floor offspring to the "lower left corner"
                            result = cNormal.compare(archiveSolution, convertSolutionToEpsilonGridVectorFLOOR(offSpring[0], epsilonGridWidth));
                            Solution archivePointMovedAway = moveEpsilonGridVectorOnceFurtherAway(archiveSolution, epsilonGridWidth);
                            offspringVsArchiveByAtMostEpsilon = cNormal.compare(archivePointMovedAway, offSpring[0]);
                            
                            
                            if (debugPrintEpsilon) System.out.println(" " +archiveSolution.toString() + " " + convertSolutionToEpsilonGridVectorFLOOR(offSpring[0], epsilonGridWidth) + " -> "+result);
                            if (debugPrintEpsilon) System.out.println(" " +archivePointMovedAway + " " + offSpring[0] + " -> "+offspringVsArchiveByAtMostEpsilon);
                        }
                        
                        
                        if (result==-1 || result==2) {
                            // break if an archive point dominates the new point
                            newPointIsDominatedByOldArchive = true;
//                            newPointIsDominatedByOldArchiveCounter++;
//                            System.out.println("dominated by i="+i+" of "+archive.size());
                            
                            if (useEpsilonBoxesArchive 
                                     && offspringVsArchiveByAtMostEpsilon==1
                                    ) {
                                newPointIsDominatedByOldArchiveTakeNevertheless = true;
                                if (debugPrintEpsilon) System.out.println("newPointIsDominatedByOldArchiveTakeNevertheless");
//                                System.out.println("newPointIsDominatedByOldArchiveTakeNevertheless");
                            }
                            
                            break;
                        }
                        if (result==1) {
                            // remove archive point if new point dominates that one
//                            System.out.println("remove i="+i+" of "+archive.size());
                            archive.remove(i);
                            i--;
                        }
                    }
                    /* END check if new offSpring is not (epsilon) dominated by an archive point */

                    
                    
                    /* the following can be done unconditionally, as we would have 
                       - stopped if offSpring[0] is dominated by the current population
                       - 
                     */
                    if (newPointIsDominatedByOldArchive) {
                        // in case of !useEpsilonBoxesArchive: forget this point
                        
                        if (/*useEpsilonBoxesArchive &&*/ newPointIsDominatedByOldArchiveTakeNevertheless) {
                            // use if the domination by the archive is not too bad (within the epsilon box)
                            offSpringSolutionSet.add(offSpring[0]);
//                            newPointIsDominatedByOldArchiveTakeNeverthelessCounter++;
                        }
                        
                        continue;
                        
                    } else {
                        offSpringSolutionSet.add(offSpring[0]);

                        if (useEpsilonBoxesArchive) {
                            offSpring[0] = convertSolutionToEpsilonGridVectorFLOOR(offSpring[0], epsilonGridWidth);
                        }
                        
                        SolutionSet temp = new SolutionSet(1);
                        temp.add(offSpring[0]);
                        archive = archive.union(temp);
                    }  
                    
                    
                    
                    
                    
              }
            /*END generate lambda invididuals*/

//System.out.print("post offspring generation "+offSpringSolutionSet.size()+" "+offSpringSolutionSetForArchive.size());
   

/* merge population with offSpringSolutionSet */
            solutionSet = solutionSet.union(offSpringSolutionSet);              // would it be neccessary to take just the first subfront?

             // pAGE START (XXXX): archive post processing, collaboration with Shahriar Mahbub
            /* pAGE_online */
            if (useEpsilonBoxesArchive) {
                for(int i=0; i<archive.size(); i++) {
                    Solution s = archive.get(i);
                    if (true && 
                            /* the following needs a magic number for the exponent: higher value means lower pressure in the first part of the optimisation
                             */
                            Math.random()< Math.pow( (evaluations/(double)maxEvaluations) , 5) &&  // 1: linear, 10: quite low first phase and then at 90% a chance of 35% to get kicked
                            !isSolutionInPreferred(s)) {
                        archive.remove(i);
                        i--;// needs to be decremented since the previous item at i was removed
                    }
                }
            }
            // pAGE END
            
    
            /* START select mu auf of mu+lambda */
            reducePopulationToSize(solutionSet, archive, populationSize);
            /* END select mu auf of mu+lambda */

            
            
//            if (evaluations%10000==5000) System.out.println(evaluations+": newPointIsDominatedByOldArchiveCounter="+newPointIsDominatedByOldArchiveCounter + 
//                " newPointIsDominatedByOldArchiveTakeNeverthelessCounter="+newPointIsDominatedByOldArchiveTakeNeverthelessCounter);
            
        } // end of main loop

     /*   if (doOnMPICluster) {
            infoPrinter.printLotsOfValuesToFile(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, true);
            infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, true); //correct line
//            infoPrinter.printLotsOfValues(this, problem_, operators_, inputParameters_, evaluations, solutionSet, archive, indicators, true, false);
        }*/

        
//        if (false) { 
//            for (int i=0; i<1; i++) {System.out.print("sample archive point: ");
//                Solution s = archive.get(i);
//                for (int j = 0; j<s.numberOfObjectives(); j++) {
//                    System.out.print(s.getObjective(j)+" ");
//                }
//                System.out.println();
//            }
//        }
        
//        System.out.println("newPointIsDominatedByOldArchiveCounter="+newPointIsDominatedByOldArchiveCounter + 
//                " newPointIsDominatedByOldArchiveTakeNeverthelessCounter="+newPointIsDominatedByOldArchiveTakeNeverthelessCounter);
        
System.out.println(" #### GECCO 2013 Version + preferred regions ####");


        // pAGE START
        if (useEpsilonBoxesArchive) { 
            System.out.println("this is pAGE_online (reduction of archive to the preferred region DURING the runtime)");
        } else {
            System.out.println("this is pAGE_offline (reduction of archive to the preferred region AFTER the runtime)");
        }
        /* oAGE_online is */
        
        /* pAGE_offline:
        1. cut it down so that only the preferred regions are left
        2. reduce it with reducePopulationToSize so that only the best populationSize many are left for that reduced set
        Note that the population effectively only plays the role of growing the archive.
        */
        // 1. step
        if (!useEpsilonBoxesArchive) {
            boolean debugPageOffline = !true;
            if (true) // check and remove
                for(int i=0; i<archive.size(); i++) {
                    Solution s = archive.get(i);
                    if (debugPageOffline) System.out.println("offline: c="+i+" size="+archive.size()+" --> "+archive.get(i).toString());
                    if (!isSolutionInPreferred(s)) {
                        if (debugPageOffline) System.out.println("offline: -="+i+" size="+archive.size()+" --> "+archive.get(i).toString());
                        archive.remove(i);
                        if (debugPageOffline && i<archive.size()) System.out.println("offline: -    size="+archive.size()+" --> "+archive.get(i).toString());
                        i--;// needs to be decremented since the previous item at i was removed
                    }
//                    i++;
                }
            // 2. step
            if (debugPageOffline) System.out.println("DeepCopy done");
            solutionSet = (SolutionSet)DeepCopy.copy(archive);
            if (debugPageOffline) {
                System.out.println("new solutionSet (based on archive): size="+solutionSet.size());
                printDouble2DArray(solutionSet.writeObjectivesToMatrix());
            }
            reducePopulationToSize(solutionSet, archive, populationSize);
            if (debugPageOffline) {
                System.out.println("new solutionSet (reduced to size): size="+solutionSet.size());
                printDouble2DArray(solutionSet.writeObjectivesToMatrix());
            }
        }
               
        
        
        //print objective values
        System.out.println("solutionSet: size="+solutionSet.size());
        printDouble2DArray(solutionSet.writeObjectivesToMatrix());
//        System.out.println("archive: size="+archive.size());
//        printDouble2DArray(archive.writeObjectivesToMatrix());

        return solutionSet;
    } // execute

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /// Here follow the important AGE functions. For details, contact wagner@acrocon.com ///////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
        
    
    public static Solution convertSolutionToEpsilonGridVectorFLOOR(Solution s, double epsilonGridWidth) {
        Solution result = new Solution(s.getNumberOfObjectives());
        for (int i=0; i<s.getNumberOfObjectives(); i++) {
                double v = s.getObjective(i);
                result.setObjective(i, epsilonGridWidth*Math.floor( v/epsilonGridWidth)  );
        }
        return result;
    }
    public static Solution convertSolutionToEpsilonGridVectorCEILING(Solution s, double epsilonGridWidth) {
        Solution result = new Solution(s.getNumberOfObjectives());
        for (int i=0; i<s.getNumberOfObjectives(); i++) {
                double v = s.getObjective(i);
                result.setObjective(i, epsilonGridWidth*Math.ceil( v/epsilonGridWidth)  );
        }
        return result;
    }
    
    public static Solution moveEpsilonGridVectorOnceFurtherAway(Solution s, double epsilonGridWidth) {
        Solution result = new Solution(s.getNumberOfObjectives());
        for (int i=0; i<s.getNumberOfObjectives(); i++) {
                double v = s.getObjective(i);
                result.setObjective(i, (v/epsilonGridWidth + 1) * epsilonGridWidth  );
        }
        return result;
    }
    
    
    /** Compute how well a solution p approximates a solution s. This is done by
     * determining the approximation for all decision variables and objective
     * variables and then taking the maximum of these approximations (i.e. the worst).
     * @param p
     * @param s
     * @return max(all approximations "p app s")
     */
    public static double computeApproximationSolutionForSolution(Solution p, Solution s) {
        boolean debugPrint = true && debugPrintGlobal;
      if (debugPrint) System.out.print("computeApproximationSolutionForSolution/2: ");

        double delta = 0; // in order to maximize delta

        /* DV = decision variable, O = objective variable */
        double[] sOV = solutionObjectivesToDoubleArray(s);
        double[] pOV = solutionObjectivesToDoubleArray(p);


        double temp;

        // 2. compute the approximation of s by p in image (problem-specific)
        for (int i = 0; i < sOV.length; i++) {
            temp = pOV[i] - sOV[i];
            if (temp > delta) delta = temp;
            if (debugPrint) System.out.printf(temp + " ");
        }

        if (debugPrint) System.out.print(" maxDelta:" + delta + "\n");
//      if (debugPrint) System.out.printf(" maxDelta:%8.5f\n",delta);
        return delta;
    }

    /** Returns the decision variables of a solution
     * @param s
     * @return double[]
     */
    public static double[] solutionDecisionVariablesToDoubleArray(Solution s) {
        Variable[] varsV = s.getDecisionVariables();
        double[] result = new double[varsV.length];
        for (int j = 0; j < varsV.length; j++) {
            try {
                result[j] = varsV[j].getValue();
            } catch (JMException ex) {
                System.out.println("problem in solutionDecisionVariablesToDoubleArray");
            }
        }
        return result;
    }

    /** Returns the objective variables of a solution
     * @param s
     * @return
     */
    public static double[] solutionObjectivesToDoubleArray(Solution s) {
        int numberOfObjectives = s.getNumberOfObjectives();
        double[] result = new double[numberOfObjectives];
        for (int j = 0; j < numberOfObjectives; j++) {
            result[j] = s.getObjective(j);
        }
        return result;
    }

    /** Prints a double[][] on screen.
     * @param a
     */
    public static void printDouble2DArray(double[][] a) {
        for (double[] row : a) {
            for (double d : row) {
                System.out.printf("%8.5f ", d);
            }
            System.out.println();
        }
    }




    // the followin is based on computeFitnesses/2
    /** Compute how well the population approximates the archive, see Karl's email from 07.12.2010
     *
     * @param population
     * @param archive
//     * @return
     */
    public void reducePopulationToSize(SolutionSet population, SolutionSet archive, int targetSize) {
        
        if (population.size()==1) return;
        
        boolean debugPrint = true && debugPrintGlobal;
        boolean debugPrintAdditional = false;
        if (debugPrint) System.out.println("computeFitnesses/2: ");

        // the following array stores the maximum approximations for which a population point is responsible
        int[] whichPopPointIsResponsible = new int[archive.size()];
        int[] whichPopPointIsResponsibleSecondBest = new int[archive.size()];

        SolutionSet archiveFront = archive;

        double[] results = new double[archiveFront.size()];

        if (debugPrint) System.out.println("population.size=" + population.size() + " archive.size=" + archive.size() + " front.size=" + archiveFront.size());

        /* store all the "pop app arc" in the following array:
         * for each archive point it is stored how well a population point approximates it
         */

        /* store all minimums of approximations per archive point in this array */
        double[] eps1 = new double[archiveFront.size()];
        double[] eps2 = new double[archiveFront.size()];

        // compute approximation for each non-dominated point of the archive
        // to find out how well it is approximated by the population points
        for (int i = 0; i < archiveFront.size(); i++) {
            Solution s = archiveFront.get(i);   // non-dominated archive-point

            double deltaForThisSolution = Double.MAX_VALUE; //minimize this value
            double deltaForThisSolutionSecondBest = Double.MAX_VALUE; //minimize this value
            for (int j = 0; j < population.size(); j++) {
                Solution p = population.get(j);

                // compute how an element p of the front(population) approximates an element s of the archive in domain and image
                double deltaForThisSolutionCurrent = computeApproximationSolutionForSolution(p, s); // population element p, archiveFront element s

                // if better then currentBest
                if (deltaForThisSolutionCurrent < deltaForThisSolution) {
                    deltaForThisSolutionSecondBest = deltaForThisSolution;      // old best becomes secondBest
                    deltaForThisSolution = deltaForThisSolutionCurrent;         // new best becomes best
                    whichPopPointIsResponsibleSecondBest[i] = whichPopPointIsResponsible[i];
                    whichPopPointIsResponsible[i] = j;
                } else {
                // if better then currentSecondBest
                    if (deltaForThisSolutionCurrent < deltaForThisSolutionSecondBest) {
                        deltaForThisSolutionSecondBest = deltaForThisSolutionCurrent;
                        whichPopPointIsResponsibleSecondBest[i] = j;
                    }
                }
            }

            // save the minimal approximation
            results[i] = deltaForThisSolution;
            eps1[i] = deltaForThisSolution;
            eps2[i] = deltaForThisSolutionSecondBest;

            if (debugPrint) {
                System.out.printf("-deltaForThisSolution: %8.5f ", deltaForThisSolution);
                System.out.println("(=" + deltaForThisSolution + ")"); // same output as line before, just more precise
            }
        }

        // keep track whether some point is still in the population
        boolean[] pIsInCurrentPop = new boolean[population.size()];
        for (int i = 0; i<pIsInCurrentPop.length; i++)
            pIsInCurrentPop[i] = true;

        //now determine the val(p)
        double[] val = new double[population.size()];
        double minVal = Double.MAX_VALUE;
        int minValIndex = 0;
        for (int i = 0; i<eps1.length; i++) {
            double eps2a = eps2[i];
            int p = whichPopPointIsResponsible[i];
                if (eps2a>val[p])
                    val[p] = eps2a;            
        }

        for (int i=0; i<val.length; i++) {
                if (val[i]<minVal) {
                        minVal = val[i];
                        minValIndex = i;
                }
        }
        
        int popCounter=population.size();
        if (debugPrintAdditional) System.out.println("popCounter"+popCounter);
        
        pIsInCurrentPop[minValIndex] = false;
        popCounter--;
        if (debugPrintAdditional) System.out.println("popCounter"+popCounter);

        if (debugPrintAdditional) if (true) {
            System.out.println("unsorted eps1:    " + Arrays.toString(eps1));
            System.out.println("responsible:      " + Arrays.toString(whichPopPointIsResponsible));
            System.out.println("unsorted eps2:    " + Arrays.toString(eps2));
            System.out.println("responsible:      " + Arrays.toString(whichPopPointIsResponsibleSecondBest));
            System.out.println("val        :      " + Arrays.toString(val));
            System.out.println("pIsInCurrentPop:  " + Arrays.toString(pIsInCurrentPop));
//            System.out.println("pIsP1a     :      " + Arrays.toString(pIsP1a));
            System.out.println("minVal="+minVal + " index=" +minValIndex + " population.size()now="+ (population.size()));
//            System.out.println("maxAppForPopPoint:" + Arrays.toString(maxAppForPopPoint));
        }

        while (popCounter > targetSize) {
            HashSet whichEpsToUpdate = new HashSet();

            for (int i = 0; i<whichPopPointIsResponsible.length; i++) {
                if (minValIndex == whichPopPointIsResponsible[i]) whichEpsToUpdate.add(i);
                if (minValIndex == whichPopPointIsResponsibleSecondBest[i]) whichEpsToUpdate.add(i);
            }

            Iterator it = whichEpsToUpdate.iterator();
            while (it.hasNext()) {
                int i = (Integer)it.next();
                if (debugPrintAdditional) System.out.print(i+" ");
            }
            if (debugPrintAdditional) System.out.println("");


            it = whichEpsToUpdate.iterator();
            while (it.hasNext()) {
                int i = (Integer)it.next();

                Solution s = archiveFront.get(i);   // non-dominated archive-point

                double deltaForThisSolution = Double.MAX_VALUE; //minimize this value
                double deltaForThisSolutionSecondBest = Double.MAX_VALUE; //minimize this value
                for (int j = 0; j < population.size(); j++) {

                    // skip those of the old population that are no longer in the current population
                    if (!pIsInCurrentPop[j]) continue;

                    Solution p = population.get(j);

                    // compute how an element p of the front(population) approximates an element s of the archive in domain and image
                    double deltaForThisSolutionCurrent = computeApproximationSolutionForSolution(p, s); // population element p, archiveFront element s

                    // if better then currentBest
                    if (deltaForThisSolutionCurrent < deltaForThisSolution) {
                        deltaForThisSolutionSecondBest = deltaForThisSolution;      // old best becomes secondBest
                        deltaForThisSolution = deltaForThisSolutionCurrent;         // new best becomes best
                        whichPopPointIsResponsibleSecondBest[i] = whichPopPointIsResponsible[i];
                        whichPopPointIsResponsible[i] = j;
                    } else {
                    // if better then currentSecondBest
                        if (deltaForThisSolutionCurrent < deltaForThisSolutionSecondBest) {
                            deltaForThisSolutionSecondBest = deltaForThisSolutionCurrent;
                            whichPopPointIsResponsibleSecondBest[i] = j;
                        }
                    }

                }

                // save the minimal approximation
                results[i] = deltaForThisSolution;
                eps1[i] = deltaForThisSolution;
                eps2[i] = deltaForThisSolutionSecondBest;

                if (debugPrint) {
                    System.out.printf("-deltaForThisSolution: %8.5f ", deltaForThisSolution);
                    System.out.println("(=" + deltaForThisSolution + ")"); // same output as line before, just more precise
                }
            }

            //now determine the val(p)
            minVal = Double.MAX_VALUE;  // do not use 0 here as one point might dominate another
            minValIndex = 0;
            for (int i = 0; i<eps1.length; i++) {
                double eps2a = eps2[i];
                int p = whichPopPointIsResponsible[i];
                    if (eps2a>val[p])
                        val[p] = eps2a;
            }

            for (int i=0; i<val.length; i++) {
                if (pIsInCurrentPop[i])
                    if (val[i]<minVal) {
                            minVal = val[i];
                            minValIndex = i;
                    }

            }
            popCounter--;
            if (debugPrintAdditional) System.out.println("popCounter"+popCounter);
            pIsInCurrentPop[minValIndex] = false;

            if (debugPrintAdditional) if (true) {
                System.out.println("unsorted eps1:    " + Arrays.toString(eps1));
                System.out.println("responsible:      " + Arrays.toString(whichPopPointIsResponsible));
                System.out.println("unsorted eps2:    " + Arrays.toString(eps2));
                System.out.println("responsible:      " + Arrays.toString(whichPopPointIsResponsibleSecondBest));
                System.out.println("val        :      " + Arrays.toString(val));
                System.out.println("pIsInCurrentPop:  " + Arrays.toString(pIsInCurrentPop));
    //            System.out.println("pIsP1a     :      " + Arrays.toString(pIsP1a));
                System.out.println("minVal="+minVal + " index=" +minValIndex + " population.size()now="+ (population.size()));
    //            System.out.println("maxAppForPopPoint:" + Arrays.toString(maxAppForPopPoint));
            }
        }
       

        // form new population:
        for (int i = pIsInCurrentPop.length - 1; i>=0; i--) {
            if (!pIsInCurrentPop[i]) population.remove(i);
            if (debugPrintAdditional) System.out.println("population.size()"+population.size());
        }
    } //end
    
    
    
    
    /** Minimal function to quickly execute the algorithm
     * @param args
     */
    public static void main(String[] args) {
        
        
        try {
            if (args.length == 0) {


                File test = new File(".");
                System.out.println(test.getAbsolutePath());
              
                //problem, Pareto front, pop size, doCrossover, doMutation, epsilonGrid parameter, runOnCluster true/false (less/more outputs [a powerfull switch!]), evaluations, parentSelection, outputfolder
//              pAGE_main.main(new String[]{"DTLZ3_2D", "originalParetoFronts\\DTLZ2.2D.1000.pf", "30", "true", "true", "0.01", "true","100000","BinaryTournament","foo"});   
                pAGE_main.main(new String[]{"DTLZ4", "paretoFronts\\DTLZ4.3D.pf", "40", "true", "true", "0.00", "false","50000","BinaryTournament","foo"});   
//                pAGE_main.main(new String[]{"ZDT2", "originalParetoFronts\\ZDT2.pf", "30", "true", "true", "0.01", "true","10000","BinaryTournament2","foo"});   
//              pAGE_main.main(new String[]{"WFG1_3D", "originalParetoFronts\\WFG1.3D.pf", "100", "true", "true", "0.01", "false","100000","BinaryTournament","foo"});   
//              pAGE_main.main(new String[]{"DTLZ4_3D", "originalParetoFronts\\DTLZ2.3D.1000.pf", "100", "true", "true", "0.01", "false","100000","RandomSelection","foo"});   
            }
            else if (args.length == 7) {
                // in case it is called via the commandline with parameters...
                pAGE_main.mainORIG(new String[]{
                            args[0], // problem
                            "originalParetoFronts/" + args[1], // pf-file, without the directory
                            args[2], // mu
                            "true",  // doXO
                            "true",  // doMUT
                            args[5], // epsilonGridWidth
                            "true",  // do outputs for runs on MPI cluster
                            args[3], // max evaluations
                            args[4], // selection strategy
                            args[6]  // subDirectory name for infoprinter
                        });
            } else {
                System.out.println("unsuitable number of parameters. EXIT.");
            }

        } catch (Exception ex) {
        }
    }
}