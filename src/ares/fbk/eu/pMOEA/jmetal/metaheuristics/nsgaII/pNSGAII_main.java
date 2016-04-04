//  NSGAII_main.java
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

import ares.fbk.eu.pMOEA.jmetal.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal.core.Operator;
import ares.fbk.eu.pMOEA.jmetal.core.Problem;
import ares.fbk.eu.pMOEA.jmetal.core.SolutionSet;
import ares.fbk.eu.pMOEA.jmetal.operators.crossover.CrossoverFactory;
import ares.fbk.eu.pMOEA.jmetal.operators.mutation.MutationFactory;
import ares.fbk.eu.pMOEA.jmetal.operators.selection.SelectionFactory;
import ares.fbk.eu.pMOEA.jmetal.problems.ProblemFactory;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.ZDT1;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.ZDT2;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.ZDT3;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.ZDT4;
import ares.fbk.eu.pMOEA.jmetal.problems.ZDT.ZDT6;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.QualityIndicator;
import ares.fbk.eu.pMOEA.jmetal.util.Configuration;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

import ares.fbk.eu.pMOEA.jmetal.util.PseudoRandom;
import ares.fbk.eu.pMOEA.jmetal.util.RandomGenerator;

/**
 * Class to configure and execute the NSGA-II algorithm.
 * 
 * Besides the classic NSGA-II, a steady-state version (ssNSGAII) is also
 * included (See: J.J. Durillo, A.J. Nebro, F. Luna and E. Alba "On the Effect
 * of the Steady-State Selection Scheme in Multi-Objective Genetic Algorithms"
 * 5th International Conference, EMO 2009, pp: 183-197. April 2009)
 */

public class pNSGAII_main {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object

	/**
	 * @param args
	 *            Command line arguments.
	 * @throws JMException
	 * @throws IOException
	 * @throws SecurityException
	 *             Usage: three options -
	 *             ares.fbk.eu.pMOEA.jmetal.metaheuristics.nsgaII.NSGAII_main -
	 *             ares.fbk.eu.pMOEA.jmetal.metaheuristics.nsgaII.NSGAII_main
	 *             problemName -
	 *             ares.fbk.eu.pMOEA.jmetal.metaheuristics.nsgaII.NSGAII_main
	 *             problemName paretoFrontFile
	 */
	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; // Selection operator

		HashMap parameters; // Operator parameters

		QualityIndicator indicators; // Object to get quality indicators

		// seeds
		long seeds[] = {  831779, 673381, 218419, 657845, 162571, 862985, 914300,
				923266, 459019, 176711, 717128, 725473, 890286, 113706, 344376,
				853009, 813403, 300225, 871571, 467468, 172993, 908111, 742682,
				817558, 578279, 317033, 904701, 537469, 900260, 649576, 253540,
				470826, 266413, 615825, 603054, 550654, 722492, 868682, 991831,
				564213, 547800, 294078, 682367, 818048, 511910, 862665, 705631,
				285055, 157165, 654098 };

		// Logger object and file to store log messages
		logger_ = Configuration.logger_;
		fileHandler_ = new FileHandler("NSGAII_main.log");
		logger_.addHandler(fileHandler_);
		
		String path = "WithTrack\\pNSGAII";

		for (int numnerOfRun = 0; numnerOfRun < 1; numnerOfRun++) {

			PseudoRandom
					.setRandomGenerator(new RandomGenerator(seeds[numnerOfRun]));

			problem = new ZDT1("ArrayReal");

			List<Double[]> regions = new ArrayList<Double[]>(3);
			regions.add(new Double[] { 0.8, 0.95 });
			//regions.add(new Double[] { 0.40, 0.50 });
			regions.add(new Double[] { 0.15, 0.20 });

			int numberOfSolutionPerRegion[] = { 6, 6 };

			//algorithm = new pNSGAII(problem, seeds[numnerOfRun], path, regions.size());
			algorithm = new pNSGAII(problem);
			// algorithm = new ssNSGAII(problem);

			// Algorithm parameters
			int populationSize = 0;
			for (int i = 0; i < numberOfSolutionPerRegion.length; i++) {
				populationSize += numberOfSolutionPerRegion[i];
			}
			algorithm.setInputParameter("populationSize", populationSize);
			algorithm.setInputParameter("maxEvaluations", 8400);
			algorithm.setInputParameter("numberOfSolutionsPerRegion",
					numberOfSolutionPerRegion);
			algorithm.setInputParameter("regions", regions);
			algorithm.setInputParameter("paretoFrontFilePath", "C:\\Users\\mahbub\\Documents\\GitHub\\pMOEA jMetal4.5\\paretoFronts\\ZDT2.pf");

			// Mutation and Crossover for Real codification
			parameters = new HashMap();
			parameters.put("probability", 0.9);
			parameters.put("distributionIndex", 20.0);
			crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover",
					parameters);

			parameters = new HashMap();
			parameters.put("probability", 1.0 / problem.getNumberOfVariables());
			parameters.put("distributionIndex", 20.0);
			mutation = MutationFactory.getMutationOperator(
					"PolynomialMutation", parameters);

			// Selection Operator
			parameters = new HashMap();
			parameters.put("numberOfTournaments", 5);
			parameters.put("numberOfSolutionsPerRegion",
					numberOfSolutionPerRegion);
			parameters.put("numberOfRegions", 2);
			selection = SelectionFactory.getSelectionOperator(
					"PreferredTournamentSelection", parameters);

			// Add the operators to the algorithm
			algorithm.addOperator("crossover", crossover);
			algorithm.addOperator("mutation", mutation);
			algorithm.addOperator("selection", selection);

			// Add the indicator object to the algorithm
			//algorithm.setInputParameter("indicators", indicators);

			// Execute the Algorithm
			long initTime = System.currentTimeMillis();
			SolutionSet population = algorithm.execute();
			long estimatedTime = System.currentTimeMillis() - initTime;

			// Result messages
			logger_.info("Total execution time: " + estimatedTime + "ms");
			logger_.info("Variables values have been writen to file VAR");
			population.printVariablesToFile(path + "\\VAR_seed_"
					+ seeds[numnerOfRun]);
			logger_.info("Objectives values have been writen to file FUN");
			population.printObjectivesToFile(path + "\\FUN_seed_"
					+ seeds[numnerOfRun]);

			
			
		} // if
	} // main
} // NSGAII_main
