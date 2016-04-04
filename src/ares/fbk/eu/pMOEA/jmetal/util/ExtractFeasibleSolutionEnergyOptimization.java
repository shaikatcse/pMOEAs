package ares.fbk.eu.pMOEA.jmetal.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.LineNumberReader;
import java.util.StringTokenizer;

import ares.fbk.eu.pMOEA.jmetal.problems.OptimizeEnergyPLANAalborg.EnergyPLANProblemAalborg2ObjectivesWith1EnergyPLANEvolution;
import ares.fbk.eu.pMOEA.jmetal.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal.core.Operator;
import ares.fbk.eu.pMOEA.jmetal.core.Problem;
import ares.fbk.eu.pMOEA.jmetal.core.Solution;
import ares.fbk.eu.pMOEA.jmetal.core.SolutionSet;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;

public class ExtractFeasibleSolutionEnergyOptimization {

	public ExtractFeasibleSolutionEnergyOptimization() {
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) throws JMException,
			SecurityException, IOException, ClassNotFoundException {

		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator
		Operator selection; // Selection operator

		SolutionSet population;

		FileInputStream fos = new FileInputStream(args[0]);
		InputStreamReader isr = new InputStreamReader(fos);
		BufferedReader br = new BufferedReader(isr);

		LineNumberReader lnr = new LineNumberReader(new FileReader(new File(
				args[0])));
		lnr.skip(Long.MAX_VALUE);
		int totalNumberOfIndividuals = lnr.getLineNumber();
		lnr.close();

		String line;

		population = new SolutionSet(totalNumberOfIndividuals);

		problem = new EnergyPLANProblemAalborg2ObjectivesWith1EnergyPLANEvolution("Real");

		Solution newSolution;
		for (int i = 0; i < totalNumberOfIndividuals; i++) {
			line = br.readLine();
			StringTokenizer st = new StringTokenizer(line);
			newSolution = new Solution(problem);
			int j = 0;
			while (st.hasMoreTokens()) {
				newSolution.getDecisionVariables()[j].setValue(Double
						.parseDouble(st.nextToken()));
				j++;
			}

			problem.evaluate(newSolution);
			problem.evaluateConstraints(newSolution);
			population.add(newSolution);
		}

	/*	File filetmp = new File(args[0]);
		String arrayStr[]=args[0].split("\\\\");*/
		
		//File file = new File(filetmp.getParent() + "\\allParameters_"+arrayStr[arrayStr.length-1] +".txt");
		population.printFeasibleFUN("FUN0_" +"Feasible");
		//population.printFeasibleVAR(filetmp.getParent() + "\\"+arrayStr[arrayStr.length-1] +"Feasible");
		
		br.close();
	}

}
