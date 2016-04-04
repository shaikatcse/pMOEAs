//  NSGAIIStudy.java
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

package ares.fbk.eu.pMOEA.jmetal.experiments.studies;

import ares.fbk.eu.pMOEA.jmetal.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal.experiments.Experiment;
import ares.fbk.eu.pMOEA.jmetal.experiments.Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.NSGAII_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.pNSGAII_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.util.Friedman;
import ares.fbk.eu.pMOEA.jmetal.experiments.util.Statistics;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Epsilon;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.GeneralizedSpread;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.InvertedGenerationalDistance;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Spread;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.pHypervolume;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.collections4.keyvalue.MultiKey;
import org.apache.commons.collections4.map.*;



/**
 * Class implementing an example of experiment using NSGA-II as base algorithm.
 * The experiment consisting in studying the effect of the crossover probability
 * in NSGA-II.
 */
public class CSVFileGenrator extends Experiment {
	

	/*String [] row={//pAGEStudyZDT
					"pAGE.F.HV.M", "pAGE.F.EPS.M", "pAGE.F.IGD.M", "pAGE.F.HV.Md", "pAGE.F.EPS.Md", "pAGE.F.IGD.Md", "pAGE.F.HV.St", "pAGE.F.EPS.St", "pAGE.F.IGD.St", "pAGE.F.HV.IQR", "pAGE.F.EPS.IQR", "pAGE.F.IGD.IQR",
					"pAGE.O.HV.M", "pAGE.O.EPS.M", "pAGE.O.IGD.M", "pAGE.O.HV.Md", "pAGE.O.EPS.Md", "pAGE.O.IGD.Md", "pAGE.O.HV.St", "pAGE.O.EPS.St", "pAGE.O.IGD.St", "pAGE.O.HV.IQR", "pAGE.O.EPS.IQR", "pAGE.O.IGD.IQR",
					"AGE.1.HV.M", "AGE.1.EPS.M", "AGE.1.IGD.M", "AGE.1.HV.Md", "AGE.1.EPS.Md", "AGE.1.IGD.Md", "AGE.1.HV.St", "AGE.1.EPS.St", "AGE.1.IGD.St", "AGE.1.HV.IQR", "AGE.1.EPS.IQR", "AGE.1.IGD.IQR",
					"AGE.2.HV.M", "AGE.2.EPS.M", "AGE.2.IGD.M", "AGE.2.HV.Md", "AGE.2.EPS.Md", "AGE.2.IGD.Md", "AGE.2.HV.St", "AGE.2.EPS.St", "AGE.2.IGD.St", "AGE.2.HV.IQR", "AGE.2.EPS.IQR", "AGE.2.IGD.IQR",
					"AGE.3.HV.M", "AGE.3.EPS.M", "AGE.3.IGD.M", "AGE.3.HV.Md", "AGE.3.EPS.Md", "AGE.3.IGD.Md", "AGE.3.HV.St", "AGE.3.EPS.St", "AGE.3.IGD.St", "AGE.3.HV.IQR", "AGE.3.EPS.IQR", "AGE.3.IGD.IQR",
					
					//pAGEStudyZDTMaxEvo
					"pAGE.F.Max.HV.M", "pAGE.F.Max.EPS.M", "pAGE.F.Max.IGD.M", "pAGE.F.Max.HV.Md", "pAGE.F.Max.EPS.Md", "pAGE.F.Max.IGD.Md", "pAGE.F.Max.HV.St", "pAGE.F.Max.EPS.St", "pAGE.F.Max.IGD.St", "pAGE.F.Max.HV.IQR", "pAGE.F.Max.EPS.IQR", "pAGE.F.Max.IGD.IQR",
					"pAGE.O.Max.HV.M", "pAGE.O.Max.EPS.M", "pAGE.O.Max.IGD.M", "pAGE.O.Max.HV.Md", "pAGE.O.Max.EPS.Md", "pAGE.O.Max.IGD.Md", "pAGE.O.Max.HV.St", "pAGE.O.Max.EPS.St", "pAGE.O.Max.IGD.St", "pAGE.O.Max.HV.IQR", "pAGE.O.Max.EPS.IQR", "pAGE.O.Max.IGD.IQR",
					
					
					//pNSGAIIStudyZDT
					"pNSGAII.HV.M", "pNSGAII.EPS.M", "pNSGAII.IGD.M", "pNSGAII.HV.Md", "pNSGAII.EPS.Md", "pNSGAII.IGD.Md", "pNSGAII.HV.St", "pNSGAII.EPS.St", "pNSGAII.IGD.St", "pNSGAII.HV.IQR", "pNSGAII.EPS.IQR", "pNSGAII.IGD.IQR",
					"NSGAII.1.HV.M", "NSGAII.1.EPS.M", "NSGAII.1.IGD.M", "NSGAII.1.HV.Md", "NSGAII.1.EPS.Md", "NSGAII.1.IGD.Md", "NSGAII.1.HV.St", "NSGAII.1.EPS.St", "NSGAII.1.IGD.St", "NSGAII.1.HV.IQR", "NSGAII.1.EPS.IQR", "NSGAII.1.IGD.IQR",
					"NSGAII.2.HV.M", "NSGAII.2.EPS.M", "NSGAII.2.IGD.M", "NSGAII.2.HV.Md", "NSGAII.2.EPS.Md", "NSGAII.2.IGD.Md", "NSGAII.2.HV.St", "NSGAII.2.EPS.St", "NSGAII.2.IGD.St", "NSGAII.2.HV.IQR", "NSGAII.2.EPS.IQR", "NSGAII.2.IGD.IQR",
					"NSGAII.3.HV.M", "NSGAII.3.EPS.M", "NSGAII.3.IGD.M", "NSGAII.3.HV.Md", "NSGAII.3.EPS.Md", "NSGAII.3.IGD.Md", "NSGAII.3.HV.St", "NSGAII.3.EPS.St", "NSGAII.3.IGD.St", "NSGAII.3.HV.IQR", "NSGAII.3.EPS.IQR", "NSGAII.3.IGD.IQR",
					
					//pNSGAIIStudyZDTMaxEvo
					"pNSGAII.Max.HV.M", "pNSGAII.Max.EPS.M", "pNSGAII.Max.IGD.M", "pNSGAII.Max.HV.Md", "pNSGAII.Max.EPS.Md", "pNSGAII.Max.IGD.Md", "pNSGAII.Max.HV.St", "pNSGAII.Max.EPS.St", "pNSGAII.Max.IGD.St", "pNSGAII.Max.HV.IQR", "pNSGAII.Max.EPS.IQR", "pNSGAII.Max.IGD.IQR",
					
					//pSPEA2StudyZDT
					"pSPEA2.HV.M", "pSPEA2.EPS.M", "pSPEA2.IGD.M", "pSPEA2.HV.Md", "pSPEA2.EPS.Md", "pSPEA2.IGD.Md", "pSPEA2.HV.St", "pSPEA2.EPS.St", "pSPEA2.IGD.St", "pSPEA2.HV.IQR", "pSPEA2.EPS.IQR", "pSPEA2.IGD.IQR",
					"SPEA2.1.HV.M", "SPEA2.1.EPS.M", "SPEA2.1.IGD.M", "SPEA2.1.HV.Md", "SPEA2.1.EPS.Md", "SPEA2.1.IGD.Md", "SPEA2.1.HV.St", "SPEA2.1.EPS.St", "SPEA2.1.IGD.St", "SPEA2.1.HV.IQR", "SPEA2.1.EPS.IQR", "SPEA2.1.IGD.IQR",
					"SPEA2.2.HV.M", "SPEA2.2.EPS.M", "SPEA2.2.IGD.M", "SPEA2.2.HV.Md", "SPEA2.2.EPS.Md", "SPEA2.2.IGD.Md", "SPEA2.2.HV.St", "SPEA2.2.EPS.St", "SPEA2.2.IGD.St", "SPEA2.2.HV.IQR", "SPEA2.2.EPS.IQR", "SPEA2.2.IGD.IQR",
					"SPEA2.3.HV.M", "SPEA2.3.EPS.M", "SPEA2.3.IGD.M", "SPEA2.3.HV.Md", "SPEA2.3.EPS.Md", "SPEA2.3.IGD.Md", "SPEA2.3.HV.St", "SPEA2.3.EPS.St", "SPEA2.3.IGD.St", "SPEA2.3.HV.IQR", "SPEA2.3.EPS.IQR", "SPEA2.3.IGD.IQR",
	
					//pNSGAIIStudyZDTMaxEvo
					"pSPEA2.Max.HV.M", "pSPEA2.Max.EPS.M", "pSPEA2.Max.IGD.M", "pSPEA2.Max.HV.Md", "pSPEA2.Max.EPS.Md", "pSPEA2.Max.IGD.Md", "pSPEA2.Max.HV.St", "pSPEA2.Max.EPS.St", "pSPEA2.Max.IGD.St", "pSPEA2.Max.HV.IQR", "pSPEA2.Max.EPS.IQR", "pSPEA2.Max.IGD.IQR"
	};
	
	String [] col ={"ZDT1.R0", "ZDT2.R0", "ZDT3.R0", "ZDT4.R0", "ZDT6.R0",
					"ZDT1.R1", "ZDT2.R1", "ZDT3.R1", "ZDT4.R1", "ZDT6.R1",
					"ZDT1.R2", "ZDT2.R2", "ZDT3.R2", "ZDT4.R2", "ZDT6.R2"
				};*/
	
	String [] row={//pAGEStudyDTLZ
			"pAGE.F.1.HV.M", "pAGE.F.1.EPS.M", "pAGE.F.1.IGD.M", "pAGE.F.1.HV.Md", "pAGE.F.1.EPS.Md", "pAGE.F.1.IGD.Md", "pAGE.F.1.HV.St", "pAGE.F.1.EPS.St", "pAGE.F.1.IGD.St", "pAGE.F.1.HV.IQR", "pAGE.F.1.EPS.IQR", "pAGE.F.1.IGD.IQR",
			"pAGE.F.2.HV.M", "pAGE.F.2.EPS.M", "pAGE.F.2.IGD.M", "pAGE.F.2.HV.Md", "pAGE.F.2.EPS.Md", "pAGE.F.2.IGD.Md", "pAGE.F.2.HV.St", "pAGE.F.2.EPS.St", "pAGE.F.2.IGD.St", "pAGE.F.2.HV.IQR", "pAGE.F.2.EPS.IQR", "pAGE.F.2.IGD.IQR",
			
			"pAGE.O.1.HV.M", "pAGE.O.1.EPS.M", "pAGE.O.1.IGD.M", "pAGE.O.1.HV.Md", "pAGE.O.1.EPS.Md", "pAGE.O.1.IGD.Md", "pAGE.O.1.HV.St", "pAGE.O.1.EPS.St", "pAGE.O.1.IGD.St", "pAGE.O.1.HV.IQR", "pAGE.O.1.EPS.IQR", "pAGE.O.1.IGD.IQR",
			"pAGE.O.2.HV.M", "pAGE.O.2.EPS.M", "pAGE.O.2.IGD.M", "pAGE.O.2.HV.Md", "pAGE.O.2.EPS.Md", "pAGE.O.2.IGD.Md", "pAGE.O.2.HV.St", "pAGE.O.2.EPS.St", "pAGE.O.2.IGD.St", "pAGE.O.2.HV.IQR", "pAGE.O.2.EPS.IQR", "pAGE.O.2.IGD.IQR",
			
			"AGE.1.HV.M", "AGE.1.EPS.M", "AGE.1.IGD.M", "AGE.1.HV.Md", "AGE.1.EPS.Md", "AGE.1.IGD.Md", "AGE.1.HV.St", "AGE.1.EPS.St", "AGE.1.IGD.St", "AGE.1.HV.IQR", "AGE.1.EPS.IQR", "AGE.1.IGD.IQR",
			"AGE.2.HV.M", "AGE.2.EPS.M", "AGE.2.IGD.M", "AGE.2.HV.Md", "AGE.2.EPS.Md", "AGE.2.IGD.Md", "AGE.2.HV.St", "AGE.2.EPS.St", "AGE.2.IGD.St", "AGE.2.HV.IQR", "AGE.2.EPS.IQR", "AGE.2.IGD.IQR",
			
			
			
			//pNSGAIIStudyDTLZ
			"pNSGAII.1.HV.M", "pNSGAII.1.EPS.M", "pNSGAII.1.IGD.M", "pNSGAII.1.HV.Md", "pNSGAII.1.EPS.Md", "pNSGAII.1.IGD.Md", "pNSGAII.1.HV.St", "pNSGAII.1.EPS.St", "pNSGAII.1.IGD.St", "pNSGAII.1.HV.IQR", "pNSGAII.1.EPS.IQR", "pNSGAII.1.IGD.IQR",
			"pNSGAII.2.HV.M", "pNSGAII.2.EPS.M", "pNSGAII.2.IGD.M", "pNSGAII.2.HV.Md", "pNSGAII.2.EPS.Md", "pNSGAII.2.IGD.Md", "pNSGAII.2.HV.St", "pNSGAII.2.EPS.St", "pNSGAII.2.IGD.St", "pNSGAII.2.HV.IQR", "pNSGAII.2.EPS.IQR", "pNSGAII.2.IGD.IQR",
			
			"NSGAII.1.HV.M", "NSGAII.1.EPS.M", "NSGAII.1.IGD.M", "NSGAII.1.HV.Md", "NSGAII.1.EPS.Md", "NSGAII.1.IGD.Md", "NSGAII.1.HV.St", "NSGAII.1.EPS.St", "NSGAII.1.IGD.St", "NSGAII.1.HV.IQR", "NSGAII.1.EPS.IQR", "NSGAII.1.IGD.IQR",
			"NSGAII.2.HV.M", "NSGAII.2.EPS.M", "NSGAII.2.IGD.M", "NSGAII.2.HV.Md", "NSGAII.2.EPS.Md", "NSGAII.2.IGD.Md", "NSGAII.2.HV.St", "NSGAII.2.EPS.St", "NSGAII.2.IGD.St", "NSGAII.2.HV.IQR", "NSGAII.2.EPS.IQR", "NSGAII.2.IGD.IQR",
			
			
			//pSPEA2StudyDTLZ
			"pSPEA2.1.HV.M", "pSPEA2.1.EPS.M", "pSPEA2.1.IGD.M", "pSPEA2.1.HV.Md", "pSPEA2.1.EPS.Md", "pSPEA2.1.IGD.Md", "pSPEA2.1.HV.St", "pSPEA2.1.EPS.St", "pSPEA2.1.IGD.St", "pSPEA2.1.HV.IQR", "pSPEA2.1.EPS.IQR", "pSPEA2.1.IGD.IQR",
			"pSPEA2.2.HV.M", "pSPEA2.2.EPS.M", "pSPEA2.2.IGD.M", "pSPEA2.2.HV.Md", "pSPEA2.2.EPS.Md", "pSPEA2.2.IGD.Md", "pSPEA2.2.HV.St", "pSPEA2.2.EPS.St", "pSPEA2.2.IGD.St", "pSPEA2.2.HV.IQR", "pSPEA2.2.EPS.IQR", "pSPEA2.2.IGD.IQR",
			
			"SPEA2.1.HV.M", "SPEA2.1.EPS.M", "SPEA2.1.IGD.M", "SPEA2.1.HV.Md", "SPEA2.1.EPS.Md", "SPEA2.1.IGD.Md", "SPEA2.1.HV.St", "SPEA2.1.EPS.St", "SPEA2.1.IGD.St", "SPEA2.1.HV.IQR", "SPEA2.1.EPS.IQR", "SPEA2.1.IGD.IQR",
			"SPEA2.2.HV.M", "SPEA2.2.EPS.M", "SPEA2.2.IGD.M", "SPEA2.2.HV.Md", "SPEA2.2.EPS.Md", "SPEA2.2.IGD.Md", "SPEA2.2.HV.St", "SPEA2.2.EPS.St", "SPEA2.2.IGD.St", "SPEA2.2.HV.IQR", "SPEA2.2.EPS.IQR", "SPEA2.2.IGD.IQR",

};

String [] col ={"DTLZ2.R0", "DTLZ3.R0",
			"DTLZ2.R1", "DTLZ3.R1",
			"DTLZ2.R2", "DTLZ3.R2"
		};
	
	
	MultiKeyMap<MultiKey<String>, Double> table = new MultiKeyMap();
	

	public List<Double[]> regions = new ArrayList<Double[]>(3);

	public CSVFileGenrator(List<Double[]> regions) {
		this.regions = regions;
	}

	
	/**
	 * Generate the Quality Indicators
	 */

	protected void printClearPage(String fileName) throws IOException {
		FileWriter os = new FileWriter(fileName, true);
		os.write("\\clearpage" + "\n");

		os.close();
	}

	protected void printBoldText(String fileName, String text)
			throws IOException {
		FileWriter os = new FileWriter(fileName, true);
		os.write("\\textbf{" + text + "}" + "\n");

		os.close();
	}

	

	public void generateLatexTables() throws FileNotFoundException, IOException {
	

		Vector[][][] data = new Vector[indicatorList_.length][][];
		for (int regionNumber = 0; regionNumber < regions.size(); regionNumber++) {
			for (int indicator = 0; indicator < indicatorList_.length; indicator++) {
				// A data vector per problem
				data[indicator] = new Vector[problemList_.length][];

				for (int problem = 0; problem < problemList_.length; problem++) {
					data[indicator][problem] = new Vector[algorithmNameList_.length];

					for (int algorithm = 0; algorithm < algorithmNameList_.length; algorithm++) {
						data[indicator][problem][algorithm] = new Vector();

						String directory = experimentBaseDirectory_;
						directory += "/data/";
						directory += "/" + algorithmNameList_[algorithm];
						directory += "/" + problemList_[problem];
						directory += "/" + indicatorList_[indicator] + ".r."
								+ regionNumber;
						// Read values from data files
						FileInputStream fis = new FileInputStream(directory);
						InputStreamReader isr = new InputStreamReader(fis);
						BufferedReader br = new BufferedReader(isr);
						// System.out.println(directory);
						String aux = br.readLine();
						while (aux != null) {
							data[indicator][problem][algorithm].add(Double
									.parseDouble(aux));
							// System.out.println(Double.parseDouble(aux));
							aux = br.readLine();
						} // while
					} // for
				} // for
			} // for

			double[][][] mean;
			double[][][] median;
			double[][][] stdDeviation;
			double[][][] iqr;
			double[][][] max;
			double[][][] min;
			int[][][] numberOfValues;

			Map<String, Double> statValues = new HashMap<String, Double>();

			statValues.put("mean", 0.0);
			statValues.put("median", 0.0);
			statValues.put("stdDeviation", 0.0);
			statValues.put("iqr", 0.0);
			statValues.put("max", 0.0);
			statValues.put("min", 0.0);

			mean = new double[indicatorList_.length][][];
			median = new double[indicatorList_.length][][];
			stdDeviation = new double[indicatorList_.length][][];
			iqr = new double[indicatorList_.length][][];
			min = new double[indicatorList_.length][][];
			max = new double[indicatorList_.length][][];
			numberOfValues = new int[indicatorList_.length][][];

			for (int indicator = 0; indicator < indicatorList_.length; indicator++) {
				// A data vector per problem
				mean[indicator] = new double[problemList_.length][];
				median[indicator] = new double[problemList_.length][];
				stdDeviation[indicator] = new double[problemList_.length][];
				iqr[indicator] = new double[problemList_.length][];
				min[indicator] = new double[problemList_.length][];
				max[indicator] = new double[problemList_.length][];
				numberOfValues[indicator] = new int[problemList_.length][];

				for (int problem = 0; problem < problemList_.length; problem++) {
					mean[indicator][problem] = new double[algorithmNameList_.length];
					median[indicator][problem] = new double[algorithmNameList_.length];
					stdDeviation[indicator][problem] = new double[algorithmNameList_.length];
					iqr[indicator][problem] = new double[algorithmNameList_.length];
					min[indicator][problem] = new double[algorithmNameList_.length];
					max[indicator][problem] = new double[algorithmNameList_.length];
					numberOfValues[indicator][problem] = new int[algorithmNameList_.length];

					for (int algorithm = 0; algorithm < algorithmNameList_.length; algorithm++) {
						Collections.sort(data[indicator][problem][algorithm]);

						String directory = experimentBaseDirectory_;
						directory += "/" + algorithmNameList_[algorithm];
						directory += "/" + problemList_[problem];
						directory += "/" + indicatorList_[indicator];

						// System.out.println("----" + directory + "-----");
						// calculateStatistics(data[indicator][problem][algorithm],
						// meanV, medianV, minV, maxV, stdDeviationV, iqrV) ;
						calculateStatistics(
								data[indicator][problem][algorithm], statValues);
						/*
						 * System.out.println("Mean: " +
						 * statValues.get("mean"));
						 * System.out.println("Median : " +
						 * statValues.get("median"));
						 * System.out.println("Std : " +
						 * statValues.get("stdDeviation"));
						 * System.out.println("IQR : " + statValues.get("iqr"));
						 * System.out.println("Min : " + statValues.get("min"));
						 * System.out.println("Max : " + statValues.get("max"));
						 * System.out.println("N_values: " +
						 * data[indicator][problem][algorithm].size()) ;
						 */
						mean[indicator][problem][algorithm] = statValues
								.get("mean");
						median[indicator][problem][algorithm] = statValues
								.get("median");
						stdDeviation[indicator][problem][algorithm] = statValues
								.get("stdDeviation");
						iqr[indicator][problem][algorithm] = statValues
								.get("iqr");
						min[indicator][problem][algorithm] = statValues
								.get("min");
						max[indicator][problem][algorithm] = statValues
								.get("max");
						numberOfValues[indicator][problem][algorithm] = data[indicator][problem][algorithm]
								.size();
					}
				}
			}

		
			for (int i = 0; i < indicatorList_.length; i++) {
				for(int j=0; j<problemList_.length;j++){
					for(int z=0; z<algorithmNameList_.length;z++){
						String ckey; //column key
						if(indicatorList_[i].equals("EPSILON"))
							ckey=algorithmNameList_[z]+"."+"EPS";
						else
							ckey=algorithmNameList_[z]+"."+indicatorList_[i];
					
						String cKeyC;  //common cKey
						String rKeyC; //common row key
						
						cKeyC= ckey+"."+"M";
						rKeyC = problemList_[j]+".R"+regionNumber;
						MultiKey<? extends MultiKey<String>> mk = new MultiKey(rKeyC, cKeyC);
						table.put(mk, mean[i][j][z]);
						
						cKeyC= ckey+"."+"Md";
						rKeyC = problemList_[j]+".R"+regionNumber;
						mk = new MultiKey(rKeyC, cKeyC);
						table.put(mk, median[i][j][z]);
						
						cKeyC= ckey+"."+"St";
						rKeyC = problemList_[j]+".R"+regionNumber;
						mk = new MultiKey(rKeyC, cKeyC);
						table.put(mk, stdDeviation[i][j][z]);
						
						cKeyC= ckey+"."+"IQR";
						rKeyC = problemList_[j]+".R"+regionNumber;
						mk = new MultiKey(rKeyC, cKeyC);
						table.put(mk, iqr[i][j][z]);
						
						
					}
					
					
				}
				
				
				
				
				//printMeanStdDev(latexFile, i, mean, stdDeviation, numberOfValues);
				//printMedianIQR(latexFile, i, median, iqr, numberOfValues);

			} // for
		
		}
	} // generateLatexTables

	public static void main(String[] args) throws JMException, IOException {

		List<Double[]> regions = new ArrayList<Double[]>(3);

		regions.add(new Double[] { 0.8, 0.95 });
		regions.add(new Double[] { 0.40, 0.50 });
		regions.add(new Double[] { 0.15, 0.20 });
		
		 /*regions.add(new Double[] {0.8, 0.85});
		 regions.add(new Double[] {0.35, 0.50});
		 regions.add(new Double[] {0.15, 0.20});*/

		CSVFileGenrator exp = new CSVFileGenrator(regions); // exp = experiment

		//exp.experimentName_ = "pAGEStudyZDT";
		//exp.algorithmNameList_ = new String[] { "pAGE.F", "pAGE.O", "AGE.1", "AGE.2", "AGE.3"};
		
		//exp.experimentName_ = "pAGEStudyZDTMaxEvo";
		//exp.algorithmNameList_ = new String[] {"pAGE.F.Max", "pAGE.O.Max"};
		
		//exp.experimentName_ = "pNSGAIIStudyZDT";
		//exp.algorithmNameList_ = new String[] { "pNSGAII", "NSGAII.1", "NSGAII.2", "NSGAII.3"};
		
		
		//exp.experimentName_ = "pNSGAIIStudyZDTMaxEvo";
		//exp.algorithmNameList_ = new String[] { "pNSGAII.Max"};
		
		//exp.experimentName_ = "pSPEA2StudyZDT";
		//exp.algorithmNameList_ = new String[] { "pSPEA2", "SPEA2.1", "SPEA2.2", "SPEA2.3"};
		
		
		//exp.experimentName_ = "pSPEA2StudyZDTMaxEvo";
		//exp.algorithmNameList_ = new String[] { "pSPEA2.Max"};
		
		//DTLZ
		//exp.experimentName_ = "pAGEStudyDTLZ";
		//exp.algorithmNameList_ = new String[] { "pAGE.F.1", "pAGE.F.2", "pAGE.O.1", "pAGE.O.2", "AGE.1", "AGE.2"};
		
		//exp.experimentName_ = "pNSGAIIStudyDTLZ";
		//exp.algorithmNameList_ = new String[] { "pNSGAII.1", "pNSGAII.2", "NSGAII.1", "NSGAII.2"};
		
		exp.experimentName_ = "pSPEA2StudyDTLZ";
		exp.algorithmNameList_ = new String[] { "pSPEA2.1", "pSPEA2.2", "SPEA2.1", "SPEA2.2"};
				
		
		//exp.problemList_ = new String[] { "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6" };
		exp.problemList_ = new String [] {"DTLZ2", "DTLZ3"};
		exp.indicatorList_ = new String[] { "HV", "IGD", "EPSILON" };

		/*exp.problemList_ = new String[] { "DTLZ2" };
		exp.paretoFrontFile_ = new String[] { "DTLZ2.3D.pf"};
		exp.indicatorList_ = new String[] { "HV",  "IGD", "EPSILON" };*/

		
		int numberOfAlgorithms = exp.algorithmNameList_.length;

		exp.experimentBaseDirectory_ = "C:\\Users\\mahbub\\Documents\\GitHub\\pMOEA jMetal4.5\\experiments\\"
				+ exp.experimentName_;
		
		exp.readFromCSVFile("resultsDTLZ.CSV");
			
		
		// Generate latex tables (comment this sentence is not desired)
		exp.generateLatexTables();
		
		exp.writeToFCSVFile("resultsDTLZ.CSV");

	
		
	} // main

	void readFromCSVFile(String fileName) throws IOException{
		//Delimiter used in CSV file
		final String COMMA_DELIMITER = ",";

		BufferedReader fileReader = null;

		try{
			fileReader = new BufferedReader(new FileReader(fileName));
		}catch(FileNotFoundException e){
			return;
		}
		//Read the CSV file header to skip it
         fileReader.readLine();
         String line= "";
         
         while ((line = fileReader.readLine()) != null) {
        	 String[] tokens = line.split(COMMA_DELIMITER);
        	 if (tokens.length > 1) {
        		 for(int j=1;j<tokens.length;j++){
						MultiKey<? extends MultiKey<String>> mk = new MultiKey(tokens[0], row[j-1]);
						try{
							table.put(mk, Double.parseDouble(tokens[j]));
						}catch(NumberFormatException e){
							System.out.println("exception");
						}
        		 }
        	 }
         }

         fileReader.close();
	}
	
	void writeToFCSVFile(String fileName) throws IOException{
		//Delimiter used in CSV file
		final String COMMA_DELIMITER = ",";
		final String NEW_LINE_SEPARATOR = "\n";

		FileWriter fileWriter = null;

		fileWriter = new FileWriter(fileName);
		//Write the CSV file header
		fileWriter.append("problem");
		fileWriter.append(COMMA_DELIMITER);
		for(int i=0;i<row.length;i++){ 
			fileWriter.append(row[i].toString());
			if(i<row.length-1)
				fileWriter.append(COMMA_DELIMITER);
		}
       //Add a new line separator after the header
       fileWriter.append(NEW_LINE_SEPARATOR);
       
       for(int j=0;j<col.length;j++){
    	   fileWriter.append(col[j]);
    	   fileWriter.append(COMMA_DELIMITER);
    	   for(int i=0;i<row.length;i++){
    		   MultiKey<? extends MultiKey<String>> mk = new MultiKey(col[j], row[i]);
    		   Double value = table.get(mk);
    		   if(value == null)
    			   fileWriter.append("");
    		   else
    			   fileWriter.append(value.toString());
    		   if(i<row.length-1)
    			   fileWriter.append(COMMA_DELIMITER);
    	   }
    	   fileWriter.append(NEW_LINE_SEPARATOR);
       }
       fileWriter.close();

	}

	@Override
	public void algorithmSettings(String problemName, int problemId,
			Algorithm[] algorithm) throws ClassNotFoundException {
		// TODO Auto-generated method stub
		
	}
} // NSGAIIStudy

