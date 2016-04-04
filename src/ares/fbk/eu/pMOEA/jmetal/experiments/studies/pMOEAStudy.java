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
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.AGE_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.NSGAII_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.pAGE_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.pNSGAII_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.settings.pSPEA2_Settings;
import ares.fbk.eu.pMOEA.jmetal.experiments.util.Friedman;
import ares.fbk.eu.pMOEA.jmetal.experiments.util.Statistics;
import ares.fbk.eu.pMOEA.jmetal.metaheuristics.age.pAGE;
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

import org.omg.CORBA.portable.IndirectionException;

/**
 * Class implementing an example of experiment using NSGA-II as base algorithm.
 * The experiment consisting in studying the effect of the crossover probability
 * in NSGA-II.
 */
public class pMOEAStudy extends Experiment {

	public List<Double[]> regions = new ArrayList<Double[]>(3);

	public pMOEAStudy(List<Double[]> regions) {
		this.regions = regions;
	}

	/**
	 * Configures the algorithms in each independent run
	 * 
	 * @param problemName
	 *            The problem to solve
	 * @param problemIndex
	 * @param algorithm
	 *            Array containing the algorithms to run
	 * @throws ClassNotFoundException
	 */
	public synchronized void algorithmSettings(String problemName,
			int problemIndex, Algorithm[] algorithm)
			throws ClassNotFoundException {
		try {

			int numberOfAlgorithms = algorithmNameList_.length;

			HashMap[] parameters = new HashMap[numberOfAlgorithms];

			for (int i = 0; i < numberOfAlgorithms; i++) {
				parameters[i] = new HashMap();
			} // for

			if (!paretoFrontFile_[problemIndex].equals("")) {
				for (int i = 0; i < numberOfAlgorithms; i++)
					parameters[i].put("paretoFrontFile_",
							paretoFrontFile_[problemIndex]);
			} // if

			if ((!paretoFrontFile_[problemIndex].equals(""))
					|| (paretoFrontFile_[problemIndex] == null)) {
				for (int i = 0; i < numberOfAlgorithms; i++)
					parameters[i].put("paretoFrontFile_",
							paretoFrontFile_[problemIndex]);
			} // if

			// pNSGAII configuration
			parameters[0].put("populationSize_", 30);
			parameters[0].put("maxEvaluations_", 24000);

			// pSPEA2 configuration
			parameters[1].put("populationSize_", 30);
			parameters[1].put("archiveSize", 30);
			parameters[1].put("maxEvaluations_", 24000);
			
			// AGE.F (AGE offline)  configuration
			parameters[2].put("populationSize_", 30);
			parameters[2].put("epsilonGridWidth", 0.0);
			parameters[2].put("maxEvaluations_", 24000);
			
			// AGE.O (AGE online)  configuration
			parameters[3].put("populationSize_", 30);
			parameters[3].put("epsilonGridWidth", 0.01);
			parameters[3].put("maxEvaluations_", 24000);
						
 
			// 0 -> pNSGAII, 1 -> NSGAII
			algorithm[0] = new pNSGAII_Settings(problemName)
					.configure(parameters[0]);
			algorithm[1] = new pSPEA2_Settings(problemName)
					.configure(parameters[1]);
			algorithm[2] = new pAGE_Settings(problemName)
					.configure(parameters[2]);
			algorithm[3] = new pAGE_Settings(problemName)
					.configure(parameters[3]);
	
			/*algorithm[2] = new AGE_Settings(problemName)
			.configure(parameters[2]);
			algorithm[3] = new AGE_Settings(problemName)
			.configure(parameters[3]);*/

			
		} catch (IllegalArgumentException ex) {
			Logger.getLogger(pMOEAStudy.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (IllegalAccessException ex) {
			Logger.getLogger(pMOEAStudy.class.getName()).log(Level.SEVERE,
					null, ex);
		} catch (JMException ex) {
			Logger.getLogger(pMOEAStudy.class.getName()).log(Level.SEVERE,
					null, ex);
		}
	} // algorithmSettings

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

	@Override
	public void generateQualityIndicators() {

		checkParetoFronts();

		if (indicatorList_.length > 0) {

			for (int algorithmIndex = 0; algorithmIndex < algorithmNameList_.length; algorithmIndex++) {

				String algorithmDirectory;
				algorithmDirectory = experimentBaseDirectory_ + "/data/"
						+ algorithmNameList_[algorithmIndex] + "/";

				for (int problemIndex = 0; problemIndex < problemList_.length; problemIndex++) {

					String problemDirectory = algorithmDirectory
							+ problemList_[problemIndex];
					String paretoFrontPath = frontPath_[problemIndex];

					for (String anIndicatorList_ : indicatorList_) {
						System.out.println("Experiment - Quality indicator: "
								+ anIndicatorList_);

						resetFile(problemDirectory + "/" + anIndicatorList_);

						for (int numRun = 0; numRun < independentRuns_; numRun++) {

							String outputParetoFrontFilePath;
							outputParetoFrontFilePath = problemDirectory
									+ "/FUN." + numRun;
							String solutionFrontFile = outputParetoFrontFilePath;
							String qualityIndicatorFile = problemDirectory;
							String qualityIndicatorFile__ = "";
							double value = 0;

							double[][] trueFront = new Hypervolume().utils_
									.readFront(paretoFrontPath);

							double[][] solutionFront = new Hypervolume().utils_
									.readFront(solutionFrontFile);

							// declaration: regional true Pretofront
							List<double[]>[] rTrueFronts = new List[regions
									.size()];

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

							// identify regional solution front
							for (int i = 0; i < solutionFront.length; i++) {
								for (int regionNumber = 0; regionNumber < regions
										.size(); regionNumber++) {
									if (solutionFront[i][0] >= regions
											.get(regionNumber)[0]
											&& solutionFront[i][0] <= regions
													.get(regionNumber)[1]) {
										double[] sol = solutionFront[i];
										rFronts[regionNumber].add(sol);

									}
								}
							}

							for (int i = 0; i < regions.size(); i++) {
								
								if (anIndicatorList_.equals("HV")) {
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/HV" + ".r." + i;
								}else if (anIndicatorList_.equals("SPREAD")) {
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/SPREAD" + ".r." + i;
								}else if (anIndicatorList_.equals("IGD")) {
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/IGD" + ".r." + i;
								}else if (anIndicatorList_.equals("EPSILON")) {
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/EPSILON" + ".r." + i;
								}
								
								if(rTrueFronts[i].size() != 0){

								trueFront = new double[rTrueFronts[i].size()][];
								solutionFront = new double[rFronts[i].size()][];

								for (int j = 0; j < rTrueFronts[i].size(); j++) {
									trueFront[j] = rTrueFronts[i].get(j);
								}

								for (int j = 0; j < rFronts[i].size(); j++) {
									solutionFront[j] = rFronts[i].get(j);
								}

								if (anIndicatorList_.equals("HV")) {

									// when rFront.sizr() is equal to zero,
									// there is no need to calculate indicators'
									// values
									if (rFronts[i].size() != 0) {

										pHypervolume indicators = new pHypervolume();
										// Hypervolume indicators = new
										// Hypervolume();
										value = indicators.hypervolume(
												solutionFront, trueFront,
												trueFront[0].length,
												regions.get(i)[1]);

										/*
										 * value = indicators.hypervolume(
										 * solutionFront, trueFront,
										 * trueFront[0].length);
										 */
									}
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/HV" + ".r." + i;

								}

								if (anIndicatorList_.equals("SPREAD")) {
									if (trueFront[0].length == 2) {
										if (rFronts[i].size() != 0) {
											Spread indicators = new Spread();

											// double[][] solutionFront =
											// indicators.utils_.readFront(solutionFrontFile);
											// double[][] trueFront =
											// indicators.utils_.readFront(paretoFrontPath);
											value = indicators.spread(
													solutionFront, trueFront,
													trueFront[0].length);
										}

										qualityIndicatorFile__ = qualityIndicatorFile
												+ "/SPREAD" + ".r." + i;
									} else {
										if (rFronts[i].size() != 0) {
											GeneralizedSpread indicators = new GeneralizedSpread();
											// double[][] solutionFront =
											// indicators.utils_
											// .readFront(solutionFrontFile);
											value = indicators
													.generalizedSpread(
															solutionFront,
															trueFront,
															trueFront[0].length);
										}
										qualityIndicatorFile__ = qualityIndicatorFile
												+ "/SPREAD" + ".r." + i;

									}
								}

								if (anIndicatorList_.equals("IGD")) {
									if (rFronts[i].size() != 0) {
										InvertedGenerationalDistance indicators = new InvertedGenerationalDistance();
										// double[][] solutionFront =
										// indicators.utils_
										// .readFront(solutionFrontFile);
										// double[][] trueFront =
										// indicators.utils_.readFront(paretoFrontPath);
										value = indicators
												.invertedGenerationalDistance(
														solutionFront,
														trueFront,
														trueFront[0].length);
									}
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/IGD" + ".r." + i;
								}
								if (anIndicatorList_.equals("EPSILON")) {
									if (rFronts[i].size() != 0) {
										Epsilon indicators = new Epsilon();
										// double[][] solutionFront =
										// indicators.utils_
										// .readFront(solutionFrontFile);
										// double[][] trueFront =
										// indicators.utils_.readFront(paretoFrontPath);
										value = indicators.epsilon(
												solutionFront, trueFront,
												trueFront[0].length);
									}
									qualityIndicatorFile__ = qualityIndicatorFile
											+ "/EPSILON" + ".r." + i;
								}
								}

								if (!qualityIndicatorFile__
										.equals(problemDirectory)) {
									FileWriter os;
									try {
										os = new FileWriter(
												qualityIndicatorFile__, true);
										if (rFronts[i].size() != 0) {
											os.write("" + value + "\n");
										}
										os.close();
									} catch (IOException ex) {
										Logger.getLogger(
												Experiment.class.getName())
												.log(Level.SEVERE, null, ex);
									}
								} // if
							} // for
						}
					} // for
				} // for
			} // for
		} // if
	}

	
	
	 protected void printMeanStdDev(String fileName, int indicator, double[][][] mean, double[][][] stdDev, int[][][] numberOFValues) throws IOException {
		    FileWriter os = new FileWriter(fileName, true);
		    os.write("\\" + "\n");
		    os.write("\\begin{table}" + "\n");
		    os.write("\\caption{" + indicatorList_[indicator] + ". Mean and standard deviation}" + "\n");
		    os.write("\\label{table:mean." + indicatorList_[indicator] + "}" + "\n");
		    os.write("\\centering" + "\n");
		    os.write("\\begin{scriptsize}" + "\n");
		    os.write("\\begin{tabular}{l");

		    // calculate the number of columns
		    for (String anAlgorithmNameList_ : algorithmNameList_) {
		      os.write("l");
		    }
		    os.write("}\n");

		    os.write("\\hline");
		    // write table head
		    for (int i = -1; i < algorithmNameList_.length; i++) {
		      if (i == -1) {
		        os.write(" & ");
		      } else if (i == (algorithmNameList_.length - 1)) {
		        os.write(" " + algorithmNameList_[i] + "\\\\" + "\n");
		      } else {
		        os.write("" + algorithmNameList_[i] + " & ");
		      }
		    }
		    os.write("\\hline" + "\n");

		    String m, s;
		    // write lines
		    for (int i = 0; i < problemList_.length; i++) {
		      // find the best value and second best value
		      double bestValue;
		      double bestValueIQR;
		      double secondBestValue;
		      double secondBestValueIQR;
		      int bestIndex = -1;
		      int secondBestIndex = -1;
		      if ((Boolean) indicatorMinimize_.get(indicatorList_[indicator])) {// minimize by default
		        bestValue = Double.MAX_VALUE;
		        bestValueIQR = Double.MAX_VALUE;
		        secondBestValue = Double.MAX_VALUE;
		        secondBestValueIQR = Double.MAX_VALUE;
		        for (int j = 0; j < (algorithmNameList_.length); j++) {
		          if ((mean[indicator][i][j] < bestValue) ||
		                  ((mean[indicator][i][j] == bestValue) && (stdDev[indicator][i][j] < bestValueIQR))) {
		            secondBestIndex = bestIndex;
		            secondBestValue = bestValue;
		            secondBestValueIQR = bestValueIQR;
		            bestValue = mean[indicator][i][j];
		            bestValueIQR = stdDev[indicator][i][j];
		            bestIndex = j;
		          } else if ((mean[indicator][i][j] < secondBestValue) ||
		                  ((mean[indicator][i][j] == secondBestValue) && (stdDev[indicator][i][j] < secondBestValueIQR))) {
		            secondBestIndex = j;
		            secondBestValue = mean[indicator][i][j];
		            secondBestValueIQR = stdDev[indicator][i][j];
		          } // else if
		        }
		      } // if
		      else { // indicator to maximize e.g., the HV
		        bestValue = Double.MIN_VALUE;
		        bestValueIQR = Double.MIN_VALUE;
		        secondBestValue = Double.MIN_VALUE;
		        secondBestValueIQR = Double.MIN_VALUE;
		        for (int j = 0; j < (algorithmNameList_.length); j++) {
		          if ((mean[indicator][i][j] > bestValue) ||
		                  ((mean[indicator][i][j] == bestValue) && (stdDev[indicator][i][j] < bestValueIQR))) {
		            secondBestIndex = bestIndex;
		            secondBestValue = bestValue;
		            secondBestValueIQR = bestValueIQR;
		            bestValue = mean[indicator][i][j];
		            bestValueIQR = stdDev[indicator][i][j];
		            bestIndex = j;
		          } else if ((mean[indicator][i][j] > secondBestValue) ||
		                  ((mean[indicator][i][j] == secondBestValue) && (stdDev[indicator][i][j] < secondBestValueIQR))) {
		            secondBestIndex = j;
		            secondBestValue = mean[indicator][i][j];
		            secondBestValueIQR = stdDev[indicator][i][j];
		          } // else if
		        } // for
		      } // else

		      os.write(problemList_[i].replace("_", "\\_") + " & ");
		      for (int j = 0; j < (algorithmNameList_.length - 1); j++) {
		        if (j == bestIndex) {
		          os.write("\\cellcolor{gray95}");
		        }
		        if (j == secondBestIndex) {
		          os.write("\\cellcolor{gray25}");
		        }

		        m = String.format(Locale.ENGLISH, "%10.2e", mean[indicator][i][j]);
		        s = String.format(Locale.ENGLISH, "%8.1e", stdDev[indicator][i][j]);
		        os.write("$" + m + "_{" + s + "}$" +"("+numberOFValues[indicator][i][j]+")" +" & ");
		      }
		      if (bestIndex == (algorithmNameList_.length - 1)) {
		        os.write("\\cellcolor{gray95}");
		      }
		      //shahriar
		      if (secondBestIndex == (algorithmNameList_.length - 1)) {
			        os.write("\\cellcolor{gray25}");
			      }
		      m = String.format(Locale.ENGLISH, "%10.2e", mean[indicator][i][algorithmNameList_.length - 1]);
		      s = String.format(Locale.ENGLISH, "%8.1e", stdDev[indicator][i][algorithmNameList_.length - 1]);
		      os.write("$" + m + "_{" + s + "}$" +"("+numberOFValues[indicator][i][algorithmNameList_.length - 1]+")" +"\\\\" + "\n");
		    } // for
		    //os.write("" + mean[0][problemList_.length-1][algorithmNameList_.length-1] + "\\\\"+ "\n" ) ;

		    os.write("\\hline" + "\n");
		    os.write("\\end{tabular}" + "\n");
		    os.write("\\end{scriptsize}" + "\n");
		    os.write("\\end{table}" + "\n");
		    os.close();
		  } // printMeanStdDev

		  protected void printMedianIQR(String fileName, int indicator, double[][][] median, double[][][] IQR, int [][][] numberOfValues) throws IOException {
		    FileWriter os = new FileWriter(fileName, true);
		    os.write("\\" + "\n");
		    os.write("\\begin{table}" + "\n");
		    os.write("\\caption{" + indicatorList_[indicator] + ". Median and IQR}" + "\n");
		    os.write("\\label{table:median." + indicatorList_[indicator] + "}" + "\n");
		    os.write("\\begin{scriptsize}" + "\n");
		    os.write("\\centering" + "\n");
		    os.write("\\begin{tabular}{l");

		    // calculate the number of columns
		    for (String anAlgorithmNameList_ : algorithmNameList_) {
		      os.write("l");
		    }
		    os.write("}\n");

		    os.write("\\hline");
		    // write table head
		    for (int i = -1; i < algorithmNameList_.length; i++) {
		      if (i == -1) {
		        os.write(" & ");
		      } else if (i == (algorithmNameList_.length - 1)) {
		        os.write(" " + algorithmNameList_[i] + "\\\\" + "\n");
		      } else {
		        os.write("" + algorithmNameList_[i] + " & ");
		      }
		    }
		    os.write("\\hline" + "\n");

		    String m, s;
		    // write lines
		    for (int i = 0; i < problemList_.length; i++) {
		      // find the best value and second best value
		      double bestValue;
		      double bestValueIQR;
		      double secondBestValue;
		      double secondBestValueIQR;
		      int bestIndex = -1;
		      int secondBestIndex = -1;
		      if ((Boolean) indicatorMinimize_.get(indicatorList_[indicator])) {// minimize by default
		        bestValue = Double.MAX_VALUE;
		        bestValueIQR = Double.MAX_VALUE;
		        secondBestValue = Double.MAX_VALUE;
		        secondBestValueIQR = Double.MAX_VALUE;
		        for (int j = 0; j < (algorithmNameList_.length); j++) {
		          if ((median[indicator][i][j] < bestValue) ||
		                  ((median[indicator][i][j] == bestValue) && (IQR[indicator][i][j] < bestValueIQR))) {
		            secondBestIndex = bestIndex;
		            secondBestValue = bestValue;
		            secondBestValueIQR = bestValueIQR;
		            bestValue = median[indicator][i][j];
		            bestValueIQR = IQR[indicator][i][j];
		            bestIndex = j;
		          } else if ((median[indicator][i][j] < secondBestValue) ||
		                  ((median[indicator][i][j] == secondBestValue) && (IQR[indicator][i][j] < secondBestValueIQR))) {
		            secondBestIndex = j;
		            secondBestValue = median[indicator][i][j];
		            secondBestValueIQR = IQR[indicator][i][j];
		          } // else if
		        } // for
		      } // if
		      else { // indicator to maximize e.g., the HV
		        bestValue = Double.MIN_VALUE;
		        bestValueIQR = Double.MIN_VALUE;
		        secondBestValue = Double.MIN_VALUE;
		        secondBestValueIQR = Double.MIN_VALUE;
		        for (int j = 0; j < (algorithmNameList_.length); j++) {
		          if ((median[indicator][i][j] > bestValue) ||
		                  ((median[indicator][i][j] == bestValue) && (IQR[indicator][i][j] < bestValueIQR))) {
		            secondBestIndex = bestIndex;
		            secondBestValue = bestValue;
		            secondBestValueIQR = bestValueIQR;
		            bestValue = median[indicator][i][j];
		            bestValueIQR = IQR[indicator][i][j];
		            bestIndex = j;
		          } else if ((median[indicator][i][j] > secondBestValue) ||
		                  ((median[indicator][i][j] == secondBestValue) && (IQR[indicator][i][j] < secondBestValueIQR))) {
		            secondBestIndex = j;
		            secondBestValue = median[indicator][i][j];
		            secondBestValueIQR = IQR[indicator][i][j];
		          } // else if
		        } // for
		      } // else

		      os.write(problemList_[i].replace("_", "\\_") + " & ");
		      for (int j = 0; j < (algorithmNameList_.length - 1); j++) {
		        if (j == bestIndex) {
		          os.write("\\cellcolor{gray95}");
		        }
		        if (j == secondBestIndex) {
		          os.write("\\cellcolor{gray25}");
		        }
		        m = String.format(Locale.ENGLISH, "%10.2e", median[indicator][i][j]);
		        s = String.format(Locale.ENGLISH, "%8.1e", IQR[indicator][i][j]);
		        os.write("$" + m + "_{" + s + "}$ "+"("+numberOfValues[indicator][i][j]+")" +" & ");
		      }
		      if (bestIndex == (algorithmNameList_.length - 1)) {
		        os.write("\\cellcolor{gray95}");
		      }
		      //shahriar
		      if (secondBestIndex == (algorithmNameList_.length - 1)) {
			        os.write("\\cellcolor{gray25}");
			      }
		      m = String.format(Locale.ENGLISH, "%10.2e", median[indicator][i][algorithmNameList_.length - 1]);
		      s = String.format(Locale.ENGLISH, "%8.1e", IQR[indicator][i][algorithmNameList_.length - 1]);
		      os.write("$" + m + "_{" + s + "}$" +"("+numberOfValues[indicator][i][algorithmNameList_.length - 1]+")" +" \\\\" + "\n");
		    } // for
		    //os.write("" + mean[0][problemList_.length-1][algorithmNameList_.length-1] + "\\\\"+ "\n" ) ;

		    os.write("\\hline" + "\n");
		    os.write("\\end{tabular}" + "\n");
		    os.write("\\end{scriptsize}" + "\n");
		    os.write("\\end{table}" + "\n");
		    os.close();
		  } // printMedianIQR

	public void generateLatexTables() throws FileNotFoundException, IOException {
		latexDirectory_ = experimentBaseDirectory_ + "/" + latexDirectory_;
		System.out.println("latex directory: " + latexDirectory_);

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

			File latexOutput;
			latexOutput = new File(latexDirectory_);
			if (!latexOutput.exists()) {
				boolean result = new File(latexDirectory_).mkdirs();
				System.out
						.println("Creating " + latexDirectory_ + " directory");
			}
			// System.out.println("Experiment name: " + experimentName_);
			String latexFile = latexDirectory_ + "/" + experimentName_ + ".tex";
			if (regionNumber == 0) {
				printHeaderLatexCommands(latexFile);
			}
			// print regions
			printBoldText(latexFile, "Region: " + regionNumber);
			for (int i = 0; i < indicatorList_.length; i++) {
				printMeanStdDev(latexFile, i, mean, stdDeviation, numberOfValues);
				printMedianIQR(latexFile, i, median, iqr, numberOfValues);

			} // for
			printClearPage(latexFile);
			printClearPage(latexFile);

			if (regionNumber == regions.size() - 1) {
				printEndLatexCommands(latexFile);
			}
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

		pMOEAStudy exp = new pMOEAStudy(regions); // exp = experiment

		exp.experimentName_ = "pMOEA.NSGA.SPEA.AGE.MaxEvo";
		exp.algorithmNameList_ = new String[] { "pNSGAII.Max","pSPEA2.Max","pAGE.F.Max", "pAGE.O.Max"};

		exp.problemList_ = new String[] { "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6" };
		exp.paretoFrontFile_ = new String[] { "ZDT1.pf", "ZDT2.pf", "ZDT3.pf","ZDT4.pf", "ZDT6.pf" };
		exp.indicatorList_ = new String[] { "HV", "IGD", "EPSILON" };

		//exp.problemList_ = new String[] { "DTLZ2" };
		//exp.paretoFrontFile_ = new String[] { "DTLZ2.3D.pf"};
		//exp.indicatorList_ = new String[] { "HV",  "IGD", "SPREAD","EPSILON" };

		
		int numberOfAlgorithms = exp.algorithmNameList_.length;

		exp.experimentBaseDirectory_ = "C:\\\\Users\\\\mahbub\\\\Documents\\\\GitHub\\\\pMOEA jMetal4.5\\\\experiments\\\\"
				+ exp.experimentName_;
		exp.paretoFrontDirectory_ = "C:\\Users\\mahbub\\Documents\\GitHub\\EnergyPLANDomainKnowledgeEAStep1\\paretoFronts\\";

		exp.algorithmSettings_ = new Settings[numberOfAlgorithms];

		exp.independentRuns_ = 100;

		exp.initExperiment();

		// Run the experiments
		int numberOfThreads;
		exp.runExperiment(numberOfThreads = 4);

		exp.generateQualityIndicators();

		// Generate latex tables (comment this sentence is not desired)
		exp.generateLatexTables();

		// Configure the R scripts to be generated
		int rows;
		int columns;
		String prefix;
		String[] problems;

		rows = 2;
		columns = 3;
		prefix = new String("Problems");
		problems = new String[] { "ZDT1", "ZDT2", "ZDT3", "ZDT4", "ZDT6" };
		//problems = new String[] { "DTLZ2"};

		String[] modifiedIndicatorList = new String[exp.indicatorList_.length
				* regions.size()];

		int z = 0;
		for (int i = 0; i < regions.size(); i++) {
			for (int j = 0; j < exp.indicatorList_.length; j++) {
				modifiedIndicatorList[z++] = exp.indicatorList_[j] + ".r." + i;
			}
		}
		exp.indicatorList_ = modifiedIndicatorList;

		boolean notch;
		//exp.generateRBoxplotScripts(rows, columns, problems, prefix,
			//	notch = true, exp);
		
		exp.generateRWilcoxonScripts(problems, prefix, exp);

		// Applying Friedman test
		// Friedman test = new Friedman(exp);
		 //test.executeTest("EPSILON");
		// test.executeTest("HV");
		 //test.executeTest("SPREAD");
	} // main
} // NSGAIIStudy

