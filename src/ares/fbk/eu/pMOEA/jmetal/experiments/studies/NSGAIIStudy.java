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

package ares.fbk.eu.pMOEA.jmetal
.experiments.studies;

import ares.fbk.eu.pMOEA.jmetal
.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal
.experiments.Experiment;
import ares.fbk.eu.pMOEA.jmetal
.experiments.Settings;
import ares.fbk.eu.pMOEA.jmetal
.experiments.settings.NSGAII_Settings;
import ares.fbk.eu.pMOEA.jmetal
.experiments.util.Friedman;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Epsilon;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.InvertedGenerationalDistance;
import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Spread;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;

import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class implementing an example of experiment using NSGA-II as base algorithm.
 * The experiment consisting in studying the effect of the crossover probability
 * in NSGA-II.
 */
public class NSGAIIStudy extends Experiment {
  /**
   * Configures the algorithms in each independent run
   * @param problemName The problem to solve
   * @param problemIndex
   * @param algorithm Array containing the algorithms to run
   * @throws ClassNotFoundException 
   */
  public synchronized void algorithmSettings(String problemName, 
  		                                       int problemIndex, 
  		                                       Algorithm[] algorithm) 
    throws ClassNotFoundException {  	
  	try {
      int numberOfAlgorithms = algorithmNameList_.length;

      HashMap[] parameters = new HashMap[numberOfAlgorithms];

      for (int i = 0; i < numberOfAlgorithms; i++) {
        parameters[i] = new HashMap();
      } // for

      if (!paretoFrontFile_[problemIndex].equals("")) {
        for (int i = 0; i < numberOfAlgorithms; i++)
          parameters[i].put("paretoFrontFile_", paretoFrontFile_[problemIndex]);
      } // if

      parameters[0].put("crossoverProbability_", 1.0);
      parameters[1].put("crossoverProbability_", 0.9);
      parameters[2].put("crossoverProbability_", 0.8);
      parameters[3].put("crossoverProbability_", 0.7); 

      if ((!paretoFrontFile_[problemIndex].equals("")) || 
      		(paretoFrontFile_[problemIndex] == null)) {
        for (int i = 0; i < numberOfAlgorithms; i++)
          parameters[i].put("paretoFrontFile_",  paretoFrontFile_[problemIndex]);
      } // if
 
      for (int i = 0; i < numberOfAlgorithms; i++)
        algorithm[i] = new NSGAII_Settings(problemName).configure(parameters[i]);
      
    } catch (IllegalArgumentException ex) {
      Logger.getLogger(NSGAIIStudy.class.getName()).log(Level.SEVERE, null, ex);
    } catch (IllegalAccessException ex) {
      Logger.getLogger(NSGAIIStudy.class.getName()).log(Level.SEVERE, null, ex);
    } catch (JMException ex) {
      Logger.getLogger(NSGAIIStudy.class.getName()).log(Level.SEVERE, null, ex);
    }
  } // algorithmSettings
  
  
  
  /**
   * Generate the Quality Indicators
   */
  @Override
  public void generateQualityIndicators() {

    checkParetoFronts();

    if (indicatorList_.length > 0) {

      for (int algorithmIndex = 0; algorithmIndex < algorithmNameList_.length; algorithmIndex++) {

        String algorithmDirectory;
        algorithmDirectory = experimentBaseDirectory_ + "/data/" + algorithmNameList_[algorithmIndex] + "/";

        for (int problemIndex = 0; problemIndex < problemList_.length; problemIndex++) {

          String problemDirectory = algorithmDirectory + problemList_[problemIndex];
          String paretoFrontPath = frontPath_[problemIndex];

          for (String anIndicatorList_ : indicatorList_) {
            System.out.println("Experiment - Quality indicator: " + anIndicatorList_);

            resetFile(problemDirectory + "/" + anIndicatorList_);

            for (int numRun = 0; numRun < independentRuns_; numRun++) {

              String outputParetoFrontFilePath;
              outputParetoFrontFilePath = problemDirectory + "/FUN." + numRun;
              String solutionFrontFile = outputParetoFrontFilePath;
              String qualityIndicatorFile = problemDirectory;
              double value = 0;

              double[][] trueFront =  new Hypervolume().utils_.readFront(paretoFrontPath);

              if (anIndicatorList_.equals("HV")) {

                Hypervolume indicators = new Hypervolume();
                double[][] solutionFront =
                        indicators.utils_.readFront(solutionFrontFile);
                //double[][] trueFront =
                //        indicators.utils_.readFront(paretoFrontPath);
                value = indicators.hypervolume(solutionFront, trueFront, trueFront[0].length);

                qualityIndicatorFile = qualityIndicatorFile + "/HV";

              }
              if (anIndicatorList_.equals("SPREAD")) {
                Spread indicators = new Spread();
                double[][] solutionFront =
                        indicators.utils_.readFront(solutionFrontFile);
                //double[][] trueFront =
                //        indicators.utils_.readFront(paretoFrontPath);
                value = indicators.spread(solutionFront, trueFront, trueFront[0].length);

                qualityIndicatorFile = qualityIndicatorFile + "/SPREAD";
              }
              if (anIndicatorList_.equals("IGD")) {
                InvertedGenerationalDistance indicators = new InvertedGenerationalDistance();
                double[][] solutionFront =
                        indicators.utils_.readFront(solutionFrontFile);
                //double[][] trueFront =
                //        indicators.utils_.readFront(paretoFrontPath);
                value = indicators.invertedGenerationalDistance(solutionFront, trueFront, trueFront[0].length);

                qualityIndicatorFile = qualityIndicatorFile + "/IGD";
              }
              if (anIndicatorList_.equals("EPSILON")) {
                Epsilon indicators = new Epsilon();
                double[][] solutionFront =
                        indicators.utils_.readFront(solutionFrontFile);
                //double[][] trueFront =
                //        indicators.utils_.readFront(paretoFrontPath);
                value = indicators.epsilon(solutionFront, trueFront, trueFront[0].length);

                qualityIndicatorFile = qualityIndicatorFile + "/EPSILON";
              }


              if (!qualityIndicatorFile.equals(problemDirectory)) {
                FileWriter os;
                try {
                  os = new FileWriter(qualityIndicatorFile, true);
                  os.write("" + value + "\n");
                  os.close();
                } catch (IOException ex) {
                  Logger.getLogger(Experiment.class.getName()).log(Level.SEVERE, null, ex);
                }
              } // if
            } // for
          } // for
        } // for
      } // for
    } // if
  }
  
  public static void main(String[] args) throws JMException, IOException {
    NSGAIIStudy exp = new NSGAIIStudy() ; // exp = experiment
    
    exp.experimentName_  = "NSGAIIStudy" ;
    exp.algorithmNameList_   = new String[] {
      "NSGAIIa", "NSGAIIb", "NSGAIIc", "NSGAIId"} ;
    exp.problemList_     = new String[] {
      "ZDT1", "ZDT2", "ZDT3", "ZDT4", "DTLZ1", "WFG2"} ;
    exp.paretoFrontFile_ = new String[] {
      "ZDT1.pf", "ZDT2.pf", "ZDT3.pf","ZDT4.pf", "DTLZ1.2D.pf", "WFG2.2D.pf"} ;
    exp.indicatorList_   = new String[] {"HV", "SPREAD", "IGD", "EPSILON"} ;
    
    int numberOfAlgorithms = exp.algorithmNameList_.length ;

    exp.experimentBaseDirectory_ = "/Users/antelverde/Softw/pruebas/ares.fbk.eu.pMOEA.jmetal/" +
                                   exp.experimentName_;
    exp.paretoFrontDirectory_ = "/Users/antelverde/Softw/pruebas/data/paretoFronts";
    
    exp.algorithmSettings_ = new Settings[numberOfAlgorithms] ;
    
    exp.independentRuns_ = 30 ;

    exp.initExperiment();

    // Run the experiments
    int numberOfThreads ;
    exp.runExperiment(numberOfThreads = 6) ;

    exp.generateQualityIndicators() ;
    
    // Generate latex tables (comment this sentence is not desired)
    exp.generateLatexTables() ;
    
    // Configure the R scripts to be generated
    int rows  ;
    int columns  ;
    String prefix ;
    String [] problems ;

    rows = 2 ;
    columns = 3 ;
    prefix = new String("Problems");
    problems = new String[]{"ZDT1", "ZDT2","ZDT3", "ZDT4", "DTLZ1", "WFG2"} ;

    boolean notch ;
    exp.generateRBoxplotScripts(rows, columns, problems, prefix, notch = true, exp) ;
    exp.generateRWilcoxonScripts(problems, prefix, exp) ;

    // Applying Friedman test
    Friedman test = new Friedman(exp);
    test.executeTest("EPSILON");
    test.executeTest("HV");
    test.executeTest("SPREAD");
  } // main
} // NSGAIIStudy


