//  CMAES_main.java
//
//  Author:
//       Esteban López-Camacho <esteban@lcc.uma.es>
//
//  Copyright (c) 2013 Esteban López-Camacho
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
.metaheuristics.singleObjective.cmaes;

import ares.fbk.eu.pMOEA.jmetal
.core.Algorithm;
import ares.fbk.eu.pMOEA.jmetal
.core.Problem;
import ares.fbk.eu.pMOEA.jmetal
.core.SolutionSet;
import ares.fbk.eu.pMOEA.jmetal
.problems.singleObjective.Rosenbrock;
import ares.fbk.eu.pMOEA.jmetal
.util.JMException;

/**
 * This class runs a single-objective CMA-ES algorithm.
 */
public class CMAES_main {

  public static void main(String [] args) throws JMException, ClassNotFoundException {
    int numberOfVariables = 20;
    int populationSize = 10;
    int maxEvaluations = 1000000;

    Problem problem   ;         // The problem to solve
    Algorithm algorithm ;         // The algorithm to use

    //problem = new Sphere("Real", numberOfVariables) ;
    //problem = new Easom("Real") ;
    //problem = new Griewank("Real", populationSize) ;
    //problem = new Schwefel("Real", numberOfVariables) ;
    problem = new Rosenbrock("Real", numberOfVariables) ;
    //problem = new Rastrigin("Real", numberOfVariables) ;

    algorithm = new CMAES(problem) ;
    
    /* Algorithm parameters*/
    algorithm.setInputParameter("populationSize", populationSize);
    algorithm.setInputParameter("maxEvaluations", maxEvaluations);
 
    /* Execute the Algorithm */
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    long estimatedTime = System.currentTimeMillis() - initTime;
    System.out.println("Total execution time: " + estimatedTime);

    /* Log messages */
    System.out.println("Objectives values have been written to file FUN");
    population.printObjectivesToFile("FUN");
    System.out.println("Variables values have been written to file VAR");
    population.printVariablesToFile("VAR");

  } //main

} // CMAES_main
