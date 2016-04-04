package ares.fbk.eu.pMOEA.jmetal.operators.selection;

import  ares.fbk.eu.pMOEA.jmetal.operators.selection.Selection;
import  ares.fbk.eu.pMOEA.jmetal.core.Solution;
import ares.fbk.eu.pMOEA.jmetal.core.SolutionSet;


import com.sun.jmx.mbeanserver.JmxMBeanServerBuilder;

import ares.fbk.eu.pMOEA.jmetal.util.comparators.ObjectiveSpaceAndDominanceComparator;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;
import ares.fbk.eu.pMOEA.jmetal.util.PseudoRandom;

import java.awt.font.NumericShaper;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class PreferredTournamentSelection extends Selection {

	private Comparator objectiveSpaceComparator;
	private  int numberOfTournaments;
	private  int numberOFSolutionsPerRegion[];

	private int currentRegionNumber, numberOfRegions;

	

	/** Constructor */
	public PreferredTournamentSelection(HashMap<String, Object> parameters) {
		super(parameters);
		if(parameters == null){
			new JMException("preferred trounament selection parameters are null");
		}
		else{
			numberOfTournaments = (int) parameters.get("numberOfTournaments");
			numberOFSolutionsPerRegion = (int[]) parameters.get("numberOfSolutionsPerRegion");
			numberOfRegions = (int) parameters.get("numberOfRegions");
			objectiveSpaceComparator = new ObjectiveSpaceAndDominanceComparator();
			
		}
	}



	
	@Override
	/** Execute() method */
	public Object execute(Object object) {
		if(object == null){
			new JMException("problem in Preferred tournament execute methos");
		}
		
		currentRegionNumber = (int) getParameter("currentRegion");
			
		SolutionSet solutionList = (SolutionSet)object;
		
		 
		Solution result = null;

	if (PseudoRandom.randDouble() > 0.20) {
			solutionList = findSolutionsOfSameRegion(solutionList);

			if (solutionList.size() == 1) {
				result = solutionList.get(0);
			} else {
				result = selectNRandomDifferentSolutions(1, solutionList)
						.get(0);
				int count = 1; // at least 2 solutions are compared
				do {
					Solution candidate = selectRandomDifferentSolutions(solutionList,
							result);
					result = getBestSolution(result, candidate,	objectiveSpaceComparator);
				} while (++count < this.numberOfTournaments);
			}
		} else {
			if (solutionList.size() == 1) {
				result = solutionList.get(0);
			} else {
				
				/*result = selectNRandomDifferentSolutions(1, solutionList)
						.get(0);
				int count = 1; // at least 2 solutions are compared
				do {
					S candidate = selectRandomDifferentSolutions(solutionList,
							result);
					result = SolutionUtils.getBestSolution(result, candidate,
							dominanceComparator);
				} while (++count < this.numberOfTournaments);*/
				if(numberOfRegions>1){ // check: is it ok?
				int count=0;
				while(count<=0){
					int position = PseudoRandom.randInt(0, solutionList.size() - 1); 
					Integer regionNumber = solutionList.get(position).getAssignedRegionNumber();
					if(regionNumber == this.currentRegionNumber )
						continue;
					else{
						result = solutionList.get(position);
						count++;
					}
						
				}
				}else{
					result = solutionList.get(PseudoRandom.randInt(0, solutionList.size() - 1));
				}
				
			}
		}

		return result;
	}
	
	 public static Solution getBestSolution(Solution solution1, Solution solution2, Comparator comparator) {
		    Solution result ;
		    int flag = comparator.compare(solution1, solution2);
		    if (flag == -1) {
		      result = solution1;
		    } else if (flag == 1) {
		      result = solution2;
		    } else {
		      if (PseudoRandom.randDouble() < 0.5) {
		        result = solution1;
		      } else {
		        result = solution2;
		      }
		    }

		    return result ;
		  }
		  

	public SolutionSet findSolutionsOfSameRegion(
			SolutionSet solutionSet) {
		SolutionSet resultSolutionSet = new SolutionSet(numberOFSolutionsPerRegion[currentRegionNumber]);
		for (int i=0;i<solutionSet.size();i++) {
			Solution sol = solutionSet.get(i);
			if ( sol.getAssignedRegionNumber() == currentRegionNumber) {
				resultSolutionSet.add(sol);
			}
		}
		return resultSolutionSet;
	}

	public SolutionSet selectNRandomDifferentSolutions(
			int numberOfSolutionsToBeReturned, SolutionSet solutionList) {
		if (null == solutionList) {
			try {
				throw new JMException("The solution list is null");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (solutionList.size() == 0) {
			try {
				throw new JMException("The solution list is empty");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (solutionList.size() < numberOfSolutionsToBeReturned) {
			try {
				throw new JMException("The solution list size ("
						+ solutionList.size() + ") is less than "
						+ "the number of requested solutions ("
						+ numberOfSolutionsToBeReturned + ")");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		SolutionSet resultList = new SolutionSet(numberOfSolutionsToBeReturned);

		if (solutionList.size() == 1) {
			resultList.add(solutionList.get(0));
		} else {
			Collection<Integer> positions = new HashSet<>(
					numberOfSolutionsToBeReturned);
			while (positions.size() < numberOfSolutionsToBeReturned) {
				int nextPosition = PseudoRandom.randInt(0,
						solutionList.size() - 1);
				if (!positions.contains(nextPosition)) {
					positions.add(nextPosition);
					resultList.add(solutionList.get(nextPosition));
				}
			}
		}

		return resultList;
	}

	public  Solution selectRandomDifferentSolutions(
			SolutionSet solutionList, Solution solution) {
		if (null == solutionList) {
			try {
				throw new JMException("The solution list is null");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (solutionList.size() == 0) {
			try {
				throw new JMException("The solution list is empty");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		
		SolutionSet resultList = new SolutionSet(1);

		if (solutionList.size() == 1) {
			resultList.add(solutionList.get(0));
		} else {
			Collection<Integer> positions = new HashSet<>(1);
			int numberOfAttempt=0;
			while (positions.size() < 1 && numberOfAttempt < solutionList.size()  ) {
				
				int nextPosition = PseudoRandom.randInt(0, solutionList.size() - 1);
				if (!positions.contains(nextPosition)
						&& !solutionList.get(nextPosition).equals(solution)) {
					positions.add(nextPosition);
					resultList.add(solutionList.get(nextPosition));
					break;
					
				}
				numberOfAttempt++;
			}
			if(numberOfAttempt>= solutionList.size()){
				resultList.add(solution);
			}
		}

		return resultList.get(0);
	}

}
