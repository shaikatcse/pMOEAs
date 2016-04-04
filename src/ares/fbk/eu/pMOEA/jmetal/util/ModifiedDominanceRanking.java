//  Ranking.java
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

package ares.fbk.eu.pMOEA.jmetal.util;

import ares.fbk.eu.pMOEA.jmetal.core.Solution;
import ares.fbk.eu.pMOEA.jmetal.core.SolutionSet;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.ObjectiveSpaceAndDominanceComparator;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.OverallConstraintViolationComparator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;

import javax.xml.bind.annotation.XmlElementDecl.GLOBAL;

/**
 * This class implements some facilities for ranking solutions. Given a
 * <code>SolutionSet</code> object, their solutions are ranked according to
 * scheme proposed in NSGA-II; as a result, a set of subsets are obtained. The
 * subsets are numbered starting from 0 (in NSGA-II, the numbering starts from
 * 1); thus, subset 0 contains the non-dominated solutions, subset 1 contains
 * the non-dominated solutions after removing those belonging to subset 0, and
 * so on.
 */
public class ModifiedDominanceRanking {
	
	private SolutionSet solutionSet_;

	/**
	 * The <code>SolutionSet</code> to rank
	 */
	private SolutionSet localPOpulation_, globalPopulation_;

	/**
	 * An array containing all the fronts found during the search
	 */
	private SolutionSet[] ranking_;

	/**
	 * stores a <code>Comparator</code> for dominance checking
	 */
	private static final Comparator dominance_ = new DominanceComparator();

	/**
	 * stores a <code>Comparator</code> for Overal Constraint Violation
	 * Comparator checking
	 */
	private static final Comparator constraint_ = new OverallConstraintViolationComparator();

	private static final Comparator Objective_ = new ObjectiveSpaceAndDominanceComparator();

	Random rm = new Random();
	
	/**
	 * Constructor.
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code> to be ranked.
	 */
	public ModifiedDominanceRanking(SolutionSet localPopulation, SolutionSet globalPopulation, int evaluations, int maxEvaluations) {
		//solutionSet_ = solutionSet;
		this.localPOpulation_ = localPopulation;
		this.globalPopulation_ = globalPopulation;
		
		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[localPopulation.size()];

		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[localPopulation.size()];

		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[localPopulation.size() + 1];
		
		//dominatedByGp[i] contains the number of individual in global population dominate ith individual of local population
		int[] dominateMeByGp = new int[localPopulation.size()];

		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();

		/*
		 * //-> Fast non dominated sorting algorithm for (int p = 0; p <
		 * solutionSet_.size(); p++) { // Initialice the list of individuals
		 * that i dominate and the number // of individuals that dominate me
		 * iDominate[p] = new LinkedList<Integer>(); dominateMe[p] = 0; // For
		 * all q individuals , calculate if p dominates q or vice versa for (int
		 * q = 0; q < solutionSet_.size(); q++) { flagDominate
		 * =constraint_.compare(solutionSet.get(p),solutionSet.get(q)); if
		 * (flagDominate == 0) { flagDominate
		 * =dominance_.compare(solutionSet.get(p),solutionSet.get(q)); }
		 * 
		 * if (flagDominate == -1) { iDominate[p].add(new Integer(q)); } else if
		 * (flagDominate == 1) { dominateMe[p]++; } }
		 * 
		 * // If nobody dominates p, p belongs to the first front if
		 * (dominateMe[p] == 0) { front[0].add(new Integer(p));
		 * solutionSet.get(p).setRank(0); } }
		 */

		// -> Fast non dominated sorting algorithm
		// Contribution of Guillaume Jacquenot
		for (int p = 0; p < localPopulation.size(); p++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[p] = new LinkedList<>();
			dominateMe[p] = 0;
		}
		
		//adapted from pAGE strategy
		boolean trackStrategy = false; 
		if( rm.nextDouble() < Math.pow( (evaluations/(double)maxEvaluations) , 10)){
				trackStrategy = true;
		}
		
		int flagDominate = 0;
		
		for (int p = 0; p < (localPOpulation_.size() - 1); p++) {
			// For all q individuals , calculate if p dominates q or vice versa
			for (int q = p + 1; q < localPOpulation_.size(); q++) {
				flagDominate=0;
				if(trackStrategy){
				 flagDominate = Objective_.compare(
							 localPopulation.get(p), localPopulation.get(q));
				}else{
					if (flagDominate == 0) {
						flagDominate = constraint_.compare(localPopulation.get(p),
								localPopulation.get(q));
					}
					if (flagDominate == 0) {
						flagDominate = dominance_.compare(localPopulation.get(p),
								localPopulation.get(q));
					}
				}
				
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
			// If nobody dominates p, p belongs to the first front
		}
		
		for (int p = 0; p < localPOpulation_.size(); p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
				//check
				//solutionSet.get(p).setRank(0);
			}
		}
		
		// global dominance
		for (int p = 0; p < localPopulation.size() ; p++) {
			for (int z=0;z<globalPopulation_.size();z++) {
				Solution globalSolution = globalPopulation_.get(z);
				int flagP = dominance_.compare(globalSolution,
						localPopulation.get(p));
				if (flagP == -1)
					dominateMeByGp[p]++;
					//localPopulation.get(p).setDominateMeByGP(localPopulation.get(p).getDominateMeByGP()+1);
			}
		}


		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
						//check
						//solutionSet_.get(index).setRank(i);
					}
				}
			}
		}
		// <-

		//Shahriar modification
		class SortingClass {
			int index;
			int dominateMeByGp;
			
			SortingClass(int index, int dominateMeByGp){
				this.index = index;
				this.dominateMeByGp = dominateMeByGp;
			}
			
			int getIndex(){
				return index;
			}
			
			int getDominateMeByGp(){
				return this.dominateMeByGp;
			}
		}
		

		int rankingnumber=0;
		List<ArrayList<Solution>> rankedSubpopulations = new ArrayList<>();
		
		for (int j = 0; j < i; j++) {
			List<SortingClass> sortingClass = new ArrayList<SortingClass>(front[j].size());
			
			for(int z=0;z<front[j].size();z++){
				sortingClass.add(new SortingClass(front[j].get(z), dominateMeByGp[front[j].get(z)]));
				//sortingClass.add(new SortingClass(front[j].get(z), localPopulation.get(front[j].get(z)).getDominateMeByGP() ));
			}
			//sort sortingClass according to dominateMeByGp
			Collections.sort(sortingClass, new Comparator<SortingClass>() {
				@Override public int compare(final SortingClass o1, final SortingClass o2) {
					if(o1.getDominateMeByGp() < o2.getDominateMeByGp())
						return -1;
					else if(o1.getDominateMeByGp() > o2.getDominateMeByGp())
						return 1;
					else
						return 0;
					
				}
			});
		
			//find number of different dominateMeByGp value exits in sortingClass
			Set<Integer> setUniqueNumbers = new LinkedHashSet<Integer>();
			for(SortingClass sc: sortingClass){
				setUniqueNumbers.add(sc.getDominateMeByGp());
			}
			int numberOfDifferentDominateMeByGpValues = setUniqueNumbers.size();
			
			//add all different DominateMeByGp in a list
			List<Integer> list = new ArrayList<Integer>(sortingClass.size());
			for(SortingClass sc: sortingClass){
				list.add(sc.getDominateMeByGp());
			}
			
			int sortingClassCounter=0;
			 
			for(int z=0; z<numberOfDifferentDominateMeByGpValues; z++){
				
				//count number of frequency of different DominateMeByGp Values
				int numberOfindvInSameRank = Collections.frequency(list,(Integer)setUniqueNumbers.toArray()[z] );
							
				rankedSubpopulations.add(rankingnumber, new ArrayList<Solution>(numberOfindvInSameRank));
				//ranking_[rankingnumber] = new SolutionSet (numberOfindvInSameRank);
				
				for(int i1=0;i1<numberOfindvInSameRank;i1++){
					
					localPOpulation_.get(sortingClass.get(sortingClassCounter).getIndex()).setRank(rankingnumber);
					rankedSubpopulations.get(rankingnumber).add(localPopulation.get(sortingClass.get(sortingClassCounter).getIndex()));
					
					sortingClassCounter++;					
											
				}
				
				rankingnumber++;
				
			}
			
		}
		
	
			
		/*ranking_ = new SolutionSet[i];
		// 0,1,2,....,i-1 are front, then i fronts
		for (int j = 0; j < i; j++) {
			ranking_[j] = new SolutionSet(front[j].size());
			it1 = front[j].iterator();
			while (it1.hasNext()) {
				ranking_[j].add(solutionSet.get(it1.next()));
			}
		}*/
		
		ranking_ = new SolutionSet[rankedSubpopulations.size()];
		for(int zz=0;zz<rankedSubpopulations.size();zz++){
			ranking_[zz] = new SolutionSet(rankedSubpopulations.get(zz).size());
			for(int jj=0;jj<rankedSubpopulations.get(zz).size();jj++){
				ranking_[zz].add(rankedSubpopulations.get(zz).get(jj));
			}
		}

	} // Ranking
	
	
	public ModifiedDominanceRanking(SolutionSet solutionSet, int evaluations, int maxEvaluations) {
		solutionSet_ = solutionSet;

		// dominateMe[i] contains the number of solutions dominating i
		int[] dominateMe = new int[solutionSet_.size()];

		// iDominate[k] contains the list of solutions dominated by k
		List<Integer>[] iDominate = new List[solutionSet_.size()];

		// front[i] contains the list of individuals belonging to the front i
		List<Integer>[] front = new List[solutionSet_.size() + 1];

		// flagDominate is an auxiliar encodings.variable
		int flagDominate;

		// Initialize the fronts
		for (int i = 0; i < front.length; i++)
			front[i] = new LinkedList<Integer>();

		/*
		 * //-> Fast non dominated sorting algorithm for (int p = 0; p <
		 * solutionSet_.size(); p++) { // Initialice the list of individuals
		 * that i dominate and the number // of individuals that dominate me
		 * iDominate[p] = new LinkedList<Integer>(); dominateMe[p] = 0; // For
		 * all q individuals , calculate if p dominates q or vice versa for (int
		 * q = 0; q < solutionSet_.size(); q++) { flagDominate
		 * =constraint_.compare(solutionSet.get(p),solutionSet.get(q)); if
		 * (flagDominate == 0) { flagDominate
		 * =dominance_.compare(solutionSet.get(p),solutionSet.get(q)); }
		 * 
		 * if (flagDominate == -1) { iDominate[p].add(new Integer(q)); } else if
		 * (flagDominate == 1) { dominateMe[p]++; } }
		 * 
		 * // If nobody dominates p, p belongs to the first front if
		 * (dominateMe[p] == 0) { front[0].add(new Integer(p));
		 * solutionSet.get(p).setRank(0); } }
		 */

		// -> Fast non dominated sorting algorithm
		// Contribution of Guillaume Jacquenot
		for (int p = 0; p < solutionSet_.size(); p++) {
			// Initialize the list of individuals that i dominate and the number
			// of individuals that dominate me
			iDominate[p] = new LinkedList<Integer>();
			dominateMe[p] = 0;
		}
		//adapted from pAGE strategy
		boolean trackStrategy = false; 
		if( rm.nextDouble() < Math.pow( (evaluations/(double)maxEvaluations) , 5)){
			trackStrategy = true;
		}
		
		 for (int p = 0; p < (solutionSet_.size() - 1); p++) {
			// For all q individuals , calculate if p dominates q or vice versa
			for (int q = p + 1; q < solutionSet_.size(); q++) {
				flagDominate=0;
				if(trackStrategy){
				 flagDominate = Objective_.compare(
							 solutionSet.get(p), solutionSet.get(q));
					 
					 //for debug:
					 if (evaluations==4830){
						 System.out.println("check for erro");
						 solutionSet_.printObjectivesToFile("FunDebug");
					 }
				}else{
				if (flagDominate == 0) {
					flagDominate = constraint_.compare(solutionSet.get(p),
							solutionSet.get(q));
				}
				if (flagDominate == 0) {
					flagDominate = dominance_.compare(solutionSet.get(p),
							solutionSet.get(q));
				}}
				if (flagDominate == -1) {
					iDominate[p].add(q);
					dominateMe[q]++;
				} else if (flagDominate == 1) {
					iDominate[q].add(p);
					dominateMe[p]++;
				}
			}
			// If nobody dominates p, p belongs to the first front
		}
		for (int p = 0; p < solutionSet_.size(); p++) {
			if (dominateMe[p] == 0) {
				front[0].add(p);
				solutionSet.get(p).setRank(0);
			}
		}

		// Obtain the rest of fronts
		int i = 0;
		Iterator<Integer> it1, it2; // Iterators
		while (front[i].size() != 0) {
			i++;
			it1 = front[i - 1].iterator();
			while (it1.hasNext()) {
				it2 = iDominate[it1.next()].iterator();
				while (it2.hasNext()) {
					int index = it2.next();
					dominateMe[index]--;
					if (dominateMe[index] == 0) {
						front[i].add(index);
						solutionSet_.get(index).setRank(i);
					}
				}
			}
		}
		// <-

		ranking_ = new SolutionSet[i];
		// 0,1,2,....,i-1 are front, then i fronts
		for (int j = 0; j < i; j++) {
			ranking_[j] = new SolutionSet(front[j].size());
			it1 = front[j].iterator();
			while (it1.hasNext()) {
				ranking_[j].add(solutionSet.get(it1.next()));
			}
		}

	} // Ranking

	/**
	 * Returns a <code>SolutionSet</code> containing the solutions of a given
	 * rank.
	 * 
	 * @param rank
	 *            The rank
	 * @return Object representing the <code>SolutionSet</code>.
	 */
	public SolutionSet getSubfront(int rank) {
		return ranking_[rank];
	} // getSubFront

	/**
	 * Returns the total number of subFronts founds.
	 */
	public int getNumberOfSubfronts() {
		return ranking_.length;
	} // getNumberOfSubfronts
} // Ranking
