//  Distance.java
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
import ares.fbk.eu.pMOEA.jmetal.util.comparators.ObjectiveComparator;
import ares.fbk.eu.pMOEA.jmetal.util.wrapper.XReal;

/**
 * This class implements some utilities for calculating distances
 */
public class Distance {

	/**
	 * Constructor.
	 */
	public Distance() {
		// do nothing.
	} // Distance

	/**
	 * Returns a matrix with distances between solutions in a
	 * <code>SolutionSet</code>.
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code>.
	 * @return a matrix with distances.
	 */
	public double[][] distanceMatrix(SolutionSet solutionSet) {
		Solution solutionI, solutionJ;

		// The matrix of distances
		double[][] distance = new double[solutionSet.size()][solutionSet.size()];
		// -> Calculate the distances
		for (int i = 0; i < solutionSet.size(); i++) {
			distance[i][i] = 0.0;
			solutionI = solutionSet.get(i);
			for (int j = i + 1; j < solutionSet.size(); j++) {
				solutionJ = solutionSet.get(j);
				distance[i][j] = this.distanceBetweenObjectives(solutionI,
						solutionJ);
				distance[j][i] = distance[i][j];
			} // for
		} // for

		// ->Return the matrix of distances
		return distance;
	} // distanceMatrix

	/**
	 * Returns the minimum distance from a <code>Solution</code> to a
	 * <code>SolutionSet according to the objective values</code>.
	 * 
	 * @param solution
	 *            The <code>Solution</code>.
	 * @param solutionSet
	 *            The <code>SolutionSet</code>.
	 * @return The minimum distance between solution and the set.
	 * @throws JMException
	 */
	public double distanceToSolutionSetInObjectiveSpace(Solution solution,
			SolutionSet solutionSet) throws JMException {
		// At start point the distance is the max
		double distance = Double.MAX_VALUE;

		// found the min distance respect to population
		for (int i = 0; i < solutionSet.size(); i++) {
			double aux = this.distanceBetweenObjectives(solution,
					solutionSet.get(i));
			if (aux < distance)
				distance = aux;
		} // for

		// ->Return the best distance
		return distance;
	} // distanceToSolutionSetinObjectiveSpace

	/**
	 * Returns the minimum distance from a <code>Solution</code> to a
	 * <code>SolutionSet according to the encodings.variable values</code>.
	 * 
	 * @param solution
	 *            The <code>Solution</code>.
	 * @param solutionSet
	 *            The <code>SolutionSet</code>.
	 * @return The minimum distance between solution and the set.
	 * @throws JMException
	 */
	public double distanceToSolutionSetInSolutionSpace(Solution solution,
			SolutionSet solutionSet) throws JMException {
		// At start point the distance is the max
		double distance = Double.MAX_VALUE;

		// found the min distance respect to population
		for (int i = 0; i < solutionSet.size(); i++) {
			double aux = this.distanceBetweenSolutions(solution,
					solutionSet.get(i));
			if (aux < distance)
				distance = aux;
		} // for

		// ->Return the best distance
		return distance;
	} // distanceToSolutionSetInSolutionSpace

	/**
	 * Returns the distance between two solutions in the search space.
	 * 
	 * @param solutionI
	 *            The first <code>Solution</code>.
	 * @param solutionJ
	 *            The second <code>Solution</code>.
	 * @return the distance between solutions.
	 * @throws JMException
	 */
	public double distanceBetweenSolutions(Solution solutionI,
			Solution solutionJ) throws JMException {
		/*
		 * double distance = 0.0; if ((solutionI.getDecisionVariables() != null)
		 * && (solutionJ.getDecisionVariables() != null)) { Variable[]
		 * decisionVariableI = solutionI.getDecisionVariables(); Variable[]
		 * decisionVariableJ = solutionJ.getDecisionVariables();
		 * 
		 * double diff; //Auxiliar var //-> Calculate the Euclidean distance for
		 * (int i = 0; i < decisionVariableI.length; i++){ diff =
		 * decisionVariableI[i].getValue() - decisionVariableJ[i].getValue();
		 * distance += Math.pow(diff,2.0); } // for } //-> Return the euclidean
		 * distance return Math.sqrt(distance);
		 */
		double distance = 0.0;
		XReal solI = new XReal(solutionI);
		XReal solJ = new XReal(solutionJ);

		double diff; // Auxiliar var
		// -> Calculate the Euclidean distance
		for (int i = 0; i < solI.getNumberOfDecisionVariables(); i++) {
			diff = solI.getValue(i) - solJ.getValue(i);
			distance += Math.pow(diff, 2.0);
		} // for
		// -> Return the euclidean distance
		return Math.sqrt(distance);
	} // distanceBetweenSolutions

	/**
	 * Returns the distance between two solutions in objective space.
	 * 
	 * @param solutionI
	 *            The first <code>Solution</code>.
	 * @param solutionJ
	 *            The second <code>Solution</code>.
	 * @return the distance between solutions in objective space.
	 */
	public double distanceBetweenObjectives(Solution solutionI,
			Solution solutionJ) {
		double diff; // Auxiliar var
		double distance = 0.0;
		// -> Calculate the euclidean distance
		for (int nObj = 0; nObj < solutionI.getNumberOfObjectives(); nObj++) {
			diff = solutionI.getObjective(nObj) - solutionJ.getObjective(nObj);
			distance += Math.pow(diff, 2.0);
		} // for

		// Return the euclidean distance
		return Math.sqrt(distance);
	} // distanceBetweenObjectives.

	/**
	 * Return the index of the nearest solution in the solution set to a given
	 * solution
	 * 
	 * @param solution
	 * @param solutionSet
	 * @return The index of the nearest solution; -1 if the solutionSet is empty
	 */
	public int indexToNearestSolutionInSolutionSpace(Solution solution,
			SolutionSet solutionSet) {
		int index = -1;
		double minimumDistance = Double.MAX_VALUE;
		try {
			for (int i = 0; i < solutionSet.size(); i++) {
				double distance = 0;
				distance = distanceBetweenSolutions(solution,
						solutionSet.get(i));
				if (distance < minimumDistance) {
					minimumDistance = distance;
					index = i;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return index;
	}

	/**
	 * Assigns crowding distances to all solutions in a <code>SolutionSet</code>
	 * .
	 * 
	 * @param solutionSet
	 *            The <code>SolutionSet</code>.
	 * @param nObjs
	 *            Number of objectives.
	 */
	public void crowdingDistanceAssignment(SolutionSet solutionSet, int nObjs) {
		int size = solutionSet.size();

		if (size == 0)
			return;

		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		if (size == 2) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		// Use a new SolutionSet to evite alter original solutionSet
		SolutionSet front = new SolutionSet(size);
		for (int i = 0; i < size; i++) {
			front.add(solutionSet.get(i));
		}

		for (int i = 0; i < size; i++)
			front.get(i).setCrowdingDistance(0.0);

		double objetiveMaxn;
		double objetiveMinn;
		double distance;

		for (int i = 0; i < nObjs; i++) {
			// Sort the population by Obj n
			front.sort(new ObjectiveComparator(i));
			objetiveMinn = front.get(0).getObjective(i);
			objetiveMaxn = front.get(front.size() - 1).getObjective(i);

			// Set de crowding distance
			front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			front.get(size - 1).setCrowdingDistance(Double.POSITIVE_INFINITY);

			for (int j = 1; j < size - 1; j++) {
				distance = front.get(j + 1).getObjective(i)
						- front.get(j - 1).getObjective(i);
				distance = distance / (objetiveMaxn - objetiveMinn);
				distance += front.get(j).getCrowdingDistance();
				front.get(j).setCrowdingDistance(distance);
			} // for
		} // for
	} // crowdingDistanceAssing

	// new: for AGE

	public void crowdingDistanceAssignment(SolutionSet solutionSet, int nObjs,
			int typeSwitch, SolutionSet archive) {

		if (archive == null
				&& (typeSwitch == 3 || typeSwitch == 4 || typeSwitch == 6)) {
			System.out
					.println("\n\nchosen crowdingDistanceAssignment requires 'archive', but called with 'archive==null' --> please check");
			System.exit(1);
		}

		switch (typeSwitch) {
		case -1:
			// original
			crowdingDistanceAssignmentAGE(solutionSet, nObjs);
			break;

		case 0:
			// original
			crowdingDistanceAssignment(solutionSet, nObjs);
			break;

		case 1:
			// sum of distances to max values per objective
			crowdingDistanceAssignmentToExtremeValues(solutionSet, nObjs,
					false, null);
			break;
		case 2:
			// sum of distances to max values per objective (scaled)
			crowdingDistanceAssignmentToExtremeValues(solutionSet, nObjs, true,
					null);
			break;

		case 3:
			// sum of distances to max values per objective (w.r.t. archive)
			crowdingDistanceAssignmentToExtremeValues(solutionSet, nObjs,
					false, archive);
			break;
		case 4:
			// sum of distances to max values per objective (w.r.t. archive,
			// scaled)
			crowdingDistanceAssignmentToExtremeValues(solutionSet, nObjs, true,
					archive);
			break;

		case 5:
			// sum of 'size-ownRank' per objective
			crowdingDistanceAssignmentRanking(solutionSet, nObjs, null);
			break;
		case 6:
			// sum of 'size-ownRank' per objective (w.r.t. archive)
			crowdingDistanceAssignmentRanking(archive, nObjs, null);
			break;

		}
	}

	public void crowdingDistanceAssignmentAGE(SolutionSet solutionSet, int nObjs) {
		int size = solutionSet.size();
		boolean debugPrint = !true;

		if (size == 0)
			return;

		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		if (size == 2) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		// Use a new SolutionSet to avoid altering the original solutionSet
		SolutionSet front = new SolutionSet(size);
		for (int i = 0; i < size; i++) {
			front.add(solutionSet.get(i));
		}

		for (int i = 0; i < size; i++)
			front.get(i).setCrowdingDistance(0.0);

		double objetiveMaxn;
		double objetiveMinn;
		double distance;

		for (int i = 0; i < nObjs; i++) {
			// Sort the population by Obj n
			front.sort(new ObjectiveComparator(i));
			objetiveMinn = front.get(0).getObjective(i);
			objetiveMaxn = front.get(front.size() - 1).getObjective(i);
			// System.out.println(objetiveMinn + " " + objetiveMaxn); // min is
			// 0th element, max is nth element

			// Set the crowding distance
			front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY); // original
			front.get(size - 1).setCrowdingDistance(Double.POSITIVE_INFINITY); // original

			int middle = size / 2;
			// double middle = size/2d;
			for (int j = 1; j < size - 1; j++) {

				int pos = j;
				// if (pos>=middle) { /*do not do this! we go through each
				// objectiv: why include a bad one again?*/ // why this in
				// combination with falling probability?
				// // System.out.print("!");
				// pos = size-j-1;
				// }
				if (debugPrint)
					System.out.print(" pos=" + pos + " ");

				double rand = Math.random();

				if (debugPrint)
					System.out.print(front.get(j).getObjective(i));
				// if (Math.random()< (1d/ (pos) ) )
				// front.get(j).setCrowdingDistance(Double.POSITIVE_INFINITY);//falling
				// probability toward the center of the front
				/* use */if (rand < (1d / (pos + 1))) {
					// if (rand< (1d) ) {
					front.get(j).setCrowdingDistance(Double.POSITIVE_INFINITY); // falling
																				// probability
																				// toward
																				// the
																				// center
																				// of
																				// the
																				// front
					if (debugPrint)
						System.out.println(" {" + rand + "<" + 1d / (pos + 1)
								+ "} y");
				} else {
					if (front.get(j).getCrowdingDistance() != Double.POSITIVE_INFINITY)
						front.get(j).setCrowdingDistance(0); // reset to zero if
																// not already
																// INFINITY
					if (debugPrint)
						System.out.println();
				}
				// if (Math.random()< 1d/nObjs )
				// front.get(j).setCrowdingDistance(Double.POSITIVE_INFINITY);//1/d
				// probability DTLZ3_3D=eps3
				// if (Math.random()< .5 )
				// front.get(j).setCrowdingDistance(Double.POSITIVE_INFINITY);//50%
				// probability DTLZ3_3D=eps2
				// if (Math.random()< 1 )
				// front.get(j).setCrowdingDistance(size-pos +
				// front.get(j).getCrowdingDistance());

				// if (Math.random()< 1 )
				// front.get(j).setCrowdingDistance(Double.POSITIVE_INFINITY);//100%
				// probability DTLZ3_3D=eps4

				// distance = front.get(j+1).getObjective(i) -
				// front.get(j-1).getObjective(i);
				// distance = distance / (objetiveMaxn - objetiveMinn);
				// distance += front.get(j).getCrowdingDistance();
				// front.get(j).setCrowdingDistance(distance);
				// System.out.print(front.get(j).getCrowdingDistance() + " ");
			} // for
			if (debugPrint)
				System.out.println();
		} // for
		if (debugPrint)
			System.out.println();
	} // crowdingDistanceAssing

	/*
	 * x_i = max { (max{y\inY} f_i(y)) - f_i(x) ,0} (code based on
	 * crowdingDistanceAssignment/2)
	 */
	public void crowdingDistanceAssignmentToExtremeValues(
			SolutionSet solutionSet, int nObjs, boolean scale,
			SolutionSet archive) {
		int size = solutionSet.size();

		if (size == 0)
			return;

		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		if (size == 2) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		// Use a new SolutionSet to evite alter original solutionSet
		SolutionSet front = new SolutionSet(size);
		for (int i = 0; i < size; i++) {
			front.add(solutionSet.get(i));
		}

		for (int i = 0; i < size; i++)
			front.get(i).setCrowdingDistance(0.0);

		double objectiveMaxn, archiveMaxn;
		double objectiveMinn, archiveMinn;
		double distance;

		for (int i = 0; i < nObjs; i++) {
			// Sort the population by Obj n
			front.sort(new ObjectiveComparator(i));
			objectiveMinn = front.get(0).getObjective(i);
			objectiveMaxn = front.get(front.size() - 1).getObjective(i);

			front.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY); // original
			front.get(size - 1).setCrowdingDistance(Double.POSITIVE_INFINITY); // original

			// modifications...

			if (archive != null) {
				archive.sort(new ObjectiveComparator(i));
				archiveMinn = archive.get(0).getObjective(i);
				archiveMaxn = archive.get(archive.size() - 1).getObjective(i);
				objectiveMinn = Math.min(objectiveMinn, archiveMinn);
				objectiveMaxn = Math.max(objectiveMaxn, archiveMaxn);
			}

			for (int j = 1; j < size - 1; j++) {
				// for (int j = 0; j < size; j++) {

				// distance = Math.max(
				// front.get(front.size()-1).getObjective(i) /
				// front.get(j).getObjective(i) , 1); // trial
				distance = Math.max(front.get(front.size() - 1).getObjective(i)
						- front.get(j).getObjective(i), 0);

				if (scale) {
					distance = distance / (objectiveMaxn - objectiveMinn);
				}

				// if (distance==0) distance =1; // trial

				distance += front.get(j).getCrowdingDistance();
				front.get(j).setCrowdingDistance(distance);
			} // for
		} // for
	} // crowdingDistanceAssing

	public void crowdingDistanceAssignmentRanking(SolutionSet solutionSet,
			int nObjs, SolutionSet archive) {
		int size = solutionSet.size();

		if (size == 0)
			return;

		if (size == 1) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		if (size == 2) {
			solutionSet.get(0).setCrowdingDistance(Double.POSITIVE_INFINITY);
			solutionSet.get(1).setCrowdingDistance(Double.POSITIVE_INFINITY);
			return;
		} // if

		// Use a new SolutionSet to evite alter original solutionSet
		SolutionSet front = new SolutionSet(size);
		for (int i = 0; i < size; i++) {
			front.add(solutionSet.get(i));
		}

		for (int i = 0; i < size; i++)
			front.get(i).setCrowdingDistance(0.0);

		double objectiveMaxn, archiveMaxn;
		double objectiveMinn, archiveMinn;
		double distance;

		for (int i = 0; i < nObjs; i++) {
			// Sort the population by Obj n
			front.sort(new ObjectiveComparator(i));

			// modifications...

			for (int j = 0; j < size; j++) {
				distance = size - j - 1;// j is better than size-j-1

				distance += front.get(j).getCrowdingDistance();
				front.get(j).setCrowdingDistance(distance);
			} // for
		} // for
	} // crowdingDistanceAssing

} // Distance

