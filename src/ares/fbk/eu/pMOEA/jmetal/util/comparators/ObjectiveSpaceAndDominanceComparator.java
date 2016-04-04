package ares.fbk.eu.pMOEA.jmetal.util.comparators;

import java.util.Comparator;

import ares.fbk.eu.pMOEA.jmetal.core.Solution;
import ares.fbk.eu.pMOEA.jmetal.util.JMException;
import ares.fbk.eu.pMOEA.jmetal.util.comparators.DominanceComparator;

public class ObjectiveSpaceAndDominanceComparator implements Comparator {
	private DominanceComparator dominanceComparator;

	public ObjectiveSpaceAndDominanceComparator() {
		this.dominanceComparator = new DominanceComparator();
	}

	/**
	 * Compares two solutions.
	 * 
	 * @param solution1
	 *            Object representing the first <code>Solution</code>.
	 * @param solution2
	 *            Object representing the second <code>Solution</code>.
	 * @return -1, or 0, or 1 if solution1 dominates solution2, both are
	 *         non-dominated, or solution1 is dominated by solution2,
	 *         respectively.
	 */

	public int compare(Object solution1, Object solution2) {
		if (solution1 == null) {
			try {
				throw new JMException("Solution1 is null");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else if (solution2 == null) {
			try {
				throw new JMException("Solution2 is null");
			} catch (JMException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		int result;
		@SuppressWarnings("unused")
		double sol1dis = ((Solution) solution1).getDistanceFromCloestRegion();
		@SuppressWarnings("unused")
		double sol2dis = ((Solution) solution2).getDistanceFromCloestRegion();
		if (sol1dis < sol2dis)
			result = -1;
		else if (sol1dis > sol2dis)
			result = 1;
		else
			result = 0;

		if (result == 0) {
			result = dominanceComparator.compare(solution1, solution2);
		}
		return result;
	}

}