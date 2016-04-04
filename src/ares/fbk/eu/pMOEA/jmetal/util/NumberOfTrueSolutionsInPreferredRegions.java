package ares.fbk.eu.pMOEA.jmetal.util;

import java.util.ArrayList;
import java.util.List;

import ares.fbk.eu.pMOEA.jmetal.qualityIndicator.Hypervolume;

public class NumberOfTrueSolutionsInPreferredRegions {
	
	static List<Double[]> regions= new ArrayList<Double[]>(3);
    
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		NumberOfTrueSolutionsInPreferredRegions numberOfTrueSolutionsInPreferredRegions = new NumberOfTrueSolutionsInPreferredRegions();
				
		regions.add(new Double[] {0.8, 0.95});
	    regions.add(new Double[] {0.40, 0.50});
	    regions.add(new Double[] {0.15, 0.20});
		
		//regions.add(new Double[] {0.8, 0.85});
		 //regions.add(new Double[] {0.35, 0.50});
		 //regions.add(new Double[] {0.15, 0.20});
	    
	    double[][] trueFront = new Hypervolume().utils_.readFront(args[0]);
	    
	    int numberOfTrueSolutions = 0;
	    
	 // identify regional true Paretofront
		for (int i = 0; i < trueFront.length; i++) {
			for (int regionNumber = 0; regionNumber < regions
					.size(); regionNumber++) {
				if (trueFront[i][0] >= regions
						.get(regionNumber)[0] /*left side*/
						&& trueFront[i][0] <= regions
								.get(regionNumber)[1]) /*right side*/ {
					numberOfTrueSolutions++;
				}
			}
		}
		
		System.out.println("Percentage: " + (double)numberOfTrueSolutions/trueFront.length*100);
				
	    
	}

}
