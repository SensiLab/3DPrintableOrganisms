package generative.components.examples;

import generative.components.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

//import java.util.Arrays;

import fr.inria.optimization.cmaes.CMAEvolutionStrategy;
import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;

import java.util.Arrays;
import java.util.Date;
import java.util.ArrayList;

class TotalFitness implements IObjectiveFunction{
	Environment e;
	int numOfUses;
	TotalFitness() {
		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, 12345);
		this.numOfUses = 0;
	}
	
	public double valueOf (double[] x) {
		this.numOfUses++;
		
		// Initialise result variable
		double res = 0;
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);

		Simulation s = new Simulation(c, env, 500,20);
//		
		s.generate();
		
		// High is good
		float printability = s.printability();

		res+=(double)(printability * printability);
		
		// High is bad
		float convexity = 1f - s.convexity();
		
		res+=(double)(convexity * convexity);
		
		// High is bad
		float compactness = 1f - s.compactness();
		res+=(double)(compactness*compactness);
		
		// High is good
		float angleDispersion = s.angleDispersionCoefficient();
		res+=(double)(angleDispersion * angleDispersion);
		
		res/=4;
		res = res * s.isComplete();
		
		return 1 - res;
	}
	
	// This method checks if the chromosome is within valid range of values
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

/*
 * A fitness function to explore different printability scores
 */

class Printability implements IObjectiveFunction{
	Environment e;
	int numOfUses;
	float target;
	Printability(float _target) {
		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, 12345);
		this.numOfUses = 0;
		this.target = _target;
	}
	
	public double valueOf (double[] x) {
		this.numOfUses++;
		
		// Initialise result variable
		double res = 0;
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);

		Simulation s = new Simulation(c, env, 500,20);
//		
		s.generate();
		
		// Get printability
		float printability = s.printability();

		return Math.abs(printability - target);
	}
	
	// This method checks if the chromosome is within valid range of values
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

public class Example1 {
	public static void main(String[] args) {
		// Instantiate fitness function
		IObjectiveFunction fitfun = new Printability(0.3f);
		
		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
		cma.setDimension(5); // overwrite some loaded properties
		cma.setInitialX(0.5); // in each dimension, also setTypicalX can be used
		cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
//		cma.options.stopFitness = 1e-14;       // optional setting
		cma.options.stopMaxIter = 40;
		
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm");
		LocalDateTime now = LocalDateTime.now();
		String date = dtf.format(now);
		String path = "../output/";
		String fitnessDataFile = "FitnessData_"+date;
		String genesDataFile = "GenesData_"+date;
		
		// sets the prefix for the output files
//		cma.options.outputFileNamesPrefix = "C:\\Users\\camil\\01_Other Work\\Code Dev\\out\\"+date+"cma_out_";
		cma.options.outputFileNamesPrefix = "../output/"+date+"cma_out_";
		
		// Initialise output data (fitness)
		String fitnessData = "(randomSeed=" + cma.getSeed() + ", " + new Date().toString() + ")\n";
		// Add header
		fitnessData+="iteration, evaluations, sigma, axisratio, bestever_fitness, best_fitness, median_fitness, worst_fitness, mindii, "
		+ "idxmaxSD, maxSD, idxminSD, minSD \n";
		
		// Initialise output data (genetic information)
		String genesData = "(randomSeed=" + cma.getSeed() + ", " + new Date().toString() + ")\n";
		// Add header
		genesData+="iteration, evaluations, sigma, void, fitness_of_recent_best, x_of_recent_best(1...dimension) \n";
		
		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

		// initial output to files
		cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files
		
		// ARRAY TO STORE BEST SOLUTION PER ITERATION
		ArrayList<Float> bestFunctionPerIter = new ArrayList<Float>();
		ArrayList<double[]> bestSolutionPerInterval = new ArrayList<double[]>();

		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {

			// --- core iteration step ---
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				// a simple way to handle constraints that define a convex feasible domain  
				// (like box constraints, i.e. variable boundaries) via "blind re-sampling" 
		        // assumes that the feasible domain is convex, the optimum is  
				while (!fitfun.isFeasible(pop[i]))     //   not located on (or very close to) the domain boundary,  
					pop[i] = cma.resampleSingle(i);    //   initialX is feasible and initialStandardDeviations are  
		                                                       //   sufficiently small to prevent quasi-infinite looping here
				// compute fitness/objective value	
				fitness[i] = fitfun.valueOf(pop[i]); // fitfun.valueOf() is to be minimized
			}
			cma.updateDistribution(fitness);         // pass fitness array to update search distribution
			// --- end core iteration step ---

//			// output to files
//			cma.writeToDefaultFiles();
			
			// Add data to data strings
			fitnessData+=cma.getDataRowFitness();
			genesData+=cma.getDataRowXRecentBest();
			
			// output to console
			int outmod = 5;
			if (cma.getCountIter() % (15*outmod) == 1) {
				cma.printlnAnnotation(); // might write file as well
			}
			if (cma.getCountIter() % outmod == 1) {
				cma.println();
				bestSolutionPerInterval.add(cma.getBestRecentX());
			}
			// ADD BEST VALUE TO STORING ARRAY
			bestFunctionPerIter.add((float)cma.getBestRecentFunctionValue());
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 

		// final output
//		cma.writeToDefaultFiles(2);
		
		// write data strings to output file
		saveData(path, fitnessDataFile, fitnessData);
		saveData(path, genesDataFile, genesData);
		
		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
					+ " at evaluation " + cma.getBestEvaluationNumber());
					
		// we might return cma.getBestSolution() or cma.getBestX()
		cma.println("best solution" + Arrays.toString(cma.getBestX()));
		System.out.println(bestFunctionPerIter);
		for(double[] sln : bestSolutionPerInterval)
			System.out.println("Iteration "+bestSolutionPerInterval.indexOf(sln)+" solution: "+Arrays.toString(sln));
		
	} //main
	
	static void saveData(String filePath, String fileName, String data) {
		// Create file
		File output = new File(filePath+fileName+".csv");
		try {
			FileWriter outputWriter = new FileWriter(filePath+fileName+".csv");
			outputWriter.write(data);
			outputWriter.close();
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}

} //class
