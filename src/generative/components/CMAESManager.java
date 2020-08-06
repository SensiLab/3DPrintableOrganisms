package generative.components;

import java.util.ArrayList;
import java.util.Arrays;

import fr.inria.optimization.cmaes.CMAEvolutionStrategy;
import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;
//import generative.components.examples.Printability;

class Evaluation implements IObjectiveFunction{
	Environment e;
	int numOfUses;
	Evaluation(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, 3);
		this.e = _e;
		this.numOfUses = 0;
	}
	
	public double valueOf (double[] x) {
		this.numOfUses++;
//		System.out.println("Fitness called "+this.numOfUses+" times");
		double res = 0;
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		System.out.println("Chromosome: "+c.chrToString());
//		System.out.println("Fitness: chromosome initialised");
		Simulation s = new Simulation(c, env, 500,100);
//		System.out.println("Fitness: simulation initialised");
		long genStart = System.currentTimeMillis();
		s.generate();
		long genDone = System.currentTimeMillis() - genStart;
//		System.out.println("Fitness: generation done ("+genDone+" ms)");
		
		long printStart = System.currentTimeMillis();
		float printability = s.printability();
		long printDone = System.currentTimeMillis() - printStart;
//		System.out.println("Fitness: printability = "+printability+" ("+genDone+" ms)");
//		System.out.println("Fitness: printability calculated: "+printability);
		res+=(double)printability;
		
		long convStart = System.currentTimeMillis();
		float convexity = 1f - s.convexity();
		long convDone = System.currentTimeMillis() - convStart;
//		System.out.println("Fitness: convexity = "+convexity+" ("+ convDone + " ms)");
		
//		System.out.println("Fitness: convexity calculated: "+convexity);
		res+=(double)convexity;
		
		long compStart = System.currentTimeMillis();
		float compactness = 1f - s.compactness();
		long compEnd = System.currentTimeMillis() - compStart;
		res+=(double)compactness;
		
		res/=3;
		res = res * s.isComplete();
		
		return 1 - res;
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

public class CMAESManager {
		
	public static void evolve(Environment e, int iterations) {
		IObjectiveFunction fitfun = new Evaluation(e);
	
		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
		cma.setDimension(5); // overwrite some loaded properties
		cma.setInitialX(0.5); // in each dimension, also setTypicalX can be used
		cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
		//	cma.options.stopFitness = 1e-14;       // optional setting
		cma.options.stopMaxIter = iterations;

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

			// output to files and console 
			cma.writeToDefaultFiles();
			int outmod = 10;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
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
		cma.writeToDefaultFiles(2);
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
	
	} //evolve
	
}//class
