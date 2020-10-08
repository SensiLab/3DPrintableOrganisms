package generative.components;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Random;

import fr.inria.optimization.cmaes.CMAEvolutionStrategy;
import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;
//import generative.components.examples.Printability;

/*
 * Extends the IObjectiveFunction interface to provide readable output
 * for partial fitness metrics
 */
interface FitnessFunction extends IObjectiveFunction{
	public long envRandomSeed(Environment _e);
}

/*
 * Measures fitness based on Printability and Complexity, weighted equally
 */
class FullFitness implements FitnessFunction{
	
	Environment e;
	public float printability;
	FullFitness(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);
		this.e = _e;
	}
	
	public long envRandomSeed(Environment _e) {
		return this.e.getRandomSeed();
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		
		// Initialise simulation
		Simulation s = new Simulation(c, env, 500,100);
		// Generate object
		s.generate();
		
		// Calculate printability
		//High is good so 1 - print is added to result
		float printability = s.printability();
		
		// add printability to overall results
		res+=(double)(printability);
		
		// calculate convexity
		//Low is good
		float convexity = (1 - s.convexity());
		
		//add convexity to overall results
		res+=(double)(convexity);
		
		// Add angle dispersion to overall result
		// High is good
		float angleDispersion = s.angleDispersionCoefficient();
		res+=(double)(angleDispersion);
		
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

class Convexity implements FitnessFunction{
	
	Environment e;
	public float printability;
	Convexity(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);
		this.e = _e;
	}
	
	public long envRandomSeed(Environment _e) {
		return this.e.getRandomSeed();
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		
		// Initialise simulation
		Simulation s = new Simulation(c, env, 500,100);
		// Generate object
		s.generate();
		
		// calculate convexity
		//Low is good
		return s.convexity();
		
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

class AngleDispersion implements FitnessFunction{
	
	Environment e;
	public float printability;
	AngleDispersion(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);
		this.e = _e;
	}
	
	public long envRandomSeed(Environment _e) {
		return this.e.getRandomSeed();
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		
		// Initialise simulation
		Simulation s = new Simulation(c, env, 500,100);
		// Generate object
		s.generate();
		
		// Add angle dispersion to overall result
		// High is good
		return (1 - s.angleDispersionCoefficient());
		
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

/*
 * Measures Printability
 */
class Printability implements FitnessFunction{
	
	Environment e;
	public float printability;
	Printability(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);
		this.e = _e;
	}
	
	public long envRandomSeed(Environment _e) {
		return this.e.getRandomSeed();
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		
		// Initialise simulation
		Simulation s = new Simulation(c, env, 500,100);
		// Generate object
		s.generate();
		
		// Calculate printability
		//High is good so 1 - print is added to result
		float printability = s.printability();
		
		// add printability to overall results
		res+=(double)(printability * printability);
		
		return 1 - res;
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}

/*
 * Complexity measure, based on convexity and angle dispersion
 */
class Complexity implements FitnessFunction{
	
	Environment e;
	public float printability;
	Complexity(Environment _e) {
//		this.e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);
		this.e = _e;
	}
	
	public long envRandomSeed(Environment _e) {
		return this.e.getRandomSeed();
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		// create chromosome from x
		Environment env = new Environment(this.e);
		Chromosome c = new Chromosome(x);
//		
		// Initialise simulation
		Simulation s = new Simulation(c, env, 500,100);
		// Generate object
//		System.out.println("Starting generation");
		s.generate();
		
		// add printability to overall results
//		res+=(double)(printability);
		
		// calculate convexity
		//Low is good
		float convexity = (1 - s.convexity());
		
		//add convexity to overall results
		res+=(double)(convexity);
		
		// Add angle dispersion to overall result
		// High is good
		float angleDispersion = s.angleDispersionCoefficient();
		res+=(double)(angleDispersion);
		
		res/=2;
		
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
	
	int runs;
	int maxIterations;
	
	Environment e;
	
	FitnessFunction fitfun; //0 = Full fit fun, 1 = printability, 2 = complexity
	
	double[] currentBest;
	
	
	
	public void setEnvironment(Environment _e) {
		this.e = _e;
	}
	
	public double[] getCurrentBest() {
		return this.currentBest;
	}
	
	FitnessFunction selectFitnessFunction(int _f, Environment _e) {
//		FitnessFunction fitfun;
//		if(_f < 0 || _f > 2) {
//			
//		}
		if(_f == 0) return new FullFitness(_e);
		if(_f == 1) return new Printability(_e);
		if(_f == 2) return new Complexity(_e);
		if(_f == 3) return new Convexity(_e);
		if(_f == 4) return new AngleDispersion(_e);
		System.out.println("The value selected for Fitness Function is out of range...\n"
				+ "Select one of the following:\n"
				+ "0.- Full Fitness\n"
				+ "1.- Printability\n"
				+ "2.- Complexity\n"
				+ "3.- Convexity\n"
				+ "4.- Angle Dispersion\n"
				+ "FullFitness is returned as default");
//		return;
		return new FullFitness(_e);
	
		
	}
	
	public void evolveFromParent(double[] parent, int maxIter, int _fitfun, Environment _e, String destinationFolder, String _filePrefix) {
		fitfun = selectFitnessFunction(_fitfun, _e);
		
		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
		cma.setDimension(5); // overwrite some loaded properties
		cma.setInitialX(parent); // in each dimension, also setTypicalX can be used
		cma.setInitialStandardDeviation(0.2); // also a mandatory setting 
//		cma.options.stopFitness = stopFitness;       // optional setting
		cma.options.stopMaxIter = maxIter;
		
		// initialize cma and get fitness array to fill in later
				double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];
				
				//=========FILE OUTPUT===========		
				DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm");
				LocalDateTime now = LocalDateTime.now();
				String date = dtf.format(now);
				String path = destinationFolder;
				String fitnessDataFile = _filePrefix + "_" + date+ "_FitnessData";
				String genesDataFile =  _filePrefix + "_" + date + "_GenesData_";
				Long envRandSeed = _e.getRandomSeed();
						
				// Initialise output data (fitness)
				String fitnessData = "(CMARandomSeed=" + cma.getSeed() + ", ENVRandomSeed=" + envRandSeed + ", "+ new Date().toString() + ")\n";
				// Add header
				fitnessData+="iteration, evaluations, sigma, axisratio, bestever_fitness, best_fitness, median_fitness, worst_fitness, mindii, "
				+ "idxmaxSD, maxSD, idxminSD, minSD \n";
				
				// Initialise output data (genetic information)
				String genesData = "(randomSeed=" + cma.getSeed() + ", ENVRandomSeed=" + envRandSeed + ", " + new Date().toString() + ")\n";
				// Add header
				genesData+="iteration, evaluations, sigma, void, fitness_of_recent_best, x_of_recent_best(1...dimension) \n";

				
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

					// output to files 
					// Add data to data strings
					fitnessData+=cma.getDataRowFitness();
					genesData+=cma.getDataRowXRecentBest();	

					//output to console
					int outmod = 10;
					if (cma.getCountIter() % (15*outmod) == 1)
						cma.printlnAnnotation(); // might write file as well
					if (cma.getCountIter() % outmod == 1) {
						cma.println();
//						bestSolutionPerInterval.add(cma.getBestRecentX());
					}
				
				}
				// evaluate mean value as it is the best estimator for the optimum
				cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 

				// final output		
				// write data strings to output file
				saveData(path, fitnessDataFile, fitnessData);
				saveData(path, genesDataFile, genesData);

				cma.println();
				cma.println("Terminated due to");
				for (String s : cma.stopConditions.getMessages())
					cma.println("  " + s);
				cma.println("best function value " + cma.getBestFunctionValue() 
							+ " at evaluation " + cma.getBestEvaluationNumber());
				currentBest = cma.getBestX();
	
	} //evolveFromParent
		
	public void evolve(int maxIterations, double stopFitness, int _fitfun, int _popSize, Environment _e, String destinationFolder, String _filePrefix) {
//		IObjectiveFunction fitfun = this.fitfun;
		
		fitfun = selectFitnessFunction(_fitfun, _e); 
		
		double[] initialX = new double[5];
		for(int i = 0; i < initialX.length; i++) {
			initialX[i] = new Random().nextDouble();
		}

	
		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
		cma.setDimension(5); // overwrite some loaded properties
//		cma.setInitialX(0.5); // in each dimension, also setTypicalX can be used
		cma.setInitialX(initialX); // in each dimension, also setTypicalX can be used
		cma.setInitialStandardDeviation(0.1); // also a mandatory setting 
		cma.options.stopFitness = stopFitness;       // optional setting
		cma.options.stopMaxIter = maxIterations;
		cma.parameters.setPopulationSize(_popSize);

		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];
		
		//=========FILE OUTPUT===========		
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm");
		LocalDateTime now = LocalDateTime.now();
		String date = dtf.format(now);
		String path = destinationFolder;
		String fitnessDataFile = _filePrefix + "_" + date+ "_FitnessData";
		String genesDataFile =  _filePrefix + "_" + date + "_GenesData_";
		Long envRandSeed = _e.getRandomSeed();
				
		// Initialise output data (fitness)
		String fitnessData = "(CMARandomSeed=" + cma.getSeed() + ", ENVRandomSeed=" + envRandSeed + ", "+ new Date().toString() + ")\n";
		// Add header
		fitnessData+="iteration, evaluations, sigma, axisratio, bestever_fitness, best_fitness, median_fitness, worst_fitness, mindii, "
		+ "idxmaxSD, maxSD, idxminSD, minSD \n";
		
		// Initialise output data (genetic information)
		String genesData = "(randomSeed=" + cma.getSeed() + ", ENVRandomSeed=" + envRandSeed + ", " + new Date().toString() + ")\n";
		// Add header
		genesData+="iteration, evaluations, sigma, void, fitness_of_recent_best, x_of_recent_best(1...dimension) \n";

		
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

			// output to files 
			// Add data to data strings
			fitnessData+=cma.getDataRowFitness();
			genesData+=cma.getDataRowXRecentBest();	

			//output to console
			int outmod = 10;
			if (cma.getCountIter() % (15*outmod) == 1)
				cma.printlnAnnotation(); // might write file as well
			if (cma.getCountIter() % outmod == 1) {
				cma.println();
//				bestSolutionPerInterval.add(cma.getBestRecentX());
			}
		
		}
		// evaluate mean value as it is the best estimator for the optimum
		cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution 

		// final output		
		// write data strings to output file
		saveData(path, fitnessDataFile, fitnessData);
		saveData(path, genesDataFile, genesData);

		cma.println();
		cma.println("Terminated due to");
		for (String s : cma.stopConditions.getMessages())
			cma.println("  " + s);
		cma.println("best function value " + cma.getBestFunctionValue() 
					+ " at evaluation " + cma.getBestEvaluationNumber());
		
		currentBest = cma.getBestX();
				
//		// we might return cma.getBestSolution() or cma.getBestX()
//		cma.println("best solution" + Arrays.toString(cma.getBestX()));
//		System.out.println(bestFunctionPerIter);
//		for(double[] sln : bestSolutionPerInterval)
//			System.out.println("Iteration "+bestSolutionPerInterval.indexOf(sln)+" solution: "+Arrays.toString(sln));
	
	} //evolve
	
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
		
		
	} // save data
	
	
	
	public void evolveMultiple(int _runs, int _gens, int _fitfunType, int _popSize, String _savePath) {
		System.out.println("evolveMultiple has been called to run for "+_runs+" cycles");
		
		for(int i = 0; i < _runs; i++){			
			
			//Set environment
			Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, -1);
			
			//Set file prefix
			String filePrefix = String.format("%03d", i);
			
			//call evolve
			this.evolve(_gens, 1e-6, _fitfunType, _popSize, runEnv, _savePath, filePrefix);
			
		}
		
		
	}
	
	public void evolveMultiple(long[] envSeeds, int _gens, int _fitfunType, int _popSize, String _savePath) {
//		if(envSeeds.length != _runs) {
//			System.out.println("Please make sure that the number of environment seeds ("+envSeeds.length+")\n matches the number of runs ("+_runs+")");
//			return;
//		}
//		
//		System.out.println("evolveMultiple has been called to run for "+_runs+" cycles");
		
		
		for(int i = 0; i < envSeeds.length; i++){			
			
			//Set environment
			Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10,0.2f, envSeeds[i]);
			
			//Set file prefix
			String filePrefix = String.format("%03d", i);
			
			//call evolve
			this.evolve(_gens, 1e-6, _fitfunType, _popSize, runEnv, _savePath, filePrefix);
			
		}
		
		
	}
	
}//class
