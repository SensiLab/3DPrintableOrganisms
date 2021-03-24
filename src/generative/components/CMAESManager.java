package generative.components;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
//import java.util.ArrayList;
//import java.util.Arrays;
//import java.util.Date;
import java.util.Random;

import generative.components.fitness.*;
//import generative.components.fitness.FitnessFunctionTemplate;

import fr.inria.optimization.cmaes.CMAEvolutionStrategy;
import fr.inria.optimization.cmaes.CMAParameters;
import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;
//import generative.components.examples.Printability;
import processing.core.PVector;


public class CMAESManager {
	
	int runs;
	int maxIterations;
	
	public CMAEvolutionStrategy cma;
	
	public double[] fitness;
	
	Environment e;
	int objectSize = 500;
	int objectWarmup = 20;
	
	public FitnessFunction fitfun; //0 = Full fit fun, 1 = printability, 2 = complexity
	
	double[] currentBest;
	
	// =====CONSTRUCTOR=====
	public CMAESManager() {
		cma = new CMAEvolutionStrategy();
		
	}//constructor
	
	public void setFitnessTarget(double _t) {
		this.fitfun.setTarget(_t);
	}
	
	public double getFitnessTarget() {
		double t = this.fitfun.getTarget();
		return t;
	}
	
	public void setEnvironment(Environment _e) {
		this.e = _e;
		this.fitfun.setEnvironment(this.e);
	}
	
	public void setObjectSize(int _size) {
		this.objectSize = _size;
		this.fitfun.setObjectSize(_size);
	}
	
	public void setObjectWarmup(int _warmup) {
		this.objectWarmup = _warmup;
		this.fitfun.setObjectWarmup(_warmup);
	}
	
	public double[] getCurrentBest() {
		return this.currentBest;
	}
	
	public void setFitFun(int _f) {
		switch(_f) {
		case 0:
			this.fitfun = new FullFitness();
			break;
		case 1:
			this.fitfun =  new Printability();
			break;
		case 2:
			this.fitfun =  new Complexity();
			break;
		case 3:
			this.fitfun =  new Convexity();
			break;
		case 4:
			this.fitfun =  new AngleDispersion();
			break;
		case 5:
			this.fitfun =  new Coverage();
			break;
		case 6:
			this.fitfun =  new RefinePrint();
			break;
		default:
			System.out.println("The value selected for Fitness Function is out of range...\n"
					+ "Select one of the following:\n"
					+ "0.- Full Fitness\n"
					+ "1.- Printability\n"
					+ "2.- Complexity\n"
					+ "3.- Convexity\n"
					+ "4.- Angle Dispersion\n"
					+ "5.- Coverage\n"
					+ "6.- RefinePrint\n"
					+ "FullFitness is returned as default");
			this.fitfun = new FullFitness();
			
		}
		
	}
	
	public void setSimLocations(ArrayList<PVector> _locations) {
		this.fitfun.setObjectLocations(_locations);
	}
	
	public void getFitnessFunction() {
		Class ff = this.fitfun.getClass();
		System.out.println("Fitness Function type is "+ff.getName());
	}
	
	public void setFitFunWeights(double pw) {
		Class ff = this.fitfun.getClass();
		if (!ff.getName().contains("Full")) return;
		this.fitfun.setPrintWeight(pw);
		
	}


	/**
	 * Evolve objects from a parent using its genetic information.
	 * @param _fitfun: Integer to select fitness function default: full fitness (0: full fitness, 1: printability, 2: complexity, 3: convexity, 4: angle dispersion)
	 * @param _e: environment to use for simulations
	 * @param _parent: double[] genes of the first individual
	 * @param _mI: int max number of iterations
	 * @param _sF: double stop fitness (disregarded if  < 0);
	 * @param _popSize: (a.k.a. Lambda) int size of the population per generation
	 * @param _mu: int number of parents for recombination
	 * @param _rT: int recombination type default: superlinear (0: equal, 1: linear, 2: superlinear)
	 * @param _ccov: double learning rate
	 * @param destinationFolder
	 * @param _filePrefix
	 */
	public void evolveFromParent(double[] _parent, double _stDev, int _mI, double _sF, int _popSize,int _mu, int _rT, double _ccov, String destinationFolder, String _filePrefix) {
		
//		fitfun = selectFitnessFunction(_fitfun, _e);
		this.getFitnessFunction();
		// new a CMA-ES and set some initial values
		CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
		cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
		cma.setDimension(5); // overwrite some loaded properties
		cma.setInitialStandardDeviation(_stDev); // also a mandatory setting 
		
		cma.setInitialX(_parent); // in each dimension, also setTypicalX can be used
		cma.options.stopMaxIter = _mI;
		cma.options.stopFitness = _sF;       // optional setting
		cma.parameters.setPopulationSize(_popSize);
		
		CMAParameters.RecombinationType recType;
		switch(_rT){
		case 0:
			recType = CMAParameters.RecombinationType.equal;
			break;
		case 1:
			recType = CMAParameters.RecombinationType.linear;
			break;
		case 2:
			recType = CMAParameters.RecombinationType.superlinear;
			break;
		default:
			recType = CMAParameters.RecombinationType.superlinear;
		}
		// _mu sets the number of parents.
		cma.parameters.setRecombination(_mu, recType);
		cma.parameters.setCcov(_ccov);
//		cma.options.stopTolFun = 0.01;

		
		// initialize cma and get fitness array to fill in later
		double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];
				
		//=========FILE OUTPUT===========		
		DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm");
		LocalDateTime now = LocalDateTime.now();
		String date = dtf.format(now);
		String path = destinationFolder;
		String paramLogFile = _filePrefix + "_" + date + "_ParamLog";
		String fitnessDataFile = _filePrefix + "_" + date + "_FitnessData";
		String genesDataFile =  _filePrefix + "_" + date + "_GenesData";
		Long envRandSeed = this.fitfun.envRandomSeed();
			
		// Initialise paramLog
		String paramLog = "parent=" + Arrays.toString(_parent)+"\n"
				+ "CMARandomSeed="+cma.getSeed()+"\n"
				+ "init_st_dev=" + _stDev+"\n"
				+ "step_size=" + _ccov+"\n"
				+ "fitness_target="+ this.fitfun.getTarget() +"\n"
				+ "ENVRandomSeed=" + envRandSeed
				+ e.getParams(); 
		// Initialise output data (fitness)
//		String fitnessData = "parent=" + Arrays.toString(_parent) + ",init_st_dev=" + _stDev + ",step_size=" + _ccov + ",fitness_target="+ this.fitfun.getTarget() +"\n";

//		fitnessData+="(CMARandomSeed="+cma.getSeed() + ",ENVRandomSeed=" + envRandSeed + "," + new Date().toString() + ")\n";
		// Add header
		String fitnessData ="iteration,evaluations,sigma,axisratio,bestever_fitness,best_fitness,median_fitness,worst_fitness,mindii,"
				+ "idxmaxSD,maxSD,idxminSD,minSD\n";
				
		// Initialise output data (genetic information)
//		String genesData = "(randomSeed=" + cma.getSeed() + ",ENVRandomSeed=" + envRandSeed + "," + new Date().toString() + ")\n";
		// Add header
		String genesData ="iteration,evaluations,sigma,void,fitness_of_recent_best,x_of_recent_best(1...dimension)\n";

				
		// iteration loop
		while(cma.stopConditions.getNumber() == 0) {
			
			// --- core iteration step ---
			double[][] pop = cma.samplePopulation(); // get a new population of solutions
			for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
				// a simple way to handle constraints that define a convex feasible domain  
				// (like box constraints, i.e. variable boundaries) via "blind re-sampling" 
				// assumes that the feasible domain is convex, the optimum is  
				while (!this.fitfun.isFeasible(pop[i]))     //   not located on (or very close to) the domain boundary,  
					pop[i] = cma.resampleSingle(i);    //   initialX is feasible and initialStandardDeviations are  
			                                                       //   sufficiently small to prevent quasi-infinite looping here
				// compute fitness/objective value	
				fitness[i] = this.fitfun.valueOf(pop[i]); // fitfun.valueOf() is to be minimized
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
		saveData(path, paramLogFile, paramLog);
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
		
	/**
	 * Evolve objects starting with random genetic information
	 * @param _fitfun: Integer to select fitness function default: full fitness (0: full fitness, 1: printability, 2: complexity, 3: convexity, 4: angle dispersion)
	 * @param _e: environment to use for simulations
	 * @param _mI: int max number of iterations
	 * @param _sF: double stop fitness (disregarded if  < 0);
	 * @param _popSize: (a.k.a. Lambda) int size of the population per generation
	 * @param _mu: int number of parents for recombination
	 * @param _rT: int recombination type default: superlinear (0: equal, 1: linear, 2: superlinear)
	 * @param _ccov: double learning rate (Small _ccov = slow evolution)
	 * @param destinationFolder
	 * @param _filePrefix
	 */
	public void evolve(int _fitfun, Environment _e, double _istd, int _mI, double _sF,  int _popSize, int _mu, int _rT, double _ccov, String destinationFolder, String _filePrefix) {
		IObjectiveFunction fitfun = this.fitfun;
		
//		fitfun = selectFitnessFunction(_fitfun, _e); 
		
		double[] initialX = new double[5];
		
		for(int i = 0; i < initialX.length; i++) {
			double parentMin = 0.45;
			double parentMax = 0.55;
			
			initialX[i] = parentMin + (parentMax-parentMin) * new Random().nextDouble();
		}
		
		this.evolveFromParent(initialX, _istd, _mI, _sF, _popSize, _mu, _rT, _ccov, destinationFolder, _filePrefix);

	} //evolve
	
	static void saveData(String filePath, String fileName, String data) {
		
		// search for file
		
		// Create file
		File output = new File(filePath+fileName+".csv");
		try {
			FileWriter outputWriter = new FileWriter(filePath+fileName+".csv", true);
			outputWriter.write(data);
			outputWriter.close();
		}catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	} // save data
	
	
	
//	public void evolveMultiple(int _runs, int _gens, int _fitfunType, int _popSize, String _savePath) {
//		System.out.println("evolveMultiple has been called to run for "+_runs+" cycles");
//		
//		for(int i = 0; i < _runs; i++){			
//			
//			//Set environment
//			Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 2, 7, 10, 0.1f, -1);
//			
//			//Set file prefix
//			String filePrefix = String.format("%03d", i);
//			
//			//call evolve
//			this.evolve(_fitfunType, runEnv, _gens, -1, _popSize, 4, 1, 0.001, _savePath, filePrefix);
//		}
//		
//		
//	}
	
//	public void evolveMultiple(long[] envSeeds, int _gens, int _fitfunType, int _popSize, String _savePath) {
////		if(envSeeds.length != _runs) {
////			System.out.println("Please make sure that the number of environment seeds ("+envSeeds.length+")\n matches the number of runs ("+_runs+")");
////			return;
////		}
////		
////		System.out.println("evolveMultiple has been called to run for "+_runs+" cycles");
//		
//		
//		for(int i = 0; i < envSeeds.length; i++){			
//			
//			//Set environment
//			Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 2, 7, 10, 0.1f, envSeeds[i]);
//			
//			this.setEnvironment(runEnv);
//			
//			ArrayList<PVector> locations = new ArrayList<PVector>();
//			
////			locations.add(new PVector(runEnv.getWidth() * .33f, runEnv.getHeight() * .33f));
////			locations.add(new PVector(runEnv.getWidth() * .66f, runEnv.getHeight() * .33f));
////			locations.add(new PVector(runEnv.getWidth() * .66f, runEnv.getHeight() * .66f));
////			locations.add(new PVector(runEnv.getWidth() * .33f, runEnv.getHeight() * .66f));
////			
////			this.setSimLocations(locations);
//			
//			//Set file prefix
//			String filePrefix = String.format("%03d", i);
//			
//			//call evolve
//			this.evolve(_fitfunType, runEnv, _gens, -1, _popSize, 4, 1, 0.001, _savePath, filePrefix);
////			this.evolve(_fitfunType, runEnv, _gens, -1, _popSize, 2, 1, 0.0001, _savePath, filePrefix);
//			
//			
//			
//		}
//		
//		
//	}//evolve multiple with environment seeds
	
	public void evolveMultiple(long[] envSeeds, int firstSeed, int seeds, int _gens, int _fitfunType, int _popSize, int _mu,double _istd, double _ccov, String _savePath) {
//		if(envSeeds.length != _runs) {
//			System.out.println("Please make sure that the number of environment seeds ("+envSeeds.length+")\n matches the number of runs ("+_runs+")");
//			return;
//		}
//		
//		System.out.println("evolveMultiple has been called to run for "+_runs+" cycles");
		Environment runEnv = this.e;
		
		for(int i = firstSeed; i < (firstSeed + seeds); i++){			
			
			//Set environment
//			Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 2, 7, 10, 0.1f, envSeeds[i]);
			runEnv.setRandomSeed(envSeeds[i]);
			this.setEnvironment(runEnv);
			
			ArrayList<PVector> locations = new ArrayList<PVector>();
			
//			locations.add(new PVector(runEnv.getWidth() * .33f, runEnv.getHeight() * .33f));
//			locations.add(new PVector(runEnv.getWidth() * .66f, runEnv.getHeight() * .33f));
//			locations.add(new PVector(runEnv.getWidth() * .66f, runEnv.getHeight() * .66f));
//			locations.add(new PVector(runEnv.getWidth() * .33f, runEnv.getHeight() * .66f));
//			
//			this.setSimLocations(locations);
			
			//Set file prefix
			String filePrefix = String.format("%03d", i);
			
			//call evolve
			this.evolve(_fitfunType, runEnv,_istd, _gens, -1, _popSize, _mu, 1, _ccov, _savePath, filePrefix);
//			this.evolve(_fitfunType, runEnv, _gens, -1, _popSize, 2, 1, 0.0001, _savePath, filePrefix);
			
			
			
		}
		
		
	}//evolve multiple with environment seeds
	
	public void refineParent(int _fitfun, Long _envSeed, double[] _parent, int _mI, double _sF, int _popSize,int _mu, int _rT, double _ccov, String destinationFolder, String _filePrefix){
		// check fitness function
		
		//set small stdev as starting point
		
		//set target printability
		
	} //refine parent
	
}//class
