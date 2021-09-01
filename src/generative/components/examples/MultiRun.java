package generative.components.examples;

import generative.components.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.Random;


public class MultiRun {
	public static void main(String[] args) {
		int runs = 0;
		int gens = 0;
		int fitfunType = 0;
		double printWeight = 0.5;
		int populationSize = 8;
		int mu = 4;
		double initStDev = 0.1;
		double ccov = 1e-3; //step size
		int firstSeed = 6;
		String destFolder = "";
		
		//environment parameters
		int envX = 600;
		int envY = 600;
		int envR = 15;
		float envD = 0.3f;
		int envFS = 5;
		int envFsSizeMax = 5;
		int envFsSizeMin = 2;
		float envFsGrowthRate = 3;
		float envFsDecayRate = 0.2f;
		
		
		long[] envSeeds = {-1632478226818474777L, 8558593120860338761L, 9063855533770827947L, 8141159607898132170L, -2228090770267545711L,
				5760679812470627183L, 8122425018693031776L, 7233387513096301004L, 8063979647430030046L, 1209585713320896596L, 
				724224084362998227L, -5071659569307755225L, 5823929498051207574L, -551995978363359132L, -5192310287763653706L,
				5108100186193579751L, 4676335984590537776L, 7191218674032234786L, 7152999658342319501L, 7062340219330815167L,
				4871260159946803014L, -5288381650226216036L, -1164060861243111837L, -7174953927911303402L, -9094118736910298012L,
				-3497492774682784886L, -1682527197870981019L, 476367449420437529L, -2508839713788940676L, -6172147982322942294L, 
				-8486584721678693775L, 7170492297316362611L, -7692137076181517887L, 3861397882188512846L, -3110163977026367975L,
				-8634546392007882087L, 7557350025982015301L, 4025513581029972815L, -1483865547251428910L, -309284057398871553L,
				-186292376937715704L, -4463713339575002910L, 7878392066209648151L, -4195456362027586512L, 1429266236894580111L,
				3045732330329506990L, 4597571254417314837L, 5927936720052420491L, -1593416348006868653L, 3027521254663839597L};
		
		
		// Process inputs
		// -r= for number of runs (number of different environments to be used)
		// -g= for number of gens per run
		// -f= for fitness function (0 = full fitness, 1 = ...)
		// -p= for lambda (default should be 40)
		// -mu= for mu (numebr of parents
		// -ss= step size (small means slower progress but less overshooting)
		// -fs= first seed in array of random seeds for environment
		// -df= destination folder for the .csv files produced by the evo process (/home/ubuntu/out/<name of folder>_ underscore at the end!)
		// 
		//
		//
		//
		// -ef= Number of food sources
		// -em= Max size of food sources
		// -en= Min size of food sources
		// -eg= food sources growth rate
		// -ed= food sources decay rate
		//
		//split String[]
		for(String s : args) {
			System.out.println(s);
			String[] arg = s.split("=");
			switch(arg[0]) {
				case "-r":
					runs = Integer.valueOf(arg[1]);
					break;
				case "-g":
					gens = Integer.valueOf(arg[1]);
					break;
				case "-f":
					fitfunType = Integer.valueOf(arg[1]);
					break;
				case "-pw":
					printWeight = Double.valueOf(arg[1]);
					break;
				case "-p":
					populationSize = Integer.valueOf(arg[1]);
					break;
				case "-mu":
					mu = Integer.valueOf(arg[1]);
					break;
				case "-fs":
					firstSeed = Integer.valueOf(arg[1]);
					break;
				case "-df":
					destFolder = arg[1];
					break;
				case "-sd":
					initStDev = Double.valueOf(arg[1]);
				case "-ss":
					ccov = Double.valueOf(arg[1]);
					break;
				case "-ex":
					envX = Integer.valueOf(arg[1]);
					break;
				case "-ey":
					envY = Integer.valueOf(arg[1]);
					break;
				case "-er":
					envR = Integer.valueOf(arg[1]);
					break;
				case "-ed":
					envD = Float.valueOf(arg[1]);
					break;
				case "-ef":
					envFS = Integer.valueOf(arg[1]);
					break;
				case "-em":
					envFsSizeMax = Integer.valueOf(arg[1]);
					break;
				case "-en":
					envFsSizeMin = Integer.valueOf(arg[1]);
					break;
				case "-eg":
					envFsGrowthRate = Float.valueOf(arg[1]);
					break;
				case "-ee":
					envFsDecayRate = Float.valueOf(arg[1]);
					break;	
				default:
					break;
			}
//			for(String st : arg) System.out.println(st);
//			if(arg[0] == "-r") {
//				
//				runs = Integer.valueOf(arg[1]);
//			}
//			if(arg[0] == "-g") gens = Integer.parseInt(arg[1]);
//			if(arg[0] == "-f") fitfunType = Integer.parseInt(arg[1]);
		}
		
		System.out.println("Runs: "+ runs + "\n"
						+ "Gens: " + gens +"\n"
						+ "fitfunType: " + fitfunType+"\n"
						+ "Population: " + populationSize+"\n"
						+ "Mu (# of parents): " + mu+"\n"
						+ "Ccov (step size): " + ccov+"\n"
						+ "First random seed: " + firstSeed+"\n"
						+ "destination folder: "+ destFolder);
		
		
		//Instantiate the manager
		CMAESManager cma = new CMAESManager();
		
		/* Set fitness function
		 * "0.- Full Fitness (Printability + Complexity / 2"
		 * "1.- Printability"
		 * "2.- Complexity"
		 * "3.- Convexity"
		 * "4.- Angle Dispersion"
		 * "5.- Coverage"
		 */
		cma.setFitFun(fitfunType);
		
		if(fitfunType == 0) {
			cma.setFitFunWeights(printWeight);
		}
		
		//set fitness function target (for printability)
		cma.setFitnessTarget(0);
		
		cma.setObjectSize(500);
		
		cma.setObjectWarmup(20);
		
//		cma.setSimLocations(locations);
		
//		// set environment random seed
//		long envSeed = envSeeds[0];
//		
		//set environment
		Environment e = new Environment();
		
		cma.setEnvironment(e);
		
		
		
		
		// parameters for evolve from parent
						
		// define parent
//		double[] parent = {0.9981469087156556, 0.9234891616280703, 0.61189524355323, 0.3746596492759553, 0.5079931545386839};
		
		// initial stDev
		double stDev = 1e-1;
		
//		 gens
//		int maxIter = 500;
		
		// lambda
		int lambda = 40;
		
//		// mu
//		int mu = 4;
		
		//step size
		double stepSize = 1e-3;
		
		//set destination folder and file prefix
//		String folderPath = "../out/evo_output/Ind_Refinement_001/";
						
		//set file prefix
//		String filePrefix = "IndRefinement003";
		
		//call evolve from parent
//		cma.evolveFromParent(parent, stDev, maxIter, -1, lambda, mu, 1, stepSize, folderPath, filePrefix);
		
		


		
		//Set parameters for cma
		//Environment (This environment will be used for all 
//		Environment e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);	//Environment with randomisation
//		Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, -1);	//Environment with randomisation
//		cma.setEnvironment(e);
		
//		String destFolder = "../out/evo_output/Complexity_010/";
		
		// call multirun
//		int[] popSizes = {10, 20, 40, 80};
//		int[] mus = {0, 1, 2, 4};
//		int[] rTypes = {0,1,2};

//		Environment runEnv = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, 8141159607898132170L);
		
//		for(int popSize : popSizes) {
//			for(int mu : mus) {
//				for(int rType : rTypes) {
//					String prefix = "PopSize_"+popSize+"_Mu_"+mu+"_rType_"+rType+"_";
//					cma.evolve(fitfunType, runEnv, gens, 0, popSize, mu, rType, 0.001, destFolder, prefix);
//				}
//			}
//		}
		
		//run full fitness
//		long[] envSeedsFF = Arrays.copyOfRange(envSeeds,0 , 1);
		
//		int firstSeed = 0;
//		int seeds = 6;
//		String fitnessFolder = destFolder+firstSeed+"-"+(firstSeed + runs - 1)+"/";
		String fitnessFolder = destFolder+"/";
		cma.evolveMultiple(envSeeds, firstSeed, runs, gens, fitfunType, populationSize, mu, initStDev, ccov, fitnessFolder);
//		System.out.println(fitnessFolder);
		//run complexity
//		long[] envSeedsC = Arrays.copyOfRange(envSeeds, 8, 10);
//		String complexityFolder = "../out/CMA_Testing/Complexity_018/";
//		cma.evolveMultiple(envSeeds, 250, 2, 40, complexityFolder);
//		cma.evolve(10, 1e-6, 1, e, destFolder);	
//		int totalChunks = 5;
//		int chunkSize = 50;
//		for(int seedIndex = 0; seedIndex < envSeeds.length; seedIndex++) {
//			String filePrefix = String.format("%03d", seedIndex);
//			CMAESManager manager = new CMAESManager();
//			long seed = envSeeds[seedIndex];
//			Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, seed);
//			for(int chunks = 0; chunks < totalChunks;chunks++) {
//				double[] parent = new double[5];
//				if(chunks == 0) {
//					for(int i = 0; i < parent.length; i++) {
//						parent[i] = new Random().nextDouble();
//					}
//				}else {
//					parent = manager.getCurrentBest();
//				}
//				
////				manager.evolveFromParent(2, 0.1, e, parent, chunkSize, -1, 40, 2, 1, 0.001, complexityFolder, filePrefix);
//				//this.evolve(_fitfunType, runEnv, _gens, -1, _popSize, 2, 1, 0.001, _savePath, filePrefix);
//			}
//			
//		}
		
	}
}
