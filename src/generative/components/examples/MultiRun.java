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
		int populationSize = 8;
		
		long[] envSeeds = {-1632478226818474777L, 8558593120860338761L, 9063855533770827947L, 8141159607898132170L, -2228090770267545711L,
				5760679812470627183L, 8122425018693031776L, 7233387513096301004L, 8063979647430030046L, 1209585713320896596L};
		
		
		// Process inputs
		// -r= for number of runs
		// -g= for number of gens per run
		// -f= for fitness function
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
				case "-p":
					populationSize = Integer.valueOf(arg[1]);
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
						+ "Population: " + populationSize);
		
		
		
		//Instantiate the manager
		CMAESManager cma = new CMAESManager();
		
		//Set fitness function
		cma.setFitFun(1);
		
		//set fitness function target (for printability)
		cma.setFitnessTarget(0.95);
		
		cma.setObjectSize(500);
		
		cma.setObjectWarmup(20);
		
		// set environment random seed
		long envSeed = envSeeds[0];
		
		//set environment
		Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, envSeed);
		
		cma.setEnvironment(e);
		
		
		
		
		// parameters for evolve from parent
						
		// define parent
		double[] parent = {0.9981469087156556, 0.9234891616280703, 0.61189524355323, 0.3746596492759553, 0.5079931545386839};
		
		// initial stDev
		double stDev = 1e-4;
		
		// gens
		int maxIter = 400;
		
		// lambda
		int lambda = 20;
		
		// mu
		int mu = 2;
		
		//step size
		double stepSize = 1e-6;
		
		//set destination folder and file prefix
		String folderPath = "../out/evo_output/Ind_Refinement_001/";
						
		//set file prefix
		String filePrefix = "IndRefinement003";
		
		//call evolve from parent
		cma.evolveFromParent(parent, stDev, maxIter, -1, lambda, mu, 1, stepSize, folderPath, filePrefix);
		
		


		
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
//		long[] envSeedsFF = Arrays.copyOfRange(envSeeds, 8, 10);
//		String fitnessFolder = "../out/CMA_Testing/Full_Fitness_008/";
//		cma.evolveMultiple(envSeedsFF, 250, 0, 40, fitnessFolder);
//		
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
