package generative.components.examples;

import generative.components.*;
import generative.components.fitness.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class SingleRun {
	public static void main(String[] args) {
//		double[] parent = new double[5];
		double[] parent = {0.8032763419526157,0.6776286080330413,0.4836108389762228,0.7634549845544009,0.6767341548828757};
		int gens = 100;
		int fitfunType = 6;
		int populationSize = 40;
		int mu = 4;
		double initStDev = 1e-2;
		double ccov = 1e-6; //step size
		int firstSeed = 25;
		String destFolder = "../out/evo_output/IEEE_Refine_015/";
		
		int objectSize = 500;
		int objectWarmup = 20;
		
		//environment parameters
		int envX = 600;
		int envY = 600;
		int envR = 15;
		float envD = 0.4f;
		int envFS = 5;
		int envFsSizeMax = 5;
		int envFsSizeMin = 2;
		float envFsGrowthRate = 5;
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
		// -ex= environment size in X (default 600)
		// -ey= environment size in Y (default 600)
		// -er= environment resolution (default 15)
		// -ed= environment drag (default 0.4)
		// -ef= Number of food sources (default 5)
		// -em= Max size of food sources (default 5)
		// -en= Min size of food sources (default 2)
		// -eg= food sources growth rate (default 5)
		// -ee= food sources decay rate (default 0.2)
		//
		//split String[]
		for(String s : args) {
			System.out.println(s);
			String[] arg = s.split("=");
			switch(arg[0]) {
				case "-x":
					String parentStr = String.valueOf(arg[1]);
					String[] genesStr = parentStr.replace("(","").replace(")","").split(",");
					if(genesStr.length != parent.length) {
						System.out.println("Please provide parent data using -x= followed by 5 doubles separated by commas.\n"
								+ "Example: -x=(double,double,double,double,double)");
						System.exit(0);
					}
//					ArrayList<Double> genes = new ArrayList<Double>(genesStr.length);
					for(int i = 0; i < genesStr.length; i++) {
						parent[i] = Double.parseDouble(genesStr[i]);
					}
//					parent = new Chromosome(genes);
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
			}//switch

		}//check args
		
		System.out.println("Parent            : "+ Arrays.toString(parent) + "\n"
				+ "Gens              : " + gens +"\n"
				+ "fitfunType        : " + fitfunType+"\n"
				+ "Population        : " + populationSize+"\n"
				+ "Mu (# of parents) : " + mu+"\n"
				+ "Ccov (step size)  : " + ccov+"\n"
				+ "First random seed : " + firstSeed+"\n"
				+ "destination folder: "+ destFolder);

				
		
		//Instantiate the manager
		CMAESManager cma = new CMAESManager();
		
		cma.setFitFun(fitfunType);
		
		
		//set fitness function target (for printability)
		cma.setFitnessTarget(0);
		
		cma.setObjectSize(objectSize);
		
		cma.setObjectWarmup(objectWarmup);
		

		
//		cma.setSimLocations(locations);
		
//		// set environment random seed
//		long envSeed = envSeeds[0];
//		long envSeed = CommonConsts.getEnvSeed(firstSeed);
//		
		//set environment
		Environment e = new Environment();
		e.setRandomSeed(CommonConsts.getEnvSeed(firstSeed));
		cma.setEnvironment(e);
		
		//Get reference value for complexity
		Chromosome c = new Chromosome(parent);
		
		
		
		Simulation s = new Simulation(c, e, objectSize, objectWarmup);
		s.generate();
		
		
		cma.fitfun.setReference(s.complexity());
		cma.fitfun.setTolerance(0.075);
		cma.fitfun.setMaxAttempts(1000);
		
		
		String fitnessFolder = destFolder+"/";
		String filePrefix = String.format("%03d", firstSeed);
		if(!new File(fitnessFolder).isDirectory()) {
			System.out.println("The target directory does not exist");
			return;
		}
		System.out.println("The target directory is ok");
		cma.evolveFromParent(parent, initStDev, gens, -1, populationSize, mu, 2, ccov, fitnessFolder, filePrefix);
	}
}
