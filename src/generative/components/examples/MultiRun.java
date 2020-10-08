package generative.components.examples;

import generative.components.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;


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
		
		//Set parameters for cma
		//Environment (This environment will be used for all 
//		Environment e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -1);	//Environment with randomisation
//		Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, -1);	//Environment with randomisation
//		cma.setEnvironment(e);
		
		String destFolder = "../out/evo_output/Complexity_004/";
		
		// call multirun
		cma.evolveMultiple(envSeeds, gens, fitfunType, populationSize, destFolder);
//		cma.evolve(10, 1e-6, 1, e, destFolder);		
		
	}
}
