package generative.components.examples;

//import java.util.ArrayList;
import java.util.Arrays;

import fr.inria.optimization.cmaes.CMAParameters;
import generative.components.*;
//import processing.core.*;

public class Testing{
	public static void main(String[] args) {
		
		long[] envSeeds = {-1632478226818474777L, 8558593120860338761L, 9063855533770827947L, 8141159607898132170L, -2228090770267545711L,
				5760679812470627183L, 8122425018693031776L, 7233387513096301004L, 8063979647430030046L, 1209585713320896596L};
		
		long envSeed = envSeeds[0];
		
		String folderName = "../out/evo_output/prints/";
		double[] genes = {0.9983524981822078,0.9236473791515063,0.6120141936895009,0.37475220672842674,0.5075528086203049};
		Chromosome chr = new Chromosome(genes);
		Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, envSeed);
		Simulation s = new Simulation(chr,e,500,20);
		s.generate();
		String fileName = "print-best-IndRefinement01_Run-003_ind-262";
		s.printColony(folderName+fileName);
		
	}
}
