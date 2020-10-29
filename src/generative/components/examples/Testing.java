package generative.components.examples;

//import java.util.ArrayList;
import java.util.Arrays;

import fr.inria.optimization.cmaes.CMAParameters;
import generative.components.*;
//import processing.core.*;

public class Testing{
	public static void main(String[] args) {
		String folderName = "../out/evo_output/prints/";
		long envSeed = 8141159607898132170L;
		double[] genes = {0.9629601064711307,0.9810801493988944,0.011216813853811243,0.7389555522538147,0.3225616722119252};
		Chromosome chr = new Chromosome(genes);
		Environment e = new Environment(600, 600, 15, 0.3f, 5, 5, 10, 10, 0.2f, envSeed);
		Simulation s = new Simulation(chr,e,500,20);
		s.generate();
		String fileName = "print-036_Run-003_ind-014";
		s.printColony(folderName+fileName);
		
	}
}
