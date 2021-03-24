package generative.components.examples;

import java.util.ArrayList;
//import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import fr.inria.optimization.cmaes.CMAParameters;
import generative.components.*;
//import processing.core.*;
import processing.core.PVector;

public class Testing{
	public static void main(String[] args) {
		long[] newSeeds = new long[20];
		Random rd = new Random();
		for(int i = 0; i < newSeeds.length; i++) {
			newSeeds[i] = rd.nextLong();
		}
//		long newLong = rd.nextLong();
		System.out.println(Arrays.toString(newSeeds));
		
		for(String arg : args) {
			if(arg.charAt(0) == '-') {
				System.out.println(arg);
			}
		}
		
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
//		
		long envSeed = envSeeds[25];
//
//		
////		String folderName = "/home/ubuntu/out/";
		String folderName = "../out/evo_output/prints/";
		
		
		
		double[] genes = {0.7988481744403935,0.6634500073543206,0.4733103519874558,0.8436126954530144,0.6406007942353644};
		Chromosome chr = new Chromosome(genes);
		Environment e = new Environment(600, 600, 15, 0.4f, 5, 2, 5, 5, 0.2f, envSeed);
//		ArrayList<PVector> locations = new ArrayList<PVector>();
//		locations.add(new PVector(.33f*e.getWidth(),.33f*e.getHeight()));
//		locations.add(new PVector(.66f*e.getWidth(),.33f*e.getHeight()));
//		locations.add(new PVector(.66f*e.getWidth(),.66f*e.getHeight()));
//		locations.add(new PVector(.33f*e.getWidth(),.66f*e.getHeight()));
		Simulation s = new Simulation(chr,e,500,20);
//		System.out.println("Simulation has currently "+s.getColonyLayers()+" layers");
//		s.setOutputSize(180);
//		System.out.println("Simulation has currently "+s.getColonyLayers()+" layers");
		s.generate();
////		String fileName = "013_print-best-Coverage001_Run-000_ind-96";
		String fileName = "IEEE_Complexity_025_184_01_v13";
		s.printColony(folderName+fileName);
////		s.saveColony3D(folderName+fileName);
		System.out.println("Testing has run");
		
	}
}
