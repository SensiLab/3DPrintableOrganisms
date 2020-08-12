package generative.components.examples;

import java.util.ArrayList;

import generative.components.*;
//import processing.core.*;

public class Testing{
	public static void main(String[] args) {
//		PApplet.main("Testing");
//		
		
//		
//		
//		
//		long genStart = System.currentTimeMillis();
//		s.generate();
////		s.getColony3D().saveColony("C:\\Users\\camil\\01_Other Work\\Code Dev\\JavaTests\\", "outputTest001");
//		long genDone = System.currentTimeMillis() - genStart;
//		System.out.println("Fitness: generation done ("+genDone+" ms)");
//		System.out.println(s.isComplete());
//		System.out.println(s.convexity());
//		System.out.println(s.printability());
//		System.out.println(s.compactness());
//		
//		Colony3D col3d = s.getColony3D();
//		
////		for (Colony col : col3d.getLayers()) {
////			for
////		}
//		
//		
//	}
	
		//0.272272436f, 0.030558953f, 0.007127416f, 0.001067804f, 0.096573976f
		//0.382983121f, 0.343619075f, 0.455419871f, 0.102930382f, 0.662056974f
		 // //Printability 0.3
		Chromosome c = new Chromosome(0.9948359404002058f, 0.486871528443838f, 0.900167237536065f, 0.984005146251034f, 0.11341161543211142f);
		
		// Printability 0.5
//		Chromosome c = new Chromosome(0.6769059731544329f, 0.64715821613544f, 0.831787344311771f, 0.95897969719842f, 0.10433459776762366f);
		
		// Printability 0.7
//		Chromosome c = new Chromosome(0.8688199862060835f, 0.137167605861153f, 0.807739016705708f, 0.563591815560236f, 0.4605745655195185f);
		
		Environment env = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, 12345);
//		float energy = MathLib.map(c.get(2), 0f, 1f, 5f, 20f);
	//	System.out.println();
		Simulation s = new Simulation(c, env, 500,20);
		long startTime = System.nanoTime();
		s.generate();
		long generateTime = System.nanoTime();
		System.out.println("Generate: "+ ((float)(generateTime - startTime)/1000000));
//		float p = s.printability();
//		long printabilityTime = System.nanoTime();
//		System.out.println("Printability: "+ ((float)(printabilityTime - generateTime)/1000000));
//		float conv = s.convexity();
//		long convexityTime = System.nanoTime();
//		System.out.println("Convexity: "+ ((float)(convexityTime - printabilityTime)/1000000));
//		float comp = s.compactness();
//		long compactnessTime = System.nanoTime();
//		System.out.println("Compactness: "+ ((float)(compactnessTime - convexityTime)/1000000));
//		float angles = s.angleDispersionCoefficient();
//		System.out.println(angles);
		
		String out = "../output/gCode/Printability_03.1";
		s.printColony(out);

	}
		

	
}
