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

	
		Chromosome c = new Chromosome(0.272272436f, 0.030558953f, 0.007127416f, 0.001067804f, 0.096573976f);
		Environment env = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, 12345);
		float energy = MathLib.map(c.get(2), 0f, 1f, 5f, 20f);
	//	System.out.println();
		Simulation s = new Simulation(c, env, 500,20);
		long startTime = System.nanoTime();
		s.generate();
		long generateTime = System.nanoTime();
		System.out.println("Generate: "+ ((float)(generateTime - startTime)/1000000));
		float p = s.printability();
		long printabilityTime = System.nanoTime();
		System.out.println("Printability: "+ ((float)(printabilityTime - generateTime)/1000000));
		float conv = s.convexity();
		long convexityTime = System.nanoTime();
		System.out.println("Convexity: "+ ((float)(convexityTime - printabilityTime)/1000000));
		float comp = s.compactness();
		long compactnessTime = System.nanoTime();
		System.out.println("Compactness: "+ ((float)(compactnessTime - convexityTime)/1000000));
		float angles = s.angleDispersionCoefficient();
		System.out.println(angles);
	}
		

	
}
