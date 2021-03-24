package generative.components.fitness;

import generative.components.Chromosome;
import generative.components.Environment;
import generative.components.Simulation;
//import generative.components.fitness.*;

public class Convexity extends FitnessFunctionTemplate implements FitnessFunction{
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;
		
		try {
		
			// create chromosome from x
			Environment env = new Environment(this.e);
			Chromosome c = new Chromosome(x);
	//		
			// Initialise simulation
			Simulation s;
			// Initialise simulation
			if(this.objectLocations == null) {
				s = new Simulation(c, env, this.objectSize,this.objectWarmup);
			}else {
				s = new Simulation(c, env, this.objectSize,this.objectWarmup, this.objectLocations);
			}
			// Generate object
			s.generate();
			
			// calculate convexity
			//Low is good
			res += s.convexity();
		}
		catch(Exception e) {
			System.out.println("Evironment, object size and object warmup values need to be set!");
		}
		
		return res;
		
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}
}
