package generative.components.fitness;

import generative.components.Chromosome;
import generative.components.Environment;
import generative.components.Simulation;

/*
 * Measures Printability
 */
public class Printability extends FitnessFunctionTemplate implements FitnessFunction{
	

	public void setTarget(double _t) {
		this.target = _t;
	}
	
//	public long envRandomSeed() {
//		return this.e.getRandomSeed();
//	}
	
	
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
		
			// Calculate printability
			//High is good so 1 - print is added to result
			float printability = s.printability();
		
			// add printability to overall results
			res+=printability;
		}
		catch(Exception e) {
			System.out.println("Evironment, object size and object warmup values need to be set!");
		}
		
//		return Math.abs(this.target - res);
		return 1 - res;
	}
	
}
