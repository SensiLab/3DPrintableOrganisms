package generative.components.fitness;

import generative.components.Chromosome;
import generative.components.Environment;
import generative.components.Simulation;

public class RefinePrint extends FitnessFunctionTemplate implements FitnessFunction{
	
	double res = 0;
	int attempts = 0;
	int inds = 0;

	
	public void setTarget(double _t) {
		this.target = _t;
	}
	
	@Override
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		
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
		
		float complexity = s.complexity();
		float printability = s.printability();
		
		
		if(complexity < reference && attempts < this.maxAttempts) {
			// set exit mechanism
			if(attempts == 0) {
				System.out.print("Ind: "+String.format("%03d",this.inds) +" Attempts: "+String.format("%04d", this.attempts+1));
			}else {
				for(int j = 0; j < 4; j++) {
					System.out.print('\b');
				}
				System.out.print(String.format("%04d", this.attempts+1));
			}
			attempts++;
			return false;
		}
		System.out.print('\n');
		this.attempts = 0;
		this.inds++;
		this.res = 1 - printability;
		return true;
	}
	
	@Override
	public double valueOf (double[] x) {
		
		// Initialise result variable
//		double res = 0;
		
		try {
			
		}
		catch(Exception e) {
			System.out.println("Evironment, object size and object warmup values need to be set!");
		}
		
		return this.res;
	}

}
