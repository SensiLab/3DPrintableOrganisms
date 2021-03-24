package generative.components.fitness;

import generative.components.Chromosome;
import generative.components.Environment;
import generative.components.Simulation;
//import generative.components.fitness.FitnessFunctionTemplate;

/*
 * Measures fitness based on Printability and Complexity, weighted equally
 */
public class FullFitness extends FitnessFunctionTemplate implements FitnessFunction{
	
	double printWeight = 0.5;
	
	
	public void setPrintWeight(double pw) {
		this.printWeight = pw;
		System.out.println("The fitness function is set to Full Fitness and the printability weight is set to "+this.printWeight);
	}
	
	public double valueOf (double[] x) {
		
		// Initialise result variable
		double res = 0;		
		
		try{			
			// create chromosome from x
			Environment env = new Environment(this.e);
			Chromosome c = new Chromosome(x);
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
			// High is good
			float printability = s.printability();
			
			// add printability to overall results
			res+=((double)(printability) * this.printWeight);
			
			// calculate complexity
			// High is good
			float complexity = s.complexity();
			
			//add complexity to overall results
			res+=((double)(complexity) * (1 - this.printWeight));
			
			// calculate coverage
			//High is good
//			float coverage = s.relativeCoverage();
//			res+=(double)(coverage);
			
			// average the score
			res/=2;
			
			// incorporate the rate of completion
			res = res * s.isComplete();	
			
		}
		catch(Exception e){
			System.out.println("Evironment, object size and object warmup values need to be set!");			
		}
		return 1 - res;
	}
	
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}

//	@Override
//	public void setManager(CMAESManager _m) {
//		// TODO Auto-generated method stub
//		
//	}
	
//	public double getTarget() {
//		
//	}
}
