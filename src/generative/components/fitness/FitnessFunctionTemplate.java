package generative.components.fitness;

import java.util.ArrayList;

import generative.components.*;
//import generative.components.fitness.*;
import processing.core.PVector;

public class FitnessFunctionTemplate implements FitnessFunction {
	int objectSize;
	int objectWarmup;
	ArrayList<PVector> objectLocations;
	Environment e;
	double target;
	
	double tolerance = 0.05;
	double reference = 0;
	int maxAttempts = 10;
	
	public void setPrintWeight(double _pw) {
		return;
	}
	
	public void setObjectSize(int _size) {
		this.objectSize = _size;
	}
	
	public void setObjectWarmup(int _warmup) {
		this.objectWarmup = _warmup;
	}
	
	public void setObjectLocations(ArrayList<PVector>  _locs) {
		this.objectLocations = _locs;
	}
	
	public void setEnvironment(Environment _e) {
		this.e = _e;
	}
	
	public long envRandomSeed() {
		return this.e.getRandomSeed();
	}
	
	public void setTarget(double _t) {
		this.target = _t;
	}
	
	public double getTarget() {
		return this.target;
	}
	
	@Override
	public void setReference(double _ref) {
		this.reference = _ref - (_ref * this.tolerance);
	}
	
	@Override
	public void setTolerance(double _tol) {
		this.tolerance = _tol;
		this.setReference(this.reference);
	}
	
	public void setMaxAttempts(int _a) {
		this.maxAttempts = _a;
	}
	
	@Override
	public double valueOf(double[] x) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public boolean isFeasible(double[] x) {
		for(int i = 0; i < x.length; i++) {
			if(x[i] < 0 || x[i] > 1) return false;
		}
		return true;
	}


}
