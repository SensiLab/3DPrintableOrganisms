package generative.components.fitness;

import java.util.ArrayList;

import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;
import generative.components.*;
import processing.core.PVector;

/*
 * Extends the IObjectiveFunction interface to provide readable output
 * for partial fitness metrics
 */
public interface FitnessFunction extends IObjectiveFunction{
	public long envRandomSeed();
	public void setTarget(double _t);
	public void setPrintWeight(double _pw);
	public void setEnvironment(Environment _e);
	public void setObjectSize(int _s);
	public void setObjectWarmup(int _w);
	public void setObjectLocations(ArrayList<PVector> locs);
	public double getTarget();
	public void setReference(double _ref);
	public void setTolerance(double _tol);
	public void setMaxAttempts(int _a);
	
}
