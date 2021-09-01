package generative.components;


//import processing.core.*;
import java.util.*;

public class Chromosome {
	

	double metabolicRate;
	double maxVel;
	double maxEnergy;
	double springCoef;
	double splitRate;

	/**
	 * Constructor
	 */
	public Chromosome(ArrayList<Double> newGenes){
	    this.metabolicRate = newGenes.get(0);
	    this.maxVel = newGenes.get(1);
	    this.maxEnergy = newGenes.get(2);
	    this.springCoef = newGenes.get(3);
	    this.splitRate = newGenes.get(4);
	}
	
	public Chromosome(double[] newGenes){
	    this.metabolicRate = newGenes[0];
	    this.maxVel = newGenes[1];
	    this.maxEnergy = newGenes[2];
	    this.springCoef = newGenes[3];
	    this.splitRate = newGenes[4];
	}
	
	public Chromosome(double _metabolicRate,
					double _maxVel,
					double _maxEnergy,
					double _springCoef,
					double _splitRate){
		this.metabolicRate = _metabolicRate;
		this.maxVel = _maxVel;
		this.maxEnergy = _maxEnergy;
		this.springCoef = _springCoef;
		this.splitRate = _splitRate;
	}
	
	//=====GET GENES=========
	  public ArrayList<Double> getGenes()
	  {
	     ArrayList<Double> genes = new ArrayList<Double>();
	     genes.add(this.metabolicRate);
	     genes.add(this.maxVel);
	     genes.add(this.maxEnergy);
	     genes.add(this.springCoef);
	     genes.add(this.splitRate);

	     return genes;
	  }
	  
	  public double[] getGenes2() {
		  double[] genes = {(double) this.metabolicRate, (double) this.maxVel, (double) this.maxEnergy, (double) this.springCoef, (double) this.splitRate};
		  return genes;
	  }
	  	
	  public String chrToString(){
	    ArrayList<Double> genes = this.getGenes();
	    ArrayList<String> chr = new ArrayList<String>();
	    for (Double g : genes){
	      chr.add(Double.toString(g));
	    }

	    return String.join(", ", chr);

	  }

	  public double get(int index){
	    ArrayList<Double> genes = new ArrayList<Double>();
	    genes.add(this.metabolicRate);
	    genes.add(this.maxVel);
	    genes.add(this.maxEnergy);
	    genes.add(this.springCoef);
	    genes.add(this.splitRate);

	    return genes.get(index);
	  }

}

