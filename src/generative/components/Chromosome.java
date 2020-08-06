package generative.components;


//import processing.core.*;
import java.util.*;

/**
 * This is a template class and can be used to start a new processing Library.
 * Make sure you rename this class as well as the name of the example package 'template' 
 * to your own Library naming convention.
 * 
 * (the tag example followed by the name of an example included in folder 'examples' will
 * automatically include the example in the javadoc.)
 *
 * @example Hello 
 */

public class Chromosome {
	
	// myParent is a reference to the parent sketch
//	PApplet myParent;
	 //ArrayList<Float> genes;
	float metabolicRate;
	float maxVel;
	float maxEnergy;
	float springCoef;
	float splitRate;

//	int myVariable = 0;
	
	public final static String VERSION = "##library.prettyVersion##";
	

	/**
	 * a Constructor, usually called in the setup() method in your sketch to
	 * initialize and start the Library.
	 * 
	 * @example Hello
	 * @param theParent the parent PApplet
	 */
	public Chromosome(ArrayList<Float> newGenes){
	    this.metabolicRate = newGenes.get(0);
	    this.maxVel = newGenes.get(1);
	    this.maxEnergy = newGenes.get(2);
	    this.springCoef = newGenes.get(3);
	    this.splitRate = newGenes.get(4);
	}
	
	public Chromosome(double[] newGenes){
	    this.metabolicRate = (float)newGenes[0];
	    this.maxVel = (float)newGenes[1];
	    this.maxEnergy = (float)newGenes[2];
	    this.springCoef = (float)newGenes[3];
	    this.splitRate = (float)newGenes[4];
	}
	
	public Chromosome(float _metabolicRate,
			           float _maxVel,
			           float _maxEnergy,
			           float _springCoef,
			           float _splitRate){
		this.metabolicRate = _metabolicRate;
		this.maxVel = _maxVel;
		this.maxEnergy = _maxEnergy;
		this.springCoef = _springCoef;
		this.splitRate = _splitRate;
	}
	
	//=====GET GENES=========
	  public ArrayList<Float> getGenes()
	  {
	     ArrayList<Float> genes = new ArrayList<Float>();
	     genes.add(this.metabolicRate);
	     genes.add(this.maxVel);
	     genes.add(this.maxEnergy);
	     genes.add(this.springCoef);
	     genes.add(this.splitRate);

	     return genes;
	  }

	  public String chrToString(){
	    ArrayList<Float> genes = this.getGenes();
	    ArrayList<String> chr = new ArrayList<String>();
	    for (Float g : genes){
	      chr.add(Float.toString(g));
	    }

	    return String.join("_", chr);

	  }

	  public float get(int index){
	    ArrayList<Float> genes = new ArrayList<Float>();
	    genes.add(this.metabolicRate);
	    genes.add(this.maxVel);
	    genes.add(this.maxEnergy);
	    genes.add(this.springCoef);
	    genes.add(this.splitRate);

	    return genes.get(index);
	  }

}

