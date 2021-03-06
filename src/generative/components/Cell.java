package generative.components;

import java.lang.Math;

import java.util.ArrayList;

import processing.core.PVector;

public class Cell{
	
	
	public PVector loc;
	public PVector vel;
	public PVector acc;
	
	public Chromosome chromosome;
	
	public float energy = 5f;

	// genetically determined attributes
	double metabolicRate; // ratio of the acquired nutrients that becomes energy
	float maxVel; // represents a drag coefficient: high means 
	float maxEnergy; // this variable is determined by genetic information and it is a proxy for limit mass/surface_area of the cell


	//CONSTANTS
	float costOfLiving = 0.001f;
	int seekDist = 3;
//	public float attrRad = 150f;
	float attrForce = 0.005f;
	float repRad = 5f;
	float minRepRad = 3f;
	float repForce = 0.1f;
	// velocity
	float maxVMin = 0.2f;
	float maxVMax = 2f;
	// energy
	float maxEMin = 6f;
	float maxEMax = 20f;
	
	//-----------CONSTUCTORS
	public Cell(float x, float y, Chromosome _chr) {
		this.acc = new PVector(0, 0);
		this.vel = new PVector(0, 0);
		this.loc = new PVector(x, y);
		
		this.chromosome = _chr;
		
		this.metabolicRate = this.chromosome.get(0);
//		this.maxVel = (float) MathLib.mapDouble(this.chromosome.get(1), 0f, 1f, this.maxVMin, this.maxVMax);
		this.maxVel = (float)this.chromosome.get(1);
		this.maxEnergy = (float) MathLib.mapDouble(this.chromosome.get(2), 0f, 1f, this.maxEMin, this.maxEMax);
	}

	public Cell(PVector _loc, Chromosome _chr) {
	    this.acc = new PVector(0, 0);
	    this.vel = new PVector(0, 0);
	    this.loc = new PVector(_loc.x, _loc.y);

	    this.chromosome = _chr;

	    this.metabolicRate = this.chromosome.get(0);
//	    this.maxVel = (float) MathLib.mapDouble(this.chromosome.get(1), 0f, 1f, this.maxVMin, this.maxVMax);
	    this.maxVel = (float) this.chromosome.get(1);
	    this.maxEnergy = (float) this.chromosome.get(2);
	}


	public Cell(float x, float y, float _attrRad, float _attrForce, float _repRad, float _repForce) {
	    this.acc = new PVector(0, 0);
	    this.vel = new PVector(0, 0);
	    this.loc = new PVector(x, y);

//	    this.attrRad = _attrRad;
	    this.attrForce = _attrForce;
	    this.repRad = _repRad;
	    this.repForce = _repForce;	    
	    
	}
	
	// COPY CONSTRUCTOR
	public Cell(Cell c) {
	    this.acc = c.acc;
	    this.vel = c.vel;
	    float x = c.loc.x;
	    float y = c.loc.y;
	    float z = c.loc.z;
	    this.loc = new PVector(x,y,z);

//	    this.attrRad = c.attrRad;
	    this.attrForce = c.attrForce;
	    this.repRad = c.repRad;
	    this.repForce = c.repForce;

	    this.energy = c.energy;

	    this.chromosome = c.chromosome;

	    this.metabolicRate = c.metabolicRate;
	    this.maxVel = c.maxVel;
	    this.maxEnergy = c.maxEnergy;
	}


	//------------METHODS

	//update
	public void update(Environment _e) {
		//check borders
		this.borders(_e.width,  _e.height);
		
		// update velocity
		this.vel.add(this.acc);
//		this.vel.limit(this.maxVel);

	    //update location
	    this.loc.add(this.vel);

	    //reset acceleration
	    this.acc.mult(0);

	    //update velocity based on weight
	    if(this.vel.mag() > 0) {
	    	//reduce velocity for next time step: depends on current energy.
	    	this.vel.mult(1 - this.energy * 0.0005f);
	    }
	    
//	    float currentVel = Math.max(this.vel.mag(), 0.001f);
	    float currentVel = this.vel.mag();
//	    this.energy-=(this.costOfLiving * this.maxEnergy * this.metabolicRate * currentVel * currentE); // Original formula
	    
//	    this.energy-=((this.costOfLiving) + (this.metabolicRate * 0.01f) + (((this.maxEnergy-this.maxEMin)/(this.maxEMax-this.maxEMin))*0.01f) + (currentVel * 0.01f)); // * (this.maxEMax) * this.metabolicRate * currentVel * currentE);
	    
//	    this.energy-=((this.costOfLiving) + (((this.metabolicRate) * (this.maxEnergy/this.maxEMax)) + ((currentVel/this.maxVMax)*this.costOfLiving)));
//	    this.energy-=(this.costOfLiving + ((this.metabolicRate * this.chromosome.get(2) * (currentVel/1)) * (this.costOfLiving * this.costOfLiving)));
	    this.energy-=this.costOfLiving + (this.metabolicRate + ((currentVel/this.maxVMax)/(this.maxVel*this.maxVel)))*this.costOfLiving;
	    //update split ratio
//	    this.splitThreshold = Math.max((MathLib.sigmoid((this.energy/this.maxEnergy), 1.0f, 0.75f, -15f) * this.maxEnergy), 0.01f);
	} //update

	public void update2() {
		/*
		A special update method that allows cells to move while not consuming
		energy. Used during the organism splitting process.
		*/
		
		// update velocity
		this.vel.add(this.acc);
		this.vel.limit(this.maxVel);
		
		//update location
		this.loc.add(this.vel);
		
		//reset acceleration
		this.acc.mult(0);
		
	    if(this.vel.mag() > 0) {
	    	//reduce velocity for next time step: depends on current energy.
	    	this.vel.mult(1 - this.energy * 0.0001f);
	    }
	} //update2
	
	public void update3() {
		this.loc.add(new PVector(0.1f,0.1f));
	} //update3
	
	public float getRepRadius() {
		return this.repRad;
	}
	
	public float getMaxEnergy() {
		return this.maxEnergy;
	}
	
	public void setAttrForce(float f) {
		this.attrForce = f;
	}
	
	public void setCostOfLiving(float c) {
		this.costOfLiving = c;
	}
	
//	public void setAttrRadius(float r) {
//		this.attrRad = r;
//	}
	
	public void setRepForce(float f) {
		this.repForce = f;
	}
	
	public void setRepRadius(float r) {
		this.repRad = r;
	}
	public void setMaxVel(float v) {
		this.maxVel = v;
	}
	
	public float getMaxVel() {
		return this.maxVel;
	}
	
	void resetMaxVel() {
		this.maxVel = (float) MathLib.mapDouble(this.chromosome.get(1), 0f, 1f, this.maxVMin, this.maxVMax);
	}
	
//	public void applyForce(PVector force) {
////		force.mult(1/(this.maxEnergy/this.maxEMax));
////		force.mult(1 - ((this.maxEnergy - this.maxEMin)/(this.maxEMax - this.maxEMin)));
////		force.mult(1/this.maxEnergy);
//		
//		float forceMag = force.mag();
//		
//		force.normalize();
//		force.mult((forceMag * this.maxVel));
////		force.mult((forceMag * this.maxVel)/this.maxEnergy);
//		
//		this.acc.add(force);
//	}
	
//	public void applyForce(Cell other, PVector force) {
//		//
//	}
	
	public void applyForce(PVector force) {
		// copy force
		PVector f = new PVector(force.x, force.y);
		// dissipate part of the force received because of the drag coefficient
		f.mult(1 - (this.maxVel * (this.maxEnergy/this.maxEMax)));
		
		this.acc.add(f);
	}
	
	public void applyForce2(PVector force) {
		PVector f = new PVector(force.x, force.y);
		this.acc.add(f);
	}

	//attract
	public void attract(Cell other) {
		float attrMag = this.attrForce/(this.loc.dist(other.loc) + 0.0001f);
		PVector attr = PVector.sub(this.loc, other.loc);
	    attr.setMag(attrMag);
	    other.applyForce(attr);
	}

	//repel
	public void repel(Cell other) {
		float repMag = this.repForce/(float) Math.max((double) this.loc.dist(other.loc)/this.repRad, 0.01);
	    PVector rep = PVector.sub(other.loc, this.loc);
	    rep.setMag(repMag);
	    other.applyForce(rep);
	}

	//check drag
	void addDrag(Environment e) {
		float drag = e.drag;
		this.vel.mult(1 - (drag * drag));
	}

//	  //check nutrients
	 public void eat(Environment e) {
	    //float food = e.lookupNutrients(this.loc) * sq(1+this.metabolicRate);
	    float food = e.lookupNutrients(this.loc);
	    if (this.energy < this.maxEnergy) {
	      this.energy+=(food * this.metabolicRate);
	      int[] coord = e.getCurrentTile(this.loc);
	      e.depleteNutrients(coord[0], coord[1]);
	    }
	  }

	  ArrayList<int[]> getNeighbours(Environment e) {
	    ArrayList<int[]> neighbours = new ArrayList<int[]>();
	    int[] here = e.getCurrentTile(this.loc);
	    for (int i = -this.seekDist; i < this.seekDist; i++) {
	      int col = here[0] + i;
	      for (int j = -this.seekDist; j < this.seekDist; j++) {
	        int[] n = new int[2];
	        int row = here[1] + j;
	        n[0] = col;
	        n[1] = row;
	        if ((n[0] >= 0 && n[0] <= e.cols - 1) && (n[1] >= 0 && n[1] <= e.rows - 1)) neighbours.add(n);
	      }
	    }
	    return neighbours;
	  }

	  int[] findMaxEnergyNeighbour(Environment e) {
	    ArrayList<int[]> neighbours = this.getNeighbours(e);
	    float maxFood = 0;
	    int maxFoodTile = 0;
	    for (int i = 0; i < neighbours.size(); i++) {
	      int[] tile = neighbours.get(i);
	      float food = e.lookupNutrients(tile);

	      if (food > maxFood) {
	        maxFood = food;
	        maxFoodTile = i;
	      }
	    }

	    if (maxFood == 0) {
	      return null;
	    } else {
	      return neighbours.get(maxFoodTile);
	    }
	  }

	  PVector findFoodTarget(Environment e) {
	    //ArrayList<int[]> neighbours = this.getNeighbours(e);;
	    int[] targetTile = this.findMaxEnergyNeighbour(e);

	    if (targetTile != null) {
	      float foodX = targetTile[0] * e.resolution + (e.resolution/2);
	      float foodY = targetTile[1] * e.resolution + (e.resolution/2);
	      PVector foodLoc = new PVector(foodX, foodY);
	      return foodLoc;
	    }else{
	     return null;
	    }
	  }
//
//	  //seek food
	  public void seekFood(Environment e) {
	    PVector foodLoc = this.findFoodTarget(e);

	    //apply force
	    if (foodLoc == null) return;
	    if (this.energy >= this.maxEnergy) return;
	    
	    PVector foodTarget = PVector.sub(foodLoc, this.loc);

	    foodTarget.normalize();

	    //a cell will be more capable of moving toward a food source the higher its energy and food metabolism
	    foodTarget.mult(this.energy * (float) this.metabolicRate * 0.1f);
	    this.applyForce(foodTarget);
	    
	    //subtract energy proportional to the magnitude of the force
	    
	    //drawVector(foodTarget, this.loc.x, this.loc.y, this.seekDist * e.resolution);
	  } // seek food


//	  Boolean isOccluded(PVector food, ArrayList<Spring> others) {
//	    //PVector mp = this.getMidpoint();
//	    PVector offset = PVector.sub(food, this.loc).normalize();
//
//	    float x3 = offset.x;
//	    float y3 = offset.y;
//	    float x4 = food.x;
//	    float y4 = food.y;
//
//	    for (Spring other : others) {
//
//	      float x1 = other.sp.loc.x;
//	      float y1 = other.sp.loc.y;
//	      float x2 = other.ep.loc.x;
//	      float y2 = other.ep.loc.y;
//	      float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
//	      if (den == 0) return false;
//
//	      float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
//	      float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;
//	      //println("t: "+t+" u: "+u);
//	      if (t > 1f || t < 0f || u > 1f && u < 0f) return false;
//	    }
//	    //println("true");
//
//	    //draw occlusion lines
//	    //stroke(0);
//	    //strokeWeight(0.1);
//	    //line(x3, y3, x4, y4);
//	    return true;
//	  }

	  public Boolean isInside(Organism o) {
		  float x3 = this.loc.x;
		  float y3 = this.loc.y;

		  float x4 = 9999;	//this is a constant to draw a long horizontal line
		  float y4 = this.loc.y;

		  int intersections = 0;
		  for (Spring s : o.springs) {
			  float x1 = s.sp.loc.x;
			  float y1 = s.sp.loc.y;
			  float x2 = s.ep.loc.x;
			  float y2 = s.ep.loc.y;
			  float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

			  float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
			  float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;

			  if (t <= 1f && t >= 0f && u <= 1f && u >= 0f) intersections++;
		  }
		  return intersections%2 == 1;
	  } // isInside

	  public void borders(int maxX, int maxY) {
	    if (this.loc.x <= 0) {
	      this.loc.x = 0;
	      this.vel.x *= -1;
	    }
	    if (this.loc.x >= maxX) {
	      this.loc.x = maxX;
	      this.vel.x *= -1;
	    }
	    if (this.loc.y <= 0) {
	      this.loc.y = 0;
	      this.vel.y *= -1;
	    }
	    if (this.loc.y >= maxY) {
	      this.loc.y = maxY;
	      this.vel.y *= -1;
	    }
	  } // borders
	  
	  public void noAccess(Environment e) {
		  
	  }

}
