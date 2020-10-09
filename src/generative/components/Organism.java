package generative.components;

import java.io.Serializable;
import java.util.ArrayList;
import java.lang.Math;
import processing.core.PVector;

public class Organism implements Serializable{
	
	Chromosome chromosome;
	
	ArrayList<Cell> cells;
	ArrayList<Spring> springs;
	
	PVector initialLocation;
	
	int initSubdivs;
	float initRadius;
	
	public String splitting;
	float splitThreshold;
	float splitEnergyMult;
	
	public PVector selfIntPoint = null;
    public int[] intersectingSprings = null;
	
	float newOrgEnergy;
	float orgInitEnergy;
	float diffusionRate = 0.1f; //This is currently set as a constant
	float energyLoss = 0.05f;
	int maxSize = 150; 			//The max number of cells 
	
	public Organism() {
		this.cells = new ArrayList<Cell>();
		this.springs = new ArrayList<Spring>();
	}
	
	public Organism(Chromosome chr, PVector _loc, int subdivs, float rad) {
		this.chromosome = chr;
		this.cells = new ArrayList<Cell>();
		this.springs = new ArrayList<Spring>();

		this.initialLocation = _loc;

		this.initSubdivs = subdivs;
		this.initRadius = rad;

		this.splitting = "false";

		this.newOrgEnergy = 1f;//(float) MathLib.mapDouble(this.chromosome.get(4), 0f, 1f, 0.5f, 1f);

		this.init();
	}
	
	public Organism(ArrayList<Cell> _cells) {
		ArrayList<Double> genes = _cells.get(0).chromosome.getGenes();
		this.chromosome = new Chromosome(genes);
		this.cells = _cells;
		this.springs = new ArrayList<Spring>();
		this.splitting = "false";

		this.init();
	}
	
	public Organism(Chromosome _c) {

		this.cells = new ArrayList<Cell>();
		this.springs = new ArrayList<Spring>();

		this.chromosome = _c;

		this.splitting = "false";
		
		this.init();
	}
	
	// COPY CONSTRUCTOR
	public Organism(Organism _o) {
		ArrayList<Double> genes = new ArrayList<Double>(_o.chromosome.getGenes());
		this.chromosome = new Chromosome(genes);
		
		this.cells = new ArrayList<Cell>();
				
		for(Cell c : _o.getCells()) {
			this.cells.add(new Cell(c));
					}
		this.springs = new ArrayList<Spring>();
		
		this.init();
	}
	
	//=======================METHODS==================
	public void addCell(Cell c) {
		this.cells.add(c);
	}
	
	public void addCell(int position, Cell c) {
		this.cells.add(position, c);
	}
	
	void init() {
		if (this.cells.size() < 1) {
			// If the new organism is initialised from scratch
			if (this.initialLocation == null) {
				this.initialLocation = new PVector(0,0);
			}
			for (int i=0; i<this.initSubdivs; i++) {
				float angle = (float)(2*Math.PI)/this.initSubdivs;
				float posX = this.initRadius * (float) Math.cos(angle * i) + this.initialLocation.x;
				float posY = this.initRadius * (float) Math.sin(angle * i) + this.initialLocation.y;
					
				Cell new_cell = new Cell(posX, posY, this.chromosome);
				float _repRad = this.initRadius * (float) Math.sin(angle/2);
				new_cell.repRad = _repRad;
				this.cells.add(new_cell);
				}
//			this.connectCells();
		} else {
			// If the organism is initialised from an ArrayList<Cell>
			for (int i=0; i<this.cells.size() - 1; i++) {
				// measure the repulsion radius of cells based on their current position
				Cell c = this.cells.get(i);
				for (int j = i+1; j < this.cells.size(); j++) {
					Cell other = this.cells.get(j);
					float d = c.loc.dist(other.loc);
					c.repRad = d/2;
				}
			}
//			this.connectCells();
		}
		if (this.springs.size() < this.cells.size()) {
			this.springs.clear();
			this.connectCells();
		}
		this.splitEnergyMult= (float)MathLib.mapDouble(this.chromosome.get(4), 0f, 1f, 0.66f, 1f); 
		this.splitThreshold = this.getMaxEnergy() * this.splitEnergyMult;
		
		orgInitEnergy = getTotalEnergy();
		
	}//init
	
	/**
	 * Creates springs between adjacent cells
	 */
	public void connectCells() {
		if(this.cells.size() < 1) return;
		for (int i=0; i<this.cells.size()-1; i++) {
			Cell sp = this.cells.get(i);
			Cell ep = this.cells.get(i + 1);
			Spring new_s = new Spring(sp, ep);
			new_s.setRestLen(sp.loc.dist(ep.loc));
			// new_s.lightMetabolism = 1 - this.metabolicRate;
			if(this.chromosome != null){
				new_s.setSpringCoef(this.chromosome.get(3));
			}
			this.springs.add(new_s);
		}
		Cell last = this.cells.get(this.cells.size()-1);
		Cell first = this.cells.get(0);
		Spring last_s = new Spring(last, first);
		last_s.setRestLen(last.loc.dist(first.loc));
		this.springs.add(last_s);
	}
	
	public Chromosome getChromosome() {
		return this.chromosome;
	}
	
	public void setCostOfLiving(float c) {
		for (Cell cell : this.cells) {
			cell.setCostOfLiving(c);
		}
	}
	
	public void resetInitEnergy() {
		this.orgInitEnergy = this.getTotalEnergy();
	}
	
	public float getTotalEnergy() {
		float e = 0;
		for (Cell c : this.cells) e+=c.energy;
		return e;
	}
	
	public float getMaxEnergy() {
		float e = 0;
		for(Cell c : this.cells)e+=c.getMaxEnergy();
		return e;
	}
	
	/**
	 * Getter for the arraylist of cells
	 * @return ArrayList<Cell> of cells
	 */
	public ArrayList<Cell> getCells(){
		return this.cells;
	}
	
	/**
	 * Getter for springs
	 * @return ArrayList<Spring> of all the springs
	 */
	public ArrayList<Spring> getSprings(){
		return this.springs;
	}
	
	
	/**
	 * Removes cell and its adjacent springs.
	 * Connects neighbouring cells with a new spring.
	 * @param the index of the cell to be removed
	 */
	public void killCell(int i) {
		int orgSize = this.cells.size();
		if (this.cells.size() < 2) return;

	    if (i >= this.cells.size()) {
	      i = this.cells.size() - 1;
	    }
	    
	    // remove spring after cell
	    this.springs.remove(i);
	    // remove cell
	    this.cells.remove(i);
	    
	    int prev;
	    if (i == 0) {
	      prev = this.cells.size() - 1;
	    } else {
	      prev = i-1;
	    }

	    this.springs.remove(prev);

	    if (prev == this.cells.size()-1) i = 0;

	    float newSpringLen = this.cells.get(prev).loc.dist(this.cells.get(i).loc) * 0.5f;
	    Spring newSpring = new Spring(this.cells.get(prev), this.cells.get(i));
	    newSpring.setRestLen(newSpringLen);
	    this.springs.add(prev, newSpring);
	}
	
	/**
	 * This method populates the selfIntPoint and intersectingSprings fields
	 * It may require rewriting for simplification (perhaps return a PVector point?)
	 */
	void selfIntersection() {
	    if (this.splitting == "self") return;
	    PVector intersection;
	    //check every spring against all other springs
	    for (int i = 0; i < this.springs.size() - 2; i++) {
	      Spring s = this.springs.get(i);
	      for (int j = i + 1; j < this.springs.size(); j++) {
	        Spring other = this.springs.get(j);
	        //if i and j are adjacent, continue
	        if (i == 0 && (j == 1 || j == this.springs.size() -1)) continue;
	        if (i > 0 && (j == i+1 || j == i-1)) continue;
	        //check intersection between springs
	        intersection = s.getIntersectionPoint(other);
	        if (intersection != null) {
	          this.splitting = "self";
	          this.selfIntPoint = intersection;
	          this.intersectingSprings = new int[2];
	          this.intersectingSprings[0] = i;
	          this.intersectingSprings[1] = j;
	        }
	      }
	    }
	  }
	
	/**
	 * checks if the organism should be put in "energy" splitting state
	 */
	public void checkEnergySplitting() {
//	    if (this.getTotalEnergy() >= this.splitThreshold && this.cells.size() > 10) {
		if (this.getTotalEnergy() >= (this.orgInitEnergy * 5f) && this.cells.size() > 10) {
	    	this.splitting = "energy";
	    }
	}
	
	/**
	 * checks if two cells are too close to each other and eliminates one
	 * of them.
	 * Distance threshold is set at 0.001 units.
	 */
	public void overlappingCells(){
	    for (int i = 0; i < this.cells.size()-1; i++) {
	    	Cell c = this.cells.get(i);
	    	for (int j = i+1; j < this.cells.size(); j++) {
	    		Cell other = this.cells.get(j);
	    		if (c.loc.dist(other.loc) <= 0.001) this.killCell(j);
	    	}
	    }
	}
	
	/**
	 * finds the cell with the highest amount of energy
	 * @return the index of the hottest cell
	 */
	public int getHottestCell() {
	    float maxEnergy = 0;
	    int me_index = 0;
	    for (int i = 0; i < this.cells.size(); i++) {
	      if (this.cells.get(i).energy > maxEnergy) {
	        maxEnergy = this.cells.get(i).energy;
	        me_index = i;
	      }
	    }
	    return me_index;
	  }
	
	/**
	 * Diffuse energy from the cell with the most energy in both directions.
	 */
	public void diffuseEnergy() {
		// create an empty array to store the energy that each cell will be
		// receiving
	    ArrayList<Float> cell_energy_add = new ArrayList<Float>();
	    // populate array with zeros
	    for (int i = 0; i < this.cells.size(); i++) {
	      cell_energy_add.add(0f);
	    }
	    
	    
	    int start = this.getHottestCell();
	    int current_cell = start;

	    //go  CCW
	    for (int i = 0; i < this.cells.size()/2; i++) {
	      int next_cell;
	      if (current_cell == this.cells.size() - 1) {
	        next_cell = 0;
	      } else {
	        next_cell = current_cell + 1;
	      }
	      Cell c = this.cells.get(current_cell);
	      Cell n = this.cells.get(next_cell);
	      if (c.energy > n.energy) {
	        float current_energy = cell_energy_add.get(current_cell);
	        float energy_to_pass = (c.energy * this.diffusionRate)/2;

	        current_energy -= energy_to_pass;
	        cell_energy_add.set(current_cell, current_energy);

	        float d = c.loc.dist(n.loc);

	        float next_energy = cell_energy_add.get(next_cell);
	        next_energy += Math.max((energy_to_pass - (this.energyLoss * this.springs.get(current_cell).getLen())),0f);
	        cell_energy_add.set(next_cell, next_energy);
	      }
	      if (current_cell == this.cells.size() - 1) {
	        current_cell = 0;
	      } else {
	        current_cell++;
	      }
	    }

	    //go  CW
	    current_cell =start;
	    for (int i = 0; i < this.cells.size()/2; i++) {
	      int prev_cell;
	      if (current_cell == 0) {
	        prev_cell = this.cells.size() - 1;
	      } else {
	        prev_cell = current_cell - 1;
	      }
	      Cell c = this.cells.get(current_cell);
	      Cell p = this.cells.get(prev_cell);
	      if (c.energy > p.energy) {
	    	  float current_energy = cell_energy_add.get(current_cell);
	        float energy_to_pass = (c.energy * this.diffusionRate)/2;

	        current_energy -= energy_to_pass;
	        cell_energy_add.set(current_cell, current_energy);

	        float d = c.loc.dist(p.loc);

	        float prev_energy = cell_energy_add.get(prev_cell);
	        prev_energy += (energy_to_pass * MathLib.map(d, 0, 300, 1, 0));
	        cell_energy_add.set(prev_cell, prev_energy);
	      }
	      if (current_cell == 0) {
	    	  current_cell = this.cells.size()-1;
	      } else {
	    	  current_cell--;
	      }
	    }
	}
	
	/**
	 * apply attraction force to all the cells in the organism
	 */
	public void applyAttraction() {
		for (int i = 0; i < this.cells.size() - 1; i++) {
			Cell c = this.cells.get(i);
			for (int j = i+1; j < this.cells.size(); j++) {
				Cell other = this.cells.get(j);
				if (c.loc.dist(other.loc) < c.attrRad) c.attract(other);
				if (other.loc.dist(c.loc) < other.attrRad) other.attract(c);
			}
		}
	}
	
	/**
	 * apply repulsion force to all cells in the organism
	 */
	public void applyRepulsion() {
		for (int i = 0; i < this.cells.size() - 1; i++) {
			Cell c = this.cells.get(i);
			for (int j = i+1; j < this.cells.size(); j++) {
				Cell other = this.cells.get(j);
				if (c.loc.dist(other.loc) < c.repRad) c.repel(other);
				if (other.loc.dist(c.loc) < other.repRad) other.repel(c);
			}
		}
	}
	
	public void splitCell(int c) {
		Cell cell = this.cells.get(c);
		Cell newCell = new Cell(cell);
		//split energy of cells in two
		float e = cell.energy;
		cell.energy = e/2;
		newCell.energy = e/2;
		
		//get indices of previous and next
		int prev;
		if(c == 0) {
			prev = this.cells.size()-1;
		}else {
			prev = c-1;
		}
		
		int next;
		if(c == this.cells.size()-1) {
			next = 0;
		}else {
			next = c+1;
		}
		
		//reference to previous and next cells
		Cell prevCell = this.cells.get(prev);
		Cell nextCell = this.cells.get(next);
		
		//get normals of adjacent springs
		PVector normPrev = VectorOps.getNormal(prevCell.loc,  cell.loc);
		PVector normNext = VectorOps.getNormal(cell.loc, nextCell.loc);
		
		//get normal of cell
		PVector norm = normPrev.add(normNext).normalize();
		norm.setMag(e/2);
		
		//rotate normal vector for existing cell
		cell.loc.add(norm.rotate((float)-(Math.PI)/2));
		newCell.loc.add(norm.rotate((float)Math.PI));		
		
		
		//add cell to org
		this.cells.add(c+1, newCell);
		
		//create new springs
		float coef = (float) this.springs.get(c).getSpringCoef();
		float lenMult = this.springs.get(c).newSpringLenMult;
		//remove spring
		this.springs.remove(c);
		Spring s = new Spring(cell, newCell, e/2, coef, lenMult);
		Spring nextS = new Spring(newCell, nextCell, e/2, coef, lenMult);
		
		//add new springs in their corresponding locations
		int sLoc = this.cells.indexOf(cell);
		this.springs.add(sLoc, s);
		int nextSLoc = this.cells.indexOf(nextCell);
		this.springs.add(nextSLoc, nextS);

		

		
		
	} //splitCell
	
	
	/**
	 * Adds a new cell at the midpoint of a spring.
	 * The spring becomes two new ones
	 * @param the index of the spring to split.
	 */
	public void splitSpring(int s) {
		Spring spring = this.springs.get(s);
		PVector mp = VectorOps.getMidPoint(spring.sp.loc, spring.ep.loc);
		this.splitSpring(s, mp);
//	
	}
	
	public void splitSpring(int s, 	PVector split_point) {
		if(s > this.springs.size() - 1) return;
		Cell sp = this.springs.get(s).sp;
	    Cell ep = this.springs.get(s).ep;

	    //reduce energy of sp and ep
	    float sp_transfer = sp.energy / 3;
	    float ep_transfer = ep.energy / 3;

	    sp.energy -= sp_transfer;
	    ep.energy -= ep_transfer;

	    float sc = (float) this.springs.get(s).getSpringCoef();
	    float l =  this.springs.get(s).getLen();
	    float newLenMult = this.springs.get(s).newSpringLenMult;

	    //remove spring
	    this.springs.remove(s);

	    //add cell in the middle
	    PVector mpLoc = split_point;
	    float newSprLen = l * newLenMult;
	    Cell mp = new Cell(mpLoc.x, mpLoc.y, this.chromosome);
	    mp.energy = (sp_transfer = ep_transfer);
	    //mp.setRepRadius(newSprLen/2);

	    //create new springs
	    this.cells.add(s+1, mp);
	    Spring spMp = new Spring(sp, mp, newSprLen, sc, newLenMult);
	    this.springs.add(s, spMp);

	    Spring mpEp = new Spring(mp, ep, newSprLen, sc, newLenMult);
	    this.springs.add(s+1, mpEp);
	}
	
	
	public void setSplitThreshold() {
		this.splitThreshold = this.getMaxEnergy() * this.splitEnergyMult;
	}
	
	void setSpringCoef(float sc) {
		for (Spring s : this.springs) {
			s.setSpringCoef(sc);
		}
	}
	
	public void setSpringSplitThreshold(float t) {
		for(Spring s:this.springs) {
			s.setSplitThreshold(t);
		}
	}
	
	/**
	 * Calculates the perimeter of an organism
	 * @param org
	 * @return
	 */
	public float getPerimeter() {
		float perimeter = 0f;
		for(Spring s : this.springs) perimeter+=s.getLen();
		return perimeter;
	}
	/**
	 * Calculates the area of an organism
	 * @param org
	 * @return
	 */
	public float getArea(){
		
		ArrayList<PVector[]> org = new ArrayList<PVector[]>();
		for(Spring s:this.springs){
			PVector sp = s.sp.loc;
			PVector ep = s.ep.loc;
			PVector[] spring = {sp,ep};
			org.add(spring);
		}

		float psum = 0;
		float nsum = 0;
		for(int i = 0; i< org.size(); i++){
			int sindex = (i+1) % org.size();
			float prod = org.get(i)[0].x * org.get(sindex)[0].y;
			psum+=prod;
		}

		for(int i=0; i<org.size();i++){
			int sindex = (i+1)%org.size();
			float prod = org.get(sindex)[0].x * org.get(i)[0].y;
			nsum += prod;
		}
		float area = (float) Math.abs(0.5*(psum - nsum));
		// println(area);
		return area;
	}
	
	/**
	 * Static method Calculates the area of an organism
	 * @param org
	 * @return
	 */
	public static float getArea(ArrayList<PVector[]> org) {
		float psum = 0;
		float nsum = 0;
		for(int i = 0; i< org.size(); i++){
			int sindex = (i+1) % org.size();
			float prod = org.get(i)[0].x * org.get(sindex)[0].y;
			psum+=prod;
		}

		for(int i=0; i<org.size();i++){
			int sindex = (i+1)%org.size();
			float prod = org.get(sindex)[0].x * org.get(i)[0].y;
			nsum += prod;
		}
		float area = (float) Math.abs(0.5*(psum - nsum));
		// println(area);
		return area;
	}
	
	public static Boolean leftTurn(PVector a, PVector b, PVector c) {
		// println("======== LT execution =============");
		// println("Vectors:", a, b, c);
		PVector p1 = PVector.sub(a, b);
		PVector p2 = PVector.sub(c, b);
		PVector cross = p1.cross(p2);
		// println("Cross.z = "+cross.z);
		// There seems to be a precission issue
		// With angle = 0, cross returns small negative value.
		return cross.z >= -5E-4;
	}
	
	public static PVector bottomMost(ArrayList<PVector> pts) {
		PVector winner = new PVector();
		for (PVector p : pts) {
			if (p.y > winner.y) {
				winner = p;
			}else if(p.y == winner.y){
				if(p.x > winner.x){
					winner = p;
				}
			}else{
				continue;
			}
		}
		return winner;
	}
	
	public static ArrayList<PVector> polarSort(ArrayList<PVector> pts) {
		//make a copy of pts array
		ArrayList<PVector> ptsCopy = new ArrayList<PVector>();
		for (PVector p : pts) ptsCopy.add(new PVector(p.x, p.y));
		// ArrayList<PVector> points = removeDuplicates(ptsCopy);
		//initialise sorted array
		ArrayList<PVector> sorted = new ArrayList<PVector>();

		//get bottomMost
		PVector bm = bottomMost(ptsCopy);

		//add bm to sorted
		sorted.add(bm);

		//remove bm from pts copy
		ptsCopy.remove(bm);
		double prevMinCotan = -1;
		while (ptsCopy.size() > 0) {
			double minCotan = 999999999;
			PVector winner = new PVector();
			for (PVector p : ptsCopy) {
				double c = (p.x - bm.x)/(p.y - bm.y);
				if(c < minCotan){
					minCotan = c;
					winner = p;
				}else if(c == minCotan){
					minCotan = c;
					if(bm.dist(p) <= bm.dist(winner)) winner = p;
				}else{
					continue; 
				}
			}

			float lastX = sorted.get(sorted.size() - 1).x;
			float lastY = sorted.get(sorted.size() - 1).y;
			// if(minCotan == prevMinCotan){
			//   double distToPrev = dist(bm.x, bm.y, lastX, lastY);
			//   double distToCurrent = dist(bm.x, bm.y, winner.x, winner.y);
			//   println("dist to prev = "+distToPrev+"dist to current = "+distToCurrent);
		    // }
			if(winner.x != lastX && winner.y != lastY){
				sorted.add(new PVector(winner.x, winner.y));
				ptsCopy.remove(winner);
			}else{
				ptsCopy.remove(winner);
			}
			prevMinCotan = minCotan;
		}
		return sorted;
	}
	
	public ArrayList<PVector> getConvexHull() {
		// Get all points in organism
		ArrayList<PVector> pts = new ArrayList<PVector>();
		for(Cell c : this.cells) pts.add(c.loc);
		// If an organism is a triangle, return
		if(pts.size() < 4) return pts;
		//sort points
		ArrayList<PVector> sorted = polarSort(pts);

		//traverse pts and check left turns
		ArrayList<Integer> stack = new ArrayList<Integer>();
		stack.add(0);
		stack.add(1);
		int current = stack.get(stack.size()-1);
		int prev = current - 1;
		int next = current + 1;
		while (next!=0) {
			if (leftTurn(sorted.get(prev), sorted.get(current), sorted.get(next))) {
				stack.add(next);
				current = stack.get(stack.size() - 1);
				prev = stack.get(stack.size() - 2);
				if (next == sorted.size() - 1) {
					next = 0;
				} else {
					next++;
				}
			} else {
				//remove last from stack
				stack.remove(stack.size() - 1);
				//set current as last in stack
				current = stack.get(stack.size() - 1);
				//set previous as second to last
				prev = stack.get(stack.size() - 2);

			}
		}

		//initialise convex hull container
		ArrayList<PVector> hull = new ArrayList<PVector>();
	  
		for (int i = 0; i < stack.size(); i++) {
			int index = stack.get(i);			
			PVector p = sorted.get(index);
			hull.add(p);
		}
	  
		return hull;
	}
	
	public float getConvexHullPerimeter() {
		ArrayList<PVector> convexHull = this.getConvexHull();
		float perimeter = 0f;
		for(int i = 0; i < convexHull.size()-1;i++) {
			PVector sp = convexHull.get(i);
			PVector ep = convexHull.get(i+1);
			perimeter+=sp.dist(ep);
		}
		perimeter+=convexHull.get(convexHull.size()-1).dist(convexHull.get(0));
		return perimeter;
	}
	
	public void addEnergy(float e) {
		for (Cell c : this.cells) {
			c.energy+=e/this.cells.size();
		}
	}
	
	/**
	 * Update method for organism
	 */
	public void update(int maxX, int maxY) {
		//check self intersection
		this.selfIntersection();
		//check splitting by energy excess
		this.checkEnergySplitting();
		//check overlapping cells
		this.overlappingCells();
		//apply attraction and repulsion
		this.applyAttraction();
		this.applyRepulsion();
		//update split threshold
		this.splitThreshold = this.getMaxEnergy() * this.splitEnergyMult;
		
		//check splitting state
		if (this.splitting != "false") return;
		
		// energy
		this.diffuseEnergy();

		//update springs
		for (int i = 0; i < this.springs.size(); i++) {
			Spring s = this.springs.get(i);
//	        if ((s.sp.energy + s.ep.energy) >= ((s.sp.maxEnergy + s.ep.maxEnergy) * s.getSplitThreshold()))
//	        	this.splitSpring(i);
	        if (s.getLen() >= (s.sp.maxEnergy + s.ep.maxEnergy) * s.splitThreshold)
	        	this.splitSpring(i);
		}

		for (Spring s : this.springs) {
	        //s.absorbLight(e);
			s.update();
		}
		
		//update cells
		for (int i = 0; i < this.cells.size(); i++) {
			Cell c = this.cells.get(i);
			if (c.energy <= 0.0001f) {
				this.killCell(i);
			}else if(c.energy >= c.getMaxEnergy()){
				//split cell
				this.splitCell(i);
			}else {
		
//				c.borders(maxY, maxY);
		        
				//the maxVel of c is updated based on the size of the organism using a sigmoid function
				float newMaxVel = Math.min(Math.max(MathLib.sigmoid(this.cells.size(), c.maxVel, this.maxSize, -0.1f), c.minVel), c.maxVel);
		        c.setMaxVel(newMaxVel);
		        c.update(maxX, maxY);
			}
		}
		

		

	}
}
