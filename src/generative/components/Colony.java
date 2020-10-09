package generative.components;

import java.io.Serializable;
import java.util.ArrayList;

import java.util.Collections;

import java.lang.Math;
import processing.core.PVector;

public class Colony implements Serializable {
	public ArrayList<Organism> organisms;
	public ArrayList<ArrayList<ArrayList<PVector>>> pointCloud;

	Chromosome chromosome;

	//input parameters
	//organism
	float initSubs;    //initial subdivisions of first organism
	float initRad;     //initial radius of first organism
	
	//spring
	float minRestLen;  //minimum rest length for springs
	float springCoef;  //spring coefficient
	float maxLen;      //maximum spring length (what will it do? create a new cell?)
	float orgRepRadius = 15;

	Environment e;


	  //aux variables
	int first_cell = -1;

	//=========CONSTRUCTORS================
	public Colony() {
		/*
		Generates empty colony
		*/
		this.organisms = new ArrayList<Organism>();
		
	}

	public Colony(Chromosome chr, Environment _e) {
		/*
	    Generates colony with one organism in environment _e;
	    Organism is initalised using values in chromosome chr.
	    */
	    this.chromosome = chr;

	    this.e = _e;

	    this.organisms = new ArrayList<Organism>();
	    // this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    Organism firstOrg = new Organism(this.chromosome, new PVector(this.e.width/2, this.e.height/2), 50, 100);
	    this.organisms.add(firstOrg);
	}

	public Colony(Environment _e, ArrayList<Organism> _orgs) {
	    /*
	    Generates colony in specific environment.
	    */
	    this.organisms = _orgs;
	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    this.e = _e;
	}

	public Colony(Organism o) {
	    this.organisms = new ArrayList<Organism>();
	    this.organisms.add(o);
	    this.chromosome = o.getChromosome();
	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    
	}

	public Colony(ArrayList<Organism> orgs) {
	    this.organisms = orgs;
	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	  }
	
	// COPY CONSTRUCTOR
	public Colony(Colony _c) {
		// Copy chromosome
		ArrayList<Double> genes = new ArrayList<Double>(_c.chromosome.getGenes());
		this.chromosome = new Chromosome(genes);
				
		this.organisms = new ArrayList<Organism>();
		
		//copy organisms
		for (Organism o : _c.organisms) {
			Organism newO = new Organism(o);
			this.organisms.add(newO);
		}
		
	}
		
	// METHODS
	public void setEnvironment(Environment _e) {
		this.e = _e;
	}
	
	public Environment getEnvironment() {
		return this.e;
	}
	
	public void setZPosition(float z) {
		for (Organism o : this.organisms) {
			for(Cell c : o.getCells()) {
				c.loc.z = z;
			}			
		}
	}
	
	public float getTotalEnergy() {
		float e = 0;
		for(Organism o : this.organisms) {
			e+=o.getTotalEnergy();
		}
		return e;
	}
	
	public int getTotalCells() {
		int cells = 0;
		for(Organism o : this.organisms) {
			cells+=o.getCells().size();
		}
		return cells;
	}
	
	public int getTotalSprings() {
		int springs = 0;
		for(Organism o : this.organisms) {
			springs+=o.getSprings().size();
		}
		return springs;
	}
	
	public void setSpringSplitThreshold(float t) {
		for(Organism o:this.organisms) {
			o.setSpringSplitThreshold(t);
		}
	}
	
	/*
	 * Update method
	 */
	public void update() {
		
//		System.out.println("Pre eat");
		this.eat();
		
//		System.out.println("Post eat, pre repulsion");
		this.applyOrgRepulsion();

//		System.out.println("Post repulsion, pre kill");
		//check if ant organism is smaller than 3 cells and kill it
		this.killOrganisms();

//		System.out.println("Post kill, pre overlaps");
	    this.checkFullOverlaps();

	    //this.joinOverlapping(); //This function is to be implemented in the future
	    
//	    System.out.println("Post overlaps, pre update");
	    //update organisms
	    this.updateOrganisms();
	    
//	    System.out.println("Post update");
	} //update colony


	/*
	 * Update organisms
	 */
	void updateOrganisms() {
		for (int i = 0; i< this.organisms.size(); i++) {
			Organism o = this.organisms.get(i);
			if (o.splitting == "false") {
				o.update(e.width, e.height);
			} else if (o.splitting == "energy") {
				this.splitOrganism(i);
			} else if (o.splitting == "self") {
				this.splitSelfIntersect(i);
			}
		}
	} //updateOrganisms
	  
	  void eat() {
		  if (this.e == null) return;
		  if (this.organisms.size() < 1) return;

		  for (Organism o : this.organisms) {
			  if(o.splitting != "false") {
				  return;
			  }else {
				  for (Cell c : o.cells) {
	//				  System.out.println("About to seek food");
					  c.seekFood(this.e);
	//				  System.out.println("Food found! About to eat it");
					  c.eat(this.e);
	//				  System.out.println("Ate the food");
					  c.addDrag(this.e);
			  	}
			  }
		  }
	  }

	  void killOrganisms() {
		  for (int i = 0; i < this.organisms.size(); i++) {
			  Organism o = this.organisms.get(i);
			  if (o.cells.size() < 3) this.organisms.remove(o);
			  if (o.getTotalEnergy() <= 1) this.organisms.remove(o);
		  }
	  }


	  void splitOrganism(int org) {
		  Organism o = this.organisms.get(org);
		  //o.splitting = "energy";

		  //get first point of division
		  int first_cell = o.getHottestCell();

		  //find remaining points
		  int next_cell = first_cell;
		  int last_cell;

		  float split_energy = o.cells.get(first_cell).energy;
		  int dist_counter = 0;
		  for (int i = 0; i < o.cells.size(); i++) {
			  if (split_energy < (o.getTotalEnergy()/2)) {
				  if (next_cell == o.cells.size() - 1) {
					  next_cell = 0;
				  } else {
					  next_cell+=1;
				  }
				  split_energy+=o.cells.get(next_cell).energy;
				  dist_counter++;
			  } else {
				  break;
			  }
		  }
		  
		  
		  if (dist_counter > 3 && o.cells.size() - dist_counter > 3) {
			  last_cell = next_cell;
		  } else {
			  if (first_cell + 4 < o.cells.size()) {
				  last_cell = first_cell+4;
			  } else {
				  last_cell =  first_cell+4 - o.cells.size();
			  }
		  }
		  
		  // Initialise attraction forces between splitting cells
		  PVector force_to_first = PVector.sub(o.cells.get(first_cell).loc, o.cells.get(last_cell).loc);
		  force_to_first.normalize().mult(100);

		  PVector force_to_last = PVector.sub(o.cells.get(last_cell).loc, o.cells.get(first_cell).loc);
		  force_to_last.normalize().mult(100);

		  float split_pts_dist = o.cells.get(first_cell).loc.dist(o.cells.get(last_cell).loc);
		  int cells_to_keep;
		  if (first_cell < last_cell) {
			  cells_to_keep = last_cell - first_cell + 1;
		  } else {
			  cells_to_keep = last_cell + (o.cells.size() - first_cell) + 1;
		  }

		  int cells_to_remove = o.cells.size() - cells_to_keep;
		  int next_to_remove;
		  
		  // split if the distance between split points is equal to or smaller than 2*repulsionRadius
		  if (split_pts_dist <= 2*o.cells.get(first_cell).getRepRadius()) {
			  int new_first = first_cell;
			  Organism new_o = new Organism(this.chromosome);
//			  ArrayList<Cell> newOrgCells = new ArrayList<Cell>();
			  new_o.splitting = "energy";

			  for (int i = 0; i < cells_to_remove; i++) {
				  if (new_first == 0) {
					  next_to_remove = o.cells.size() - 1;
				  } else {
					  next_to_remove = new_first - 1;
					  new_first--;
				  }
				  Cell new_c = new Cell(o.cells.get(next_to_remove));
//				  newOrgCells.add(new_c);
				  new_o.cells.add(new_c);
				  o.killCell(next_to_remove);
			  }
			  o.splitting = "false";

//			  Collections.reverse(newOrgCells);
			  Collections.reverse(new_o.cells);
//			  Organism new_o = new Organism(newOrgCells);
			  new_o.connectCells();
//			  new_o.setSplitThreshold();
			  new_o.setSpringCoef(0.5f);
			  new_o.splitting = "false";
			  o.setSpringCoef(0.5f);
			  o.resetInitEnergy();
			  new_o.resetInitEnergy();
			  this.organisms.add(new_o);
		  } else {
			  //get line between start and end
			  float x1 = o.cells.get(first_cell).loc.x;
			  float y1 = o.cells.get(first_cell).loc.y;
			  float x2 = o.cells.get(last_cell).loc.x;
			  float y2 = o.cells.get(last_cell).loc.y;

			  PVector normal1 = VectorOps.getNormal(x1, y1, x2, y2);
			  PVector normal2 = VectorOps.getNormal(x2, y2, x1, y1);

			  normal1.mult(15);
			  normal2.mult(15);

			  for (int k = 0; k < o.cells.size(); k++) {
				  Cell c = o.cells.get(k);
				  float x = c.loc.x;
				  float y = c.loc.y;

				  float d = (x - x1)*(y2 -y1) - (y -y1)*(x2 - x1);
				  if (first_cell < last_cell) {
					  if (d < 0 && (k > first_cell && k < last_cell)) c.applyForce(normal1);
					  if (d > 0 && (k < first_cell || k > last_cell)) c.applyForce(normal2);
				  } else if (first_cell > last_cell) {
					  if (d < 0 && (k > first_cell || k < last_cell)) c.applyForce(normal1);
					  if (d > 0 && (k < first_cell && k > last_cell)) c.applyForce(normal2);
				  }
			  }

			  o.cells.get(first_cell).applyForce(force_to_last);
			  //apply force_to_last to first_cell's neighbours
			  //clockwise
			  for (int i = 1; i <= cells_to_keep/2; i++) {
				  int cwNeighbour;
				  if (first_cell + i > o.cells.size() - 1) {
					  cwNeighbour =  (first_cell + i) - o.cells.size();
				  } else {
					  cwNeighbour = first_cell + i;
				  }
				  PVector f = force_to_last.mult((cells_to_keep/2-i)/(cells_to_keep/2));
				  o.cells.get(cwNeighbour).applyForce(f);
			  }
			  //counterclockwise
			  for (int i = 1; i <= cells_to_remove/2; i++) {
				  int ccwNeighbour;
				  if (first_cell - i < 0) {
					  ccwNeighbour = o.cells.size() - i;
				  } else {
					  ccwNeighbour = first_cell - i;
				  }
				  PVector f = force_to_last.mult((cells_to_remove/2-i)/(cells_to_remove/2));
				  o.cells.get(ccwNeighbour).applyForce(f);
			  }

			  o.cells.get(last_cell).applyForce(force_to_first);
			  //apply force_to_last to first_cell's neighbours
			  //clockwise
			  for (int i = 1; i <= cells_to_remove/2; i++) {
				  int cwNeighbour;
				  if (last_cell + i > o.cells.size() - 1) {
					  cwNeighbour =  (last_cell + i) - o.cells.size();
				  } else {
					  cwNeighbour = last_cell + i;
				  }
				  PVector f = force_to_first.mult((cells_to_remove/2-i)/(cells_to_remove/2));
				  o.cells.get(cwNeighbour).applyForce(f);
			  }
			  //counterclockwise
			  for (int i = 1; i <= cells_to_keep/2; i++) {
				  int ccwNeighbour;
				  if (last_cell - i < 0) {
					  ccwNeighbour = o.cells.size() - i;
				  } else {
					  ccwNeighbour = last_cell - i;
				  }
				  PVector f = force_to_first.mult((cells_to_keep/2-i)/(cells_to_keep/2));
				  o.cells.get(ccwNeighbour).applyForce(f);
			  }

			  o.applyRepulsion();
			  o.applyAttraction();
		  }

		  for (Cell c : o.cells) {
			  c.update2();
		  }
		  for (Spring s : o.springs) {
			  s.update();
		  }
	  }

	  void splitSelfIntersect(int org) {
		  if (this.organisms.get(org).cells.size() <= 4){
			  this.organisms.remove(org);
			  return;
		  }
		  Organism o = this.organisms.get(org);
		  o.splitting = "self";

		  //Point of intersection
		  PVector intPoint = o.selfIntPoint;

		  //Intersecting springs
		  int selfInt1 = o.intersectingSprings[0];
		  int selfInt2 = o.intersectingSprings[1];

		  //Add intersection point for first spring
		  o.splitSpring(selfInt1, intPoint);

		  selfInt2+=1;

		  //Add intersection point for second spring
		  o.splitSpring(selfInt2, intPoint);

		  //get first point of division
		  int first_cell = selfInt1 + 1;
		  int last_cell = selfInt2 + 1;


		  int cells_to_keep = last_cell  - first_cell;

		  int cells_to_remove = o.cells.size() - cells_to_keep;
		  int next_to_remove;

		  int new_first = first_cell;

		  Organism new_o = new Organism(this.chromosome);
		  new_o.splitting = "self";

		  for (int i = 0; i < cells_to_remove; i++) {
			  if (new_first == 0) {
				  next_to_remove = o.cells.size() - 1;
			  } else {
				  next_to_remove = new_first - 1;
				  new_first--;
			  }
			  Cell new_c = new Cell(o.cells.get(next_to_remove));
			  // new_c.metabolicRate = this.chromosome.get(0);
			  // new_c.maxVel = this.chromosome.get(1);
			  // new_c.maxEnergy = this.chromosome.get(2);
			  new_o.cells.add(new_c);
			  o.killCell(next_to_remove);
		  }

		  Collections.reverse(new_o.cells);
		  new_o.connectCells();
		  new_o.setSplitThreshold();
		  new_o.setSpringCoef((float)this.chromosome.get(3));
		  new_o.splitting = "false";
	    
		  this.organisms.add(new_o);
		  if(new_o.getArea() < o.getArea() / 100) this.organisms.remove(new_o);
		  if(o.getArea() < new_o.getArea() / 100) this.organisms.remove(o);

		  o.splitting = "false";
	  }

	  //======================OVERLAP FUNCTIONS====================
	  //Check
	  public void checkFullOverlaps() {
		  if (this.organisms.size() < 2) return;
		  for (int i = 0; i < this.organisms.size() -1; i++) {
			  Organism o = this.organisms.get(i);
			  for (int j = i+1; j < this.organisms.size(); j++) {
				  Organism other = this.organisms.get(j);
				  int cellsInside = 0;
				  for (Cell c : other.cells) {
					  if (c.isInside(o)) cellsInside++;
				  }
				  if (cellsInside == other.cells.size()) {
					  o.addEnergy(other.getTotalEnergy());
					  this.organisms.remove(other);
				  }

				  cellsInside = 0;
				  for (Cell c : o.cells) {
					  if (c.isInside(other)) cellsInside++;
				  }
				  if (cellsInside == o.cells.size()) {
					  other.addEnergy(o.getTotalEnergy());
					  this.organisms.remove(o);
				  }
			  }
		  }
	  }

	  //Execute
	  public void joinOverlapping() {
		  if (this.organisms.size() < 2) return;
		  for (int i = 0; i < this.organisms.size() -1; i++) {
			  Organism o = this.organisms.get(i);
			  for (int j = i+1; j < this.organisms.size(); j++) {
				  Organism other = this.organisms.get(j);
				  ArrayList<Intersection> intersections = findIntersections(o, other);
				  if (intersections.size() > 0) {
					  
					  //add all cells not in intersections to arraylist
					  ArrayList<Cell> newOrg = new ArrayList<Cell>();
					  for(Cell c : o.cells) {
						  if(!c.isInside(other)) newOrg.add(c);
					  }
					  
					  for(Cell c : other.cells) {
						  if(!c.isInside(o)) newOrg.add(c);
					  }
					  
					  //polar sort cells
					  
					  
//					  o.splitting = "overlapping";
//					  other.splitting = "overlapping";
//
//					  Spring first_o = intersections.get(0).a;
//					  int index_first_o = o.springs.indexOf(first_o);
//					  Spring last_other = intersections.get(0).b;
//					  int index_last_other = other.springs.indexOf(last_other);
//					  PVector firstInt = first_o.getIntersectionPoint(last_other);
//
//					  o.splitSpring(index_first_o, firstInt);
//					  //other.splitSpring(index_last_other, firstInt);
//
//					  Spring last_o = intersections.get(intersections.size() - 1).a;
//					  int index_last_o = o.springs.indexOf(last_o);
//					  Spring first_other = intersections.get(intersections.size() - 1).b;
//					  int index_first_other = other.springs.indexOf(first_other);
//					  PVector lastInt = last_o.getIntersectionPoint(first_other);
//
//					  o.splitSpring(index_last_o, lastInt);
//					  other.splitSpring(index_first_other, lastInt);
//
//					  ArrayList<Cell> newCells = new ArrayList<Cell>();
//
//					  for (int k = 0; k < o.cells.size(); k++) {
//						  Cell c = o.cells.get(k);
//						  if (c.isInside(other)) {
//							  o.killCell(k);
//						  }
//					  }
//					  for (int k = 0; k < other.cells.size(); k++) {
//						  Cell c = other.cells.get(k);
//						  if (c.isInside(o)) {
//							  other.killCell(k);
//						  }
//					  }
				  }
			  }
		  }
	  }


	  public ArrayList<Intersection> findIntersections(Organism a, Organism b) {
		  ArrayList<Intersection> foundIntersections = new ArrayList<Intersection>();
		  //find intersections
		  for (Spring sa : a.springs) {
			  for (Spring sb : b.springs) {
				  if (sa.getIntersectionPoint(sb) != null) foundIntersections.add(new Intersection(sa, sb));
			  }
		  }
		  return foundIntersections;
	  }

	  void applyAttraction() {
	    ArrayList<Cell> allCells = new ArrayList<Cell>();
	    for (Organism o : this.organisms) {
	      allCells.addAll(o.cells);
	    }
	    for (int i = 0; i < allCells.size() - 1; i++) {
	      Cell c = allCells.get(i);
	      for (int j = i+1; j < allCells.size(); j++) {
	        Cell other = allCells.get(j);
	        if (c.loc.dist(other.loc) < c.attrRad) c.attract(other);
	        if (other.loc.dist(c.loc) < other.attrRad) other.attract(c);
	      }
	    }
	  }

	  void applyRepulsion() {
	    ArrayList<Cell> allCells = new ArrayList<Cell>();
	    for (Organism o : this.organisms) {
	      allCells.addAll(o.cells);
	    }
	    for (int i = 0; i < allCells.size() - 1; i++) {
	      Cell c = allCells.get(i);
	      for (int j = i+1; j < allCells.size(); j++) {
	        Cell other = allCells.get(j);
	        if (c.loc.dist(other.loc) < c.repRad) c.repel(other);
	        if (other.loc.dist(c.loc) < other.repRad) other.repel(c);
	      }
	    }
	  }

	  void applyOrgRepulsion() {
	    if (this.organisms.size() < 2) return;
	    for (int i = 0; i < this.organisms.size() - 1; i++) {
	      Organism o = this.organisms.get(i);
	      for (Cell c : o.cells) {
	        for (int j = i+1; j < this.organisms.size(); j++) {
	          Organism other = this.organisms.get(j);
	          for (Cell otherC : other.cells) {
	            if (c.loc.dist(otherC.loc) < this.orgRepRadius) c.repel(otherC);
	            if (otherC.loc.dist(c.loc) < this.orgRepRadius) otherC.repel(c);
	          }
	        }
	      }
	    }
	  }
	  
	  public void setCostOfLiving(float c) {
		  for (Organism org : this.organisms) {
			  org.setCostOfLiving(c);
		  }
	  }

	  int getCellCount() {
	    int totalCells = 0;
	    for (Organism o : this.organisms) {
	      totalCells += o.cells.size();
	    }
	    return totalCells;
	  }
	  
	  public void setSpringRLMult(float rl) {
		  for (Organism org : this.organisms) {
			  for (Spring s : org.getSprings()) {
				  s.setEnergyMult(rl);
			  }
		  }
	  }
	  
	  public Chromosome getChromosome() {
		  return this.chromosome;
	  }

}
