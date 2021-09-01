package generative.components;

import java.util.ArrayList;

import java.util.Collections;

//import java.lang.Math;
import processing.core.PVector;

public class Colony {
	public ArrayList<Organism> organisms;
//	public ArrayList<ArrayList<ArrayList<PVector>>> pointCloud;

	Chromosome chromosome;
	
	float splitRatio;

	//input parameters
	//organism
	float initSubs;    //initial subdivisions of first organism
	float initRad;     //initial radius of first organism
	
	//spring
	float minRestLen;  //minimum rest length for springs
	float springCoef;  //spring coefficient
	float maxLen;      //maximum spring length (what will it do? create a new cell?)
	float orgRepRadius = 15;
	
	int splits;
	public int overlaps = 0;

	Environment e;
	
	double realSize;


	  //aux variables
	int first_cell = -1;

	//=========CONSTRUCTORS================
	public Colony() {
		/*
		Generates empty colony
		*/
		this.organisms = new ArrayList<Organism>();
		this.splits = 0;
		
	}

	public Colony(Chromosome chr, Environment _e) {
		/*
	    Generates colony with one organism in environment _e;
	    Organism is initalised using values in chromosome chr.
	    */
	    this.chromosome = chr;
	    
	    this.splitRatio = 0.2f + ((float) this.chromosome.get(4) * 0.3f);

	    this.e = _e;

	    this.organisms = new ArrayList<Organism>();
	    // this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    Organism firstOrg = new Organism(this.chromosome, new PVector(this.e.width/2, this.e.height/2), 100);
	    firstOrg.setIndex(0);
	    this.organisms.add(firstOrg);
	    //set cells ids
//	    int cell_count = 0;
	    for(Organism o : this.organisms) {
		    for(int i = 0; i < o.cells.size(); i++) {
		    	Cell c = o.cells.get(i);
//		    	cell_count++;
//		    	c.setID(cell_count);
		    }
	    }
	    
	    this.splits = 0;
	}
	
	public Colony(Chromosome chr, Environment _e, ArrayList<PVector> _locations) {
		/*
	    Generates colony with organisms located at every point in _locations in environment _e;
	    Organisms are initalised using values in chromosome chr.
	    */
	    this.chromosome = chr;
	    
	    this.splitRatio = 0.2f + ((float) this.chromosome.get(4) * 0.3f);

	    this.e = _e;

	    this.organisms = new ArrayList<Organism>();
	    
	    for(PVector loc : _locations) {
	    // this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    	Organism org = new Organism(this.chromosome, loc, 50);
	    	org.setIndex(_locations.indexOf(loc));
	    	this.organisms.add(org);
	    }
	    
	    //set cells ids
	    int cell_count = 0;
	    for(Organism o : this.organisms) {
		    for(int i = 0; i < o.cells.size(); i++) {
		    	Cell c = o.cells.get(i);
		    	cell_count++;
		    	c.setID(cell_count);
		    }
	    }
	    
	    this.splits = 0;
	}

	public Colony(Environment _e, ArrayList<Organism> _orgs) {
	    /*
	    Generates colony in specific environment.
	    */
	    this.organisms = _orgs;
//	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    this.e = _e;
	    
	    this.splits = 0;
	}

	public Colony(Organism o) {
	    this.organisms = new ArrayList<Organism>();
	    this.organisms.add(o);
	    this.chromosome = o.getChromosome();
//	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    
	    this.splits = 0;
	    
	}

	public Colony(ArrayList<Organism> orgs) {
	    this.organisms = orgs;
//	    this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector>>>();
	    this.splits = 0;
	  }
	
	// COPY CONSTRUCTOR
	public Colony(Colony _c) {
		this.e = new Environment(_c.e);
		// Copy chromosome
		ArrayList<Double> genes = new ArrayList<Double>(_c.chromosome.getGenes());
		this.chromosome = new Chromosome(genes);
				
		this.organisms = new ArrayList<Organism>();
		
		//copy organisms
		for (Organism o : _c.organisms) {
			Organism newO = new Organism(o);
			this.organisms.add(newO);
		}
		this.splits = _c.getSplits();
	}
		
	// METHODS
	public int getSplits() {
		return this.splits;
	}
	
	public int getOverlaps() {
		return this.overlaps;
	}
	
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
	
	public void setRealSize(double _rs) {
		this.realSize = _rs;
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
	
	public float getHeight() {
		if(this.organisms.size() == 0) return 0f;
		float h = this.organisms.get(0).cells.get(0).loc.z;
		return h;
	}
	
	/*
	 * Update method
	 */
	public void update() {

		//check if ant organism is smaller than 3 cells and kill it
		this.killOrganisms();

	    this.checkFullOverlaps();

	    this.joinOverlapping(); //This function is to be implemented in the future
	    
	    this.applyGravity();
	    //update organisms
	    this.updateOrganisms();
	    
	} //update colony
	
	/*
	 * Update method with evaluation
	 */
	public void updateAndEval() {
		
//		this.applyOrgRepulsion();

		//check if ant organism is smaller than 3 cells and kill it
		this.killOrganisms();

	    this.checkFullOverlaps();

	    this.joinOverlapping(); //This function is to be implemented in the future
	    
	    //update organisms
	    this.updateAndEvalOrganisms();
	    
	} //update colony


	/*
	 * Update organisms
	 */
	void updateOrganisms() {
		for (int i = 0; i< this.organisms.size(); i++) {
			Organism o = this.organisms.get(i);
//			o.setTotalCells(this.getCellCount());
			if (o.splitting == "false") {
				o.update(e);
			} else if (o.splitting == "energy") {
				this.splitOrganism(i);
			} else if (o.splitting == "self") {
				this.splitSelfIntersect(i);
			}
		}
	} //updateOrganisms
	
//	void updateOrganisms() {
//		for (int i = 0; i< this.organisms.size(); i++) {
//			Organism o = this.organisms.get(i);
//			o.update(e);
//		}
//	}
	
	/*
	 * Update organisms
	 */
	void updateAndEvalOrganisms() {
		for (int i = 0; i< this.organisms.size(); i++) {
			Organism o = this.organisms.get(i);
			if (o.splitting == "false") {
				o.updateAndEval(e);
			} else if (o.splitting == "energy") {
				this.splitOrganism(i);
			} else if (o.splitting == "self") {
				this.splitSelfIntersect(i);
			}
		}
	} //updateAndEvalOrganisms

	void killOrganisms() {
		for (int i = 0; i < this.organisms.size(); i++) {
			Organism o = this.organisms.get(i);
			if (o.cells.size() < 3) {
//				System.out.println("Organism "+i+" was killed by cells");
				this.organisms.remove(i);
//				continue;
			}
			if (o.getTotalEnergy() <= 0.5) {
//				System.out.println("Organism "+i+" was killed by energy ("+o.getTotalEnergy()+")");
				this.organisms.remove(o);
			}
		}
	}


	void splitOrganism(int org) {
		  		  
		Organism o = this.organisms.get(org);
		//o.splitting = "energy";
		  
		//check self intersection
		o.selfIntersection();
		if(o.splitting == "self") return;
		  
		//get first point of division
		int first_cell = o.getHottestCell();

		//find remaining points
		int next_cell = first_cell;
		int last_cell;

		float split_energy = o.cells.get(first_cell).energy;
		int dist_counter = 0;
		for (int i = 0; i < o.cells.size(); i++) {
			if (split_energy < (o.getTotalEnergy() * this.splitRatio)) {
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
		force_to_first.mult(100f);

		PVector force_to_last = PVector.sub(o.cells.get(last_cell).loc, o.cells.get(first_cell).loc);
		force_to_last.mult(100f);
		
		// measure distance between split points
		float split_pts_dist = o.cells.get(first_cell).loc.dist(o.cells.get(last_cell).loc);
		  
		  // count cells to keep
		int cells_to_keep;
		if (first_cell < last_cell) {
			cells_to_keep = last_cell - first_cell + 1;
		} else {
			cells_to_keep = last_cell + (o.cells.size() - first_cell) + 1;
		}
  
		// determine number of cells to remove
		int cells_to_remove = o.cells.size() - cells_to_keep;
		int next_to_remove;
	  
		// split if the distance between split points is equal to or smaller than 2*repulsionRadius
		if (split_pts_dist <= 2*o.cells.get(first_cell).getRepRadius()) {
			int new_first = first_cell;
			  
			  // Initialise new organism
			Organism new_o = new Organism(this.chromosome);
			new_o.splitting = "energy";
			  			  
			for (int i = 0; i < cells_to_remove; i++) {
				if (new_first == 0) {
					next_to_remove = o.cells.size() - 1;
				} else {
					next_to_remove = new_first - 1;
					new_first--;
				}
				Cell new_c = new Cell(o.cells.get(next_to_remove));
				new_o.cells.add(new_c);
				o.killCell(next_to_remove);
			}
			o.splitting = "false";


			Collections.reverse(new_o.cells);
			new_o.connectCells();
//			new_o.setSpringCoef(0.5f);
			new_o.splitting = "false";
//			o.setSpringCoef(0.5f);
			o.resetInitEnergy();
			new_o.resetInitEnergy();
			this.organisms.add(new_o);
			this.splits++;
		} else {
			//get line between start and end
			float x1 = o.cells.get(first_cell).loc.x;
			float y1 = o.cells.get(first_cell).loc.y;
			float x2 = o.cells.get(last_cell).loc.x;
			float y2 = o.cells.get(last_cell).loc.y;

			PVector normal1 = VectorOps.getNormal(x1, y1, x2, y2);
			PVector normal2 = VectorOps.getNormal(x2, y2, x1, y1);

			normal1.mult(15f);
			normal2.mult(15f);

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

			o.cells.get(first_cell).applyForce2(force_to_last);
			//apply force_to_last to first_cell's neighbours
			//clockwise
			for (int i = 1; i <= cells_to_keep/2; i++) {
				int cwNeighbour;
				if (first_cell + i > o.cells.size() - 1) {
					cwNeighbour =  (first_cell + i) - o.cells.size();
				} else {
					cwNeighbour = first_cell + i;
				}
				PVector f = new PVector(force_to_last.x, force_to_last.y);
				f.mult((cells_to_keep/2-i)/(cells_to_keep/2));
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
				PVector f = new PVector(force_to_last.x,force_to_last.y);
				f.mult((cells_to_remove/2-i)/(cells_to_remove/2));
				o.cells.get(ccwNeighbour).applyForce(f);
			}

			o.cells.get(last_cell).applyForce2(force_to_first);
			//apply force_to_last to first_cell's neighbours
			//clockwise
			for (int i = 1; i <= cells_to_remove/2; i++) {
				int cwNeighbour;
				if (last_cell + i > o.cells.size() - 1) {
					cwNeighbour =  (last_cell + i) - o.cells.size();
				} else {
					cwNeighbour = last_cell + i;
				}
				PVector f = new PVector(force_to_first.x,force_to_first.y);
				f.mult((cells_to_remove/2-i)/(cells_to_remove/2));
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
				  PVector f = new PVector(force_to_first.x,force_to_first.y);
				  f.mult((cells_to_keep/2-i)/(cells_to_keep/2));
				  o.cells.get(ccwNeighbour).applyForce(f);
			  }

//			  o.applyRepulsion();
//			  o.applyAttraction();
		  }

		  for (Cell c : o.cells) {
			  c.update2();
		  }
		  for (Spring s : o.springs) {
			  s.update();
		  }
	  } //splitOrganism

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
		  
		  this.splits++;
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
	  
	  
	  // Join overlapping
	  public void joinOverlapping() {
		  if (this.organisms.size() < 2) return;
		  int cellCount = 0;
//		  Organism newOrg;
		  for(int i = 0; i<this.organisms.size()-1;i++) {
			  Organism org = this.organisms.get(i);
			  for(int j = i+1; j < this.organisms.size(); j++) {
				  Organism other = this.organisms.get(j);
				  ArrayList<Integer> cells_inside = new ArrayList<Integer>();
				  ArrayList<Integer> cells_inside_other = new ArrayList<Integer>();
				  
				  // check all the cells from org that are inside other
				  for(Cell cell : org.getCells()) {
					  if (cell.isInside(other)) {
						  int index = org.getCells().indexOf(cell);
						  cells_inside.add(index);
						  cellCount++;
					  }					
				  }
				  
				  // check all the cells from other that are inside of org
				  for(Cell cell : other.getCells()) {
					  if (cell.isInside(org)) {
						  int index = other.getCells().indexOf(cell);
						  cells_inside_other.add(index);
						  cellCount++;
					  }					
				  }
				  
				  if(cellCount > 0) {
					  Organism newOrg = new Organism();
					  if(cells_inside.size() > 0) {
						  if(cells_inside.get(0) != 0 || cells_inside.get(cells_inside.size()-1) != org.getCells().size()-1) {
							  org.shiftCells(cells_inside.get(cells_inside.size()-1) + 1);
							  int last = org.getCells().size()-1;

						  }else {
							  int shift_point = 0;
							  for(int s=0; s<cells_inside.size()-1;s++) {
								  if(cells_inside.get(s+1) - cells_inside.get(s) != 1) shift_point = cells_inside.get(s) + 1;
							  }
							  org.shiftCells(shift_point);
						  }
							
						  int last = org.getCells().size()-1;
						  for(int r=0;r<cells_inside.size();r++) {
							  org.getCells().remove(last-r);
						  }
//						  for(int r = cells_inside.size(); r < org.getCells().size();r++) {
//							  newOrg.addCell(new Cell(org.getCells().get(r)));
//						  }
							
						
					  }
					  for(Cell cell : org.getCells()) {
						  newOrg.addCell(new Cell(cell));
					  }
//					  System.out.println("New Org: "+newOrg.getTotalEnergy());
//					  System.out.println();
					  
					  if(cells_inside_other.size() > 0) {
						  
						  if(cells_inside_other.get(0) != 0 || cells_inside_other.get(cells_inside_other.size()-1) != other.getCells().size()-1) {
							  other.shiftCells(cells_inside_other.get(cells_inside_other.size()-1) + 1);
							  int last = other.getCells().size()-1;
//								for(int r=0;r<cells_inside.size();r++) {
//									org.getCells().remove(last-r);
//								}
						  }else {
							  int shift_point = 0;
							  for(int s=0; s<cells_inside_other.size()-1;s++) {
								  if(cells_inside_other.get(s+1) - cells_inside_other.get(s) != 1) shift_point = cells_inside_other.get(s) + 1;
							  }
							  other.shiftCells(shift_point);
						  }
							
						  int last = other.getCells().size()-1;
						  for(int r=0;r<cells_inside_other.size();r++) {
							  other.getCells().remove(last-r);
						  }
							
							
					  }
					  for(Cell cell : other.getCells()) {
						  newOrg.addCell(new Cell(cell));
					  }
					  newOrg.setChromosome(this.chromosome);
					  
//					  System.out.println("Org: "+org.getTotalEnergy());
//					  System.out.println("Other: "+other.getTotalEnergy());
//					  System.out.println("New Org: "+newOrg.getTotalEnergy());
//					  System.out.println();
					  
					  int orgIndex = this.organisms.indexOf(org);
					  
//					  System.out.println("Org "+orgIndex+" will be removed");
					  this.organisms.remove(orgIndex);
//					  System.out.println("Org "+orgIndex+" has been removed");
					  
					  
					  int otherIndex = this.organisms.indexOf(other);
//					  System.out.println("Org "+otherIndex+" will be removed");
					  
						
					  
					  this.organisms.remove(otherIndex);
//					  System.out.println("Org "+otherIndex+" has been removed");
//					  this.organisms.remove(org);
//					  this.organisms.remove(other);
					  this.organisms.add(newOrg);
//					  System.out.println("New Org has been added");
					  this.overlaps++;
					  return;
			
					  
				  }
			  }
		  }
	  }//joinOverlapping
	  
	  // gravity
	  public void applyGravity() {
		  if(this.organisms.size() == 1) return;
		  
		  //find largest organism by energy
		  int largestOrg = -1;
		  float largestOrgWeight = -1;
		  for(int i = 0; i<this.organisms.size() - 1; i++) {
			  float w = this.organisms.get(i).getTotalEnergy();
			  if(w > largestOrgWeight) {
				  largestOrg = i;
				  largestOrgWeight = w;
			  }
		  }
		  
		  if(largestOrg == -1) return;
		  
		  // apply force to all cells within OrgDiameter
		  // get heaviest org radius
		  float orgDiameter = this.organisms.get(largestOrg).maxDiameter();
		  
		// get heaviest org centroid
		  PVector orgCentroid = this.organisms.get(largestOrg).getHullCentroid();
		  
		  for(int i = 0; i < this.organisms.size(); i++) {
			  if (i != largestOrg) {
				  Organism o = this.organisms.get(i);
				  for(Cell c : o.getCells()) {
					  float dist = orgCentroid.dist(c.loc);
					  if(dist < orgDiameter) {
						  PVector gravity = PVector.sub(orgCentroid, c.loc).normalize();
						  gravity.mult((largestOrgWeight/dist)/3);
						  c.applyForce(gravity);
					  }
				  }
			  }
		  }
		  
	  }//gravity
	  
	  public double getArea() {
		  double a = 0;
		  for(Organism o : this.organisms) {
			  a+=(double)o.getArea();
		  }
		  return a;
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
	  
	  public void setChromosome(Chromosome _chr) {
		  this.chromosome = _chr;
	  }
	  
	  public Chromosome getChromosome() {
		  return this.chromosome;
	  }
	  
//	  public Cell getCellByID(int _id) {
//		  for (Organism o: this.organisms) {
//			  for(Cell c : o.getCells()) {
//				  if(c.getID() == _id) return c;
//			  }
//		  }
//		  
//	  }
		  
}
