package generative.components.examples;

import java.util.ArrayList;
import java.util.Random;
import processing.core.*;
import generative.components.*;
import peasy.*;

public class VisTests extends PApplet{
	
	PeasyCam camera;
	
	Random generator;
	
	double[] genes;
	Chromosome chr;
	// Environment size set in millimeters
	Environment e;
	Cell c;
	Cell c2;
	ArrayList<Cell> cells;
	Organism org;
	Colony col;
	Simulation s;
	Colony3D col3d;
	
	
	
	public void settings() {
//		size(600,600,P3D);
		size(600,600,P2D);
		
	}
	
	public void setup() {
//		camera = new PeasyCam(this, width/2, height/2, 0, 1.5*height);
//	    camera.rotateZ(PI/4);
//	    camera.rotateX(-PI/4);
//	    camera.pan(0,-width/2.35);
		
		generator = new Random();
		genes = new double[5];
		for(int i=0; i<5; i++) {
//			genes[i] = generator.nextDouble();
			genes[i] = 0.5;
		}
		chr = new Chromosome(genes);
//		c = new Cell(width*.33f, height/2, chr);
//		c2 = new Cell(width*.66f, height/2, chr);
		e = new Environment();
		e.setRandomSeed(CommonConsts.getEnvSeed(25));
		e.setFoodSources(5);
		
		ArrayList<PVector> locations = new ArrayList<PVector>();
		locations.add(new PVector(width * .25f, height * .5f));
		locations.add(new PVector(width * .75f, height * .5f));
//		locations.add(new PVector(width * .33f, height * .66f));
//		locations.add(new PVector(width * .66f, height * .66f));
//		
//		col = new Colony(chr, e, locations);
		col = new Colony();
		col.setEnvironment(e);
		col.setChromosome(chr);
		
//		cells = new ArrayList<Cell>();
//		for(int i = 0; i < 10; i++) {
//			Cell newCell = new Cell(random(width), random(height), chr);
//			cells.add(newCell);
//		}
////		cells.add(c);
//		s = new Simulation(chr, e, 500, 20);
//		s.generate();
		
////		println(s.getColony3D().getLayer(0).organisms.get(0).getCells().get(0).getID());
//		for(Colony _col : s.getColony3D().getLayers()) {
//			for(Organism _o : _col.organisms) {
//				println(s.getColony3D().getLayers().indexOf(_col)+ ": Total Cells: "+_o.getCells().size());
////				for(Cell _c : _o.getCells()) {
////					if(_c.getSplitted()) {
////						println("Cell "+_o.getCells().indexOf(_c)+" on layer"+s.getColony3D().getLayers().indexOf(_col)+" has splitted");
////					}
////				}
//			}
//		}
		
//		org = new Organism(chr, new PVector(width/2,height/2), 20, 50f);
//		org.addCell(c);
//		org.addCell(c2);
//		org.springs.add(new Spring(c,c2));
		
	}
	
	public void draw() {
		background(220);
		text(col.organisms.size(),10,20);
		stroke(0,100);
		
//		c.update(e);
		col.update();
		e.updateNutrientField(col);
		
//		for(Organism o : col.organisms) {
//			if(isInside(o)) println("True");
//		}
		
		CommonFunctions.displayColony(this, col, 0);

		
		
//		int cellCount = 0;
//		for(int i = 0; i < col.organisms.size() - 1; i++) {
//			Organism org = col.organisms.get(i);
//			for(int j = i+1; j < col.organisms.size(); j++) {
//				Organism other = col.organisms.get(j);
//				
//				Organism orgCopy = new Organism();
//				for(Cell cell : org.getCells()) {
//					if (cell.isInside(other)) {
//						stroke(255,0,0);
//						ellipse(cell.loc.x, cell.loc.y, cell.energy, cell.energy);
//						cellCount++;
//					}					
//					orgCopy.addCell(cell);
//				}
//				
//				Organism otherCopy = new Organism();
//				for(Cell cell : other.getCells()) {
//					if (cell.isInside(org)) {
//						stroke(0,0,255);
//						ellipse(cell.loc.x, cell.loc.y, cell.energy, cell.energy);
//						cellCount++;
//					}					
//					otherCopy.addCell(cell);
//				}
//				
//				
//				if(cellCount > 0) {
//					// initialise new organism
//					Organism newOrg = new Organism();
//					
//					// find new first for org
//					int newFirst;
//					if(orgCopy.getCells().get(0).isInside) {
//						for(int k = 0; i < orgCopy.getCells().size(); k++) {
//							if(!orgCopy.getCells().get(k).isInside) {
//								newFirst = k;
//								break;
//							}
//						}
//					}else{
//						
//					}
//					int first_inside_org = -1;
//					int cells_inside_org = 0;
//					for(int k = 0; k < orgCopy.getCells().size(); k++) {
//						Cell cell = orgCopy.getCells().get(k);
//						if(cell.isInside) {
//							if(first_inside_org == -1) first_inside_org = k;
//							cells_inside_org++;
//						}
//					}
//					//shift cells in orgCopy
//					
//					
//					for(Cell cell : org.getCells()) {
//						if (!cell.isInside(other)) {
//							newOrg.addCell(cell);
//						}					
//					}
//					for(Cell cell : other.getCells()) {
//						if (!cell.isInside(org)) {
//							newOrg.addCell(cell);
//						}					
//					}
//					newOrg.polarSort();
//					newOrg.connectCells();
//					for(Spring s:newOrg.getSprings()) {
//						line(s.sp.loc.x, s.sp.loc.y, s.ep.loc.x, s.ep.loc.y);
//					}
//					noLoop();
////					col.organisms.add(newOrg);
//				}
//			}
//		}
//		text(cellCount,10,35);
		
//		ellipse(c.loc.x, c.loc.y, c.energy, c.energy);
		
//		CommonFunctions.displayColony3D(this, s.getColony3D(), 3);
//		CommonFunctions.displayColony3DMesh(this, s.getColony3D());
		
		
//		strokeWeight(8);
		
//		fill(255);
//		ellipseMode(DIAMETER);
//		CommonFunctions.displayEnvironment(this, e);
//		
//		ellipseMode(RADIUS);		
//		for(int i=0; i<org.getCells().size(); i++) {
//		for(int i =0; i < cells.size();i++) {
////			Cell cell = org.getCells().get(i);
//			Cell cell = cells.get(i);
//			if(cell.state == "seek") {
//				cell.seekFood(e);
//				cell.eat(e);				
//			}
////		c.addDrag(e);
//			cell.update(e);
////			println(cell.wander);
////			if(cell.energy >= cell.getMaxEnergy()) {
////				org.splitCell(cell);
////			}
//			fill(0,0,128);
//			ellipse(cell.loc.x, cell.loc.y, cell.energy, cell.energy);
//			
////			
//			noFill();
//			if(cell.state == "seek") {
//				stroke(255,0,0);
//			}else{
//				stroke(0,0,255);
//			}
//			ellipse(cell.loc.x, cell.loc.y, c.getMaxEnergy(), c.getMaxEnergy());
//			PVector v = new PVector(cell.vel.x, cell.vel.y);
//			PVector vel = PVector.add(cell.loc, v.mult(10));
//			line(cell.loc.x, cell.loc.y, vel.x, vel.y);
////			
//		}
//		for(int i = 0; i < org.getCells().size(); i++) {
//			Cell c = org.getCells().get(i);
//			Spring s = org.getSprings().get(i);
//			line(s.sp.loc.x,s.sp.loc.y, s.ep.loc.x, s.ep.loc.y);
////			ellipse(c.loc.x, c.loc.y,8,8);
//			
//		}
//		org.update(e);
		
		
//		for(Spring s : org.getSprings()) {
//			s.update();
//			stroke(0,120);
//			line(s.sp.loc.x, s.sp.loc.y, s.ep.loc.x, s.ep.loc.y);
//		}
//		org.applyRepulsion();
		
//		e.updateNutrientField(col);
		
		
		
//		println("Current Energy: "+c.energy+"\t Max Energy: "+c.getMaxEnergy());
		
	}
	
	public void mouseReleased() {
//		println(mouseX, mouseY);
		PVector l = new PVector(mouseX, mouseY);
		col.organisms.add(new Organism(chr, l, 40));
//		if(col.organisms.size() > 1) {
//			noLoop();
//		}
			
	}
	
	public Boolean isInside(Organism o) {
		  float x3 = mouseX;
		  float y3 = mouseY;

		  float x4 = 9999;	//this is a constant to draw a long horizontal line
		  float y4 = mouseY;

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
	
	
	public static void main(String[] args) {
		String[] processingArgs = {"VisTests"};
		VisTests vis2D = new VisTests();
		PApplet.runSketch(processingArgs, vis2D);
		
	}
}
