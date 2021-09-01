package generative.components.examples;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
	Spring spr;
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
		
		e = new Environment();
		e.setRandomSeed(CommonConsts.getEnvSeed(25));
		e.setFoodSources(5);
		
		println(e.getParams());
		
		ArrayList<PVector> locations = new ArrayList<PVector>();
		locations.add(new PVector(width * .33f, height * .33f));
		locations.add(new PVector(width * .6f, height * .32f));
		locations.add(new PVector(width * .33f, height * .66f));
		locations.add(new PVector(width * .66f, height * .67f));
//		
		col = new Colony(chr, e, locations);	
		
//		frameRate(2);
		
	}
	
	public void draw() {
//		background(220);
//		text(col.organisms.size(),10,20);
		stroke(0,50);
		
		
		
		CommonFunctions.displayColony(this, col, 0);
//		CommonFunctions.displayEnvironment(this, e);
		
		col.update();
		e.updateNutrientField();
		
		
		for(Organism o : col.organisms) {
			PVector centroid = o.getHullCentroid();
//			text(col.organisms.indexOf(o), centroid.x, centroid.y);
			for(Cell cell : o.getCells()) {
//				PVector ft = cell.findFoodTarget(e);
//				if(ft != null) {
//					PVector appliedForce = PVector.sub(cell.loc, ft).normalize();
//					appliedForce.mult(cell.energy * (float)chr.getGenes2()[0]);	
////					println(appliedForce+" "+appliedForce.mag());
//					line(cell.loc.x,cell.loc.y,cell.loc.x+appliedForce.x, cell.loc.y+appliedForce.y);
//				}
//				if(cell.food) {
//					fill(0,255,0);
//					
////					
//					
//				}else {
//					fill(255,0,0);
//				}
//				CommonFunctions.drawVector(this, cell.vel, cell.loc.x, cell.loc.y, 20);
//				ellipse(cell.loc.x, cell.loc.y, cell.energy*2, cell.energy*2);
				
			}
		}
		for(FoodSource fs : e.getFoodSources()) {
			int col = fs.getCol();
			int row = fs.getRow();
			
//			println("Col: "+col+", row: "+row+", energy: "+e.nutrientField[col][row]);
		}
		
		int[] tile = e.getCurrentTile(new PVector(mouseX, mouseY));
		float gridSize = width/e.getCols();
		float hPos = tile[0]*gridSize;
		float vPos = tile[1]*gridSize;
		
//		pushMatrix();
//		translate(hPos, vPos);
//		noFill();
//		strokeWeight(3);
//		rectMode(1);
//		rect(0,0,gridSize,gridSize);
//		strokeWeight(1);
//		popMatrix();

	}
	
	public void mouseReleased() {
//		println(mouseX, mouseY);
//		PVector l = new PVector(mouseX, mouseY);
//		col.organisms.add(new Organism(chr, l, 50));
//		if(col.organisms.size() > 1) {
//			noLoop();
//		}
//		int[] tile = e.getCurrentTile(new PVector(mouseX, mouseY));
//		FoodSource fs = new FoodSource(tile[0], tile[1], 60);
//		e.addFoodSource(fs);
//		e.updateNutrientField();
//		println("Col: "+col+", row: "+row+e.getFoodSources().length);
		
//		col.update();
//		e.updateNutrientField();
			
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
