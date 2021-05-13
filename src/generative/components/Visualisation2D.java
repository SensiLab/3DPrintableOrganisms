package generative.components;

import generative.components.CommonConsts;
import processing.core.PApplet;
import processing.core.PVector;

public class Visualisation2D extends PApplet {
	
	long envSeed = CommonConsts.getEnvSeed(25);
	
	double[] genes = {0.7971354411267045,0.679551696552745,0.43895084203454265,0.7669021860943841,0.5715081078078987};
	Chromosome chr = new Chromosome(genes);
	// Environment size set in millimeters
	Environment e;
	Colony c;
//	Simulation s = new Simulation(chr,e,100,0);
//	Colony3D col;
	int counter = 0;
	
	public void settings() {
		size(600,600, P2D);
		e = new Environment(600, 600, 15, 0.4f, 5, 2, 5, 5, 0.2f, envSeed);
		c = new Colony(chr, e);
		
		PVector v1 = new PVector(30,60);
		v1.div(3);
		println(v1);
//		s.generate();
//		col = s.colony3d;
			
		}
	
//	@Override
//	public void setup() {
//		col.setPApplet(this);
//		background(220);
//		
//	}
	
//	public void draw() {
//		if(counter > 520) {
//			fill(0);
//			text("Counter:\t"+(counter-1),10,15);
//			text("Env Size:\t"+e.getWidth(),10,30); 
//			noLoop();
//		}
////		background(220);
//        stroke(20,50);
//        strokeWeight(1f);
//        displayColony(c);
////        displayEnvironment(e);
//        e.updateNutrientField(this.c);
////        for (int i = 0; i < this.colonyLayers + this.colonyWarmupLayers; i++) {
//			
//		      
//		this.c.update();
//		counter++;
//
//	}
//	
//	public void displayColony(Colony _c){
//		// Scale colony based on drawing size
//		
//	    for(Organism o : _c.organisms){
//	        for(Spring s : o.getSprings()){
//	        	float spX = s.sp.loc.x * width/_c.e.width;
//	        	float spY = s.sp.loc.y * height/_c.e.height;
//	        	float epX = s.ep.loc.x * width/_c.e.width;
//	        	float epY = s.ep.loc.y * height/_c.e.height;
//	            line(spX,spY,epX,epY);
////	            ellipse(spX, spY, width * .01f, height * .01f);
//	        }
//	        
//	    }
//	}
	
//	public void displayColony(Colony _c){
//		// Scale colony based on drawing size
//		
//	    for(Organism o : _c.organisms){
//	        for(Spring s : o.getSprings()){
//	        	float spX = s.sp.loc.x * width/_c.e.width;
//	        	float spY = s.sp.loc.y * height/_c.e.height;
//	        	float epX = s.ep.loc.x * width/_c.e.width;
//	        	float epY = s.ep.loc.y * height/_c.e.height;
//	            line(spX,spY,epX,epY);
////	            ellipse(spX, spY, width * .01f, height * .01f);
//	        }
//	        
//	    }
//	}
	
	public void displayEnvironment(Environment e) {
		float gridSize = width/e.cols;
		//draw grid
		stroke(100f);
		strokeWeight(0.5f);
		for(int i = 1; i < e.cols; i++) {
			float x = i * gridSize;
			line(x,0,x,height);
			line(0,x,width,x);
		}
		//draw food
		for(int i = 0; i < e.cols; i++) {
			for(int j = 0; j < e.rows; j++) {
				int[] tile = {i,j};
				if(e.lookupNutrients(tile) > 0) {
					float r = e.lookupNutrients(tile) * gridSize;
					float x = tile[0] * gridSize + gridSize/2;
					float y = tile[1] * gridSize + gridSize/2;
					ellipse(x,y,r,r);
				}
			}
		}
	}
	
	public static void main(String[] args) {
		String[] processingArgs = {"Visualisation2D"};
		Visualisation2D vis2D = new Visualisation2D();
		PApplet.runSketch(processingArgs, vis2D);
		
	}
	
	

}
