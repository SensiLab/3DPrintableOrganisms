package generative.components.examples;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import processing.core.*;
import generative.components.*;
import peasy.*;

public class VisTests3D extends PApplet{
	
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
		size(600,600,P3D);
//		size(600,600,P2D);
		
	}
	
	public void setup() {
		camera = new PeasyCam(this, width/2, height/2, 0, 1.5*height);
	    camera.rotateZ(PI/4);
	    camera.rotateX(-PI/4);
	    camera.pan(0,-width/2.35);
		
		generator = new Random();
		genes = new double[5];
		for(int i=0; i<5; i++) {
			genes[i] = generator.nextDouble();
//			genes[i] = 0.5;
		}
		chr = new Chromosome(genes);
		
		e = new Environment();
		e.setRandomSeed(CommonConsts.getEnvSeed(40));
		e.setFoodSources(10);
		
		println(e.getParams());
		
		ArrayList<PVector> locations = new ArrayList<PVector>();
		locations.add(new PVector(width * .33f, height * .33f));
		locations.add(new PVector(width * .6f, height * .32f));
		locations.add(new PVector(width * .33f, height * .66f));
		locations.add(new PVector(width * .66f, height * .67f));
//		
//		col = new Colony(chr, e, locations);	
		s = new Simulation(chr,e,200,200,locations);
		s.generate();
		
		println("Total overlaps: "+s.getColony().overlaps);
		
//		frameRate(2);
		
	}
	
	public void draw() {
		background(220);
//		text(col.organisms.size(),10,20);
		stroke(0,100);
		
		CommonFunctions.displayColony3D(this, s.getColony3D(), 3);
		
		
		
		
	}
	
		
	
	public static void main(String[] args) {
		String[] processingArgs = {"VisTests"};
		VisTests3D vis2D = new VisTests3D();
		PApplet.runSketch(processingArgs, vis2D);
		
	}
}
