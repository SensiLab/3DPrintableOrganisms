package generative.components;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
//import java.util.Arrays;
import java.util.Collections;
//import java.util.Comparator;

//import fr.inria.optimization.cmaes.IntDouble;
import processing.core.PVector;
import processing.core.PApplet;


public class Simulation {
	/*
	Simulation takes care of developing a colony in time and evaluating its fitness.
	Inputs: A chromosome and an environment.
	Process:
	0. Initialise environment using resolution, number of food sources, food
	source size range (min food, max food) and food growth rate.
	1. Initialise a colony on environment using a chromosome
	2. Generate colony point cloud
	3. evaluate fitness of colony


	*/
	//Environment init parameters
	Environment e;
	Chromosome chromosome;
	int colonyHeight;
	int colonyWarmupHeight;
	
	int colonyLayers;
	int colonyWarmupLayers;
	
	int envRes = 20;
	ArrayList<ArrayList<ArrayList<PVector[]>>> pointCloud;
	Colony3D colony3d;
	Colony c;
	float costOfLiving;
	float printabilityThreshold = 0.2f;
	float minPrintDiameter = 5f;
	float outputSize = 120f;
	float outputLayerHeight = 0.2f;
	
	float angleDispersion;
	float printability;

	// MAIN CONSTRUCTOR	
	public Simulation(Chromosome _c, Environment _e, int _colonyHeight, int _colonyWarmupHeight){
		this.e = _e;
		this.chromosome = _c;
		this.colonyHeight = _colonyHeight;
		this.colonyWarmupHeight = _colonyWarmupHeight;

		this.c = new Colony(this.chromosome, this.e);

		this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector[]>>>();
		this.colony3d = new Colony3D();
		
		this.init();

		// println("Simulation initialised\n"+" Env: "+this.e+"\n");

	} //constructor
	
	public Simulation(Chromosome _c, Environment _e, int _colonyHeight, int _colonyWarmupHeight, ArrayList<PVector> _colonyLocs){
		this.e = _e;
		this.chromosome = _c;
		this.colonyHeight = _colonyHeight;
		this.colonyWarmupHeight = _colonyWarmupHeight;

		this.c = new Colony(this.chromosome, this.e, _colonyLocs);

		this.pointCloud = new ArrayList<ArrayList<ArrayList<PVector[]>>>();
		this.colony3d = new Colony3D();

		// println("Simulation initialised\n"+" Env: "+this.e+"\n");
		this.init();

	} //constructor
		
	/** INIT
	 * Method to set the number of layers based on output size
	 */
	void init() {
		float realHeight = (float)this.colonyHeight * this.outputSize/e.width;
		float warmupRealHeight = (float)this.colonyWarmupHeight * this.outputSize/e.width;
		
		
		this.colonyLayers = (int)(realHeight/this.outputLayerHeight);
		this.colonyWarmupLayers = (int)(warmupRealHeight/this.outputLayerHeight);
	} //init


	//======================GENERATE==================
	// Method to generate the geometric data of an object
	public void generate() {
		
		//copy environment
		Environment env = new Environment(this.e);
		
		// set environment for colony
		this.c.setEnvironment(env);
		
		//		this.pointCloud.clear();
		this.colony3d.clear();

		for (int i = 0; i < this.colonyLayers + this.colonyWarmupLayers; i++) {
			env.updateNutrientField(this.c);
		      
			this.c.update();

			if (i >= this.colonyWarmupLayers) {
				
				Colony layer;
				// Copy colony
				if (this.c.organisms.size() != 0) {
					layer = new Colony(this.c);
					layer.setZPosition((float) (i - this.colonyWarmupLayers));
				} else {
					layer = new Colony();
				}
				layer.setEnvironment(this.e);
				// Add layer to colony 3d
				this.colony3d.addLayer(layer);

			}
		}
	} // generate
	
	//=====================UTILITY METHODS===============
	
	// return the size of the environment in millimiters
	public float getOutputSize() {
		return this.outputSize;
	}
	
	// set size of environment in millimiters
	public void setOutputSize(float _os) {
		this.outputSize = _os;
		this.init();
	}
	
	// return the number of layers of a colony scaled to be printed
	public int getColonyLayers() {
		return this.colonyLayers;
	}
	
	// returns the colony3D object currently in memory
	public Colony3D getColony3D() {
		return this.colony3d;
	}
	
	public Environment getEnvironment() {
		return this.e;
	}
	
	public Colony getColony() {
		return this.c;
	}
	
	//====================OUTPUT METHODS==================
	
	// generate gcode to 3D print a colony	
	public void printColony(boolean _timelapse, String fileName) {
		if(this.colony3d.layers.size() < 1) return;
//		Colony3D scaledColony3d = this.colony3d.scale(this.outputSize, this.outputLayerHeight, this.e);
		
		// initialise gCode
		GCode gCode = new GCode();
		
		gCode.setTimeLapse(_timelapse);
		
		gCode.setOutputSize(this.outputSize);
		
		// iterate through layers
		for(Colony layer : this.colony3d.layers) {
			int layerNumber = this.colony3d.layers.indexOf(layer);
//			if(layerNumber < 3) gCode.setMultiplier(1.25f);
//			if(layerNumber >= 3) gCode.setMultiplier(1.05f);
			gCode.writeGcodeLayer(layerNumber, layer);
			System.out.println("Layer "+layerNumber+" added to gCode");
		}
		
		
		Colony lastLayer = this.colony3d.layers.get(this.colony3d.layers.size() - 1);
		Organism lastOrg = lastLayer.organisms.get(lastLayer.organisms.size() -1);
		PVector finalPoint = lastOrg.cells.get(lastOrg.cells.size() - 1).loc;
		
		gCode.export(fileName, finalPoint);
		
	}//printColony
	
	// Save colony as .xyz pointcloud
	public void saveColony3D(String fileName) {
		int res = 3;	// resolution of the mesh
		System.out.println("saveColony3D called");
		// Create file
		File output = new File(fileName+".xyz");
		
		// Create string
		String pc = new String();
		
		// Add points to string
		// for every layer
		for(Colony col : this.colony3d.layers) {
			int layerNumber = this.colony3d.layers.indexOf(col);
//			System.out.println("Colony "+layerNumber+" being processed");
			// for every organism
			for(Organism o : col.organisms) {
				int orgIndex = col.organisms.indexOf(o);
//				System.out.println("Organism "+layerNumber+"-"+orgIndex+" being processed");
				// for every cell
//				for(Cell c : o.getCells()) {
				for(int i=0; i < o.cells.size(); i++) {
//					System.out.println("Cell "+i+"/"+o.cells.size()+" being processed");
					Cell c = o.cells.get(i);
					Cell next;
					if(i == o.cells.size() - 1) {
						next = o.cells.get(0);
					}else {
						next = o.cells.get(i+1);
					}
					// get normals
					PVector pointer = PVector.sub(next.loc, c.loc);
					PVector normal = new PVector(pointer.x, pointer.y).normalize().rotate((float)-Math.PI/2);
					
					pc+=(c.loc.x+" "+c.loc.y+" "+c.loc.z+" "+"128"+" "+"128"+" "+"128"+" "+normal.x+" "+normal.y+" "+normal.z+"\n");
					
					// add intermediate points
					float len = c.loc.dist(next.loc);
					PVector increment = PVector.div(pointer, res);
					for(int j = 1; j < res; j++) {
						PVector additional = PVector.mult(increment, j);
						PVector np = PVector.add(c.loc,  additional);
						pc+=(np.x+" "+np.y+" "+np.z+" "+"128"+" "+"128"+" "+"128"+" "+normal.x+" "+normal.y+" "+normal.z+"\n");
					}
					
				}
			}
			System.out.println("Layer "+layerNumber+" added");
		}	
		
		// Write to file
		try {
			FileWriter outputWriter = new FileWriter(fileName+".xyz");
			outputWriter.write(pc);
			outputWriter.close();
		}catch(IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}// saveColony3D
	
	
	/**
	 * Utility methods for the evaluation of 3d colonies
	 */
	
	/**
	 * Diameter factor is the probability that an organism will print
	 * correctly based on the diameter of its convex hull.
	 * @param d: diameter of the convex hull of the organism
	 * @param threshold: constant that determines the minumum diameter for printing
	 * @return
	 */
	float diameterFactor(float d, float threshold){
		if(d >= threshold) {
			return 1;
		}else {
			return d / threshold;
		}		
	}
	
	/**
	 * Method to find the closest spring to a given cell 
	 * @param pt (location of the cell)
	 * @param layer (the colony in which the spring is)
	 * @return array with the index of the organism and the spring in the colony
	 */
	static int[] closestEdge(PVector pt, Colony layer){
		  float minDist = 99999999;
		  int closestOrg = -1;
		  int closestEdge = -1;
		  for(int i = 0; i < layer.organisms.size(); i++){
			  Organism organism = layer.organisms.get(i);
			  for(int j = 0; j < organism.springs.size(); j++){
				  PVector[] edge = {organism.springs.get(j).sp.loc, organism.springs.get(j).ep.loc};
				  float d = distPointLine(pt, edge);
				  if(d < minDist){
					  minDist = d;
					  closestOrg = i;
					  closestEdge = j;
				  }
			  }
		  }
		  int[] r = {closestOrg, closestEdge};
		  return r;
	}//closestEdge
	
	
	/**
	 * 
	 * @param d: the distance between a point and the closest edge on the layer below
	 * @param threshold: the tolerance for printability based on the dimensions of the printing filament
	 * @return a value between 0 (non printable) and 1 fully printable
	 */
	static float printabilityScore(float d, float threshold){
	    if(d < threshold) return 1;
	    if(d < (2 * threshold)) return 1 - ((d - threshold)/threshold);
	    return 0;
	  }
	
	/**
	 * Method to calculate if the 3d colony developed fully
	 */
	public float isComplete() {
		float expectedLayers = (float)this.colony3d.layers.size();
		float actualLayers = 0;
		for(Colony layer : this.colony3d.layers) {
			if(layer.organisms.size() > 0) actualLayers++;
		}
		return actualLayers/expectedLayers;
	}
	
	
	/**
	 * Method to calculate printability of a 3D colony
	 * Printability looks at the minimum diameter of the whole shape
	 * at every layer (this includes 'branches') and the support of each vertex
	 * @return a float between 0 and 1.
	 */
	public float printability() {
		float overallPerimeter = 0f;
		float supportedPerimeter = 0f;
		
//		 scale 3d colony
		Colony3D scaledColony3d = this.colony3d.scale(this.outputSize, this.outputLayerHeight, this.e);
		
		// Initialise array to store scores per vertex
		// The array is nested in order to keep track of scores of vertices in layers below
		ArrayList<ArrayList<ArrayList<Float>>> colonyScores = new ArrayList<ArrayList<ArrayList<Float>>>();
		
		// add scores for layer 0
		ArrayList<ArrayList<Float>> layer0Scores = new ArrayList<ArrayList<Float>>(); 
		// for every organism in layer 0
		for(Organism o : scaledColony3d.layers.get(0).organisms) {
			ArrayList<Float> orgLayer0Scores = new ArrayList<Float>();
			
			ArrayList<PVector> cHull = o.getConvexHull();
			float hullMinDiameter = minDiameter(cHull);
			float diameterFactor = diameterFactor(hullMinDiameter, this.minPrintDiameter);
			
			// for every spring in o
			for(Spring s : o.springs) {
				// calculate the score for every spring
				float edgeLen = s.getLen();
				float edgeScore = diameterFactor;
				
				// add these scores to the overall counters
				overallPerimeter+=edgeLen;
				supportedPerimeter+=(edgeScore*edgeLen);
				
				// add edgeScore to the record
				orgLayer0Scores.add(edgeScore);
			}
			layer0Scores.add(orgLayer0Scores);
		}
		colonyScores.add(layer0Scores);
		
		// add scores for all other layers
		for(int i = 1; i < scaledColony3d.layers.size(); i++) {
			ArrayList<ArrayList<Float>> layerScores = new ArrayList<ArrayList<Float>>();
			Colony layer = scaledColony3d.layers.get(i);
			Colony layerBelow = scaledColony3d.layers.get(i-1);
			
			//for every organism in layer
			for(int j = 0; j < layer.organisms.size(); j++) {
				Organism o = layer.organisms.get(j);
				
				// calculate the diameter factor here
				ArrayList<PVector> cHull = o.getConvexHull();
				float hullMinDiameter = minDiameter(cHull);
				float diameterFactor = diameterFactor(hullMinDiameter, this.minPrintDiameter);				
				
				ArrayList<Float> orgScores = new ArrayList<Float>();
				//for every spring in organism
				for(int k = 0; k < o.springs.size(); k++) {
					Spring s = o.springs.get(k);
					float edgeLen = s.getLen();

					// calculate the score for this edge here
					
					//check support for start point
					PVector sp = s.sp.loc;
					// Get the index of the organism and spring in the colony layerBelow
					int[] spClosestEdgeCoordinates = closestEdge(sp, layerBelow);
					// get the support score of the edge closest to sp on layer below
					float closestEdgeScoreSp = colonyScores.get(i-1).get(spClosestEdgeCoordinates[0]).get(spClosestEdgeCoordinates[1]);
					// Get the actual spring
					Spring closestSpringSp = layerBelow.organisms.get(spClosestEdgeCoordinates[0]).springs.get(spClosestEdgeCoordinates[1]);
					// make a reference to the positions of the closest spring's sp and ep
					PVector[] closestEdgeSp = {closestSpringSp.sp.loc, closestSpringSp.ep.loc};
					
					float distSP = distPointLine(sp, closestEdgeSp);
					float scoreSp = printabilityScore(distSP,this.printabilityThreshold) * closestEdgeScoreSp;
					
					//check support for mid point
					PVector mp = VectorOps.getMidPoint(s.sp.loc, s.ep.loc);
					int[] mpClosestEdgeCoordinates = closestEdge(mp, layerBelow);
					// get the support score of the edge closest to mp on layer below
					float closestEdgeScoreMp = colonyScores.get(i-1).get(mpClosestEdgeCoordinates[0]).get(mpClosestEdgeCoordinates[1]);
					Spring closestSpringMp = layerBelow.organisms.get(mpClosestEdgeCoordinates[0]).springs.get(mpClosestEdgeCoordinates[1]);
					PVector[] closestEdgeMp = {closestSpringMp.sp.loc, closestSpringMp.ep.loc};
					float distMP = distPointLine(mp, closestEdgeMp);
					float scoreMp = printabilityScore(distMP,this.printabilityThreshold) * closestEdgeScoreMp;
					
					//check support for end point
					PVector ep = s.ep.loc;
					int[] epClosestEdgeCoordinates = closestEdge(ep, layerBelow);
					// get the support score of the edge closest to ep on layer below
					float closestEdgeScoreEp = colonyScores.get(i-1).get(epClosestEdgeCoordinates[0]).get(epClosestEdgeCoordinates[1]);
					Spring closestSpringEp = layerBelow.organisms.get(epClosestEdgeCoordinates[0]).springs.get(epClosestEdgeCoordinates[1]);
					PVector[] closestEdgeEp = {closestSpringEp.sp.loc, closestSpringEp.ep.loc};
					float distEP = distPointLine(ep, closestEdgeEp);
					float scoreEp = printabilityScore(distEP,this.printabilityThreshold) * closestEdgeScoreEp;
					
					float edgeScore = (diameterFactor * ((scoreSp + (2 * scoreMp) + scoreEp )/4));
					
					// add scores to overall counters
					overallPerimeter+=edgeLen;
					supportedPerimeter+=(edgeScore*edgeLen);
					orgScores.add(edgeScore);
				}
				layerScores.add(orgScores);
			}
			colonyScores.add(layerScores);
		}
		float result = (supportedPerimeter/overallPerimeter);
		if(Float.isNaN(result)) {
			return 0f;
		}else {
			return result;
		}
	}
	
	
	/**
	 * Calculates the distance between a point and a line on the xy plane
	 * @param pt
	 * @param ln
	 * @return distance as a float
	 */
	
	static float distPointLine(PVector pt, PVector[] ln) {
		  float x = pt.x;
		  float y = pt.y;

		  float x1 = ln[0].x;
		  float y1 = ln[0].y;

		  float x2 = ln[1].x;
		  float y2 = ln[1].y;

		  float d = (float) (Math.abs(((y2 - y1)*x)-((x2-x1)*y) + (x2*y1) - (y2*x1))/Math.sqrt(Math.pow(y2 - y1, 2) + Math.pow(x2 - x1,2)));
		  return d;
		
	}
	
	/**
	 * Finds the farthest point to given line
	 * @param line as an PVector[] with a start point and an end point
	 * @param array of points
	 * @return the distance between the line and the point farthest away from it
	 */
	static float maxLinePointDist(PVector[] ln, ArrayList<PVector> points){
		float maxDist = 0;
		for(PVector pt : points){
			float d = distPointLine(pt, ln);
			if(d > maxDist) maxDist = d;
		}
		return maxDist;
	} 

	/**
	 * 
	 * @param organism
	 * @return
	 */
	public static float minDiameter(ArrayList<PVector> organism){
		  // println(" called");
		  ArrayList<PVector[]> edges = new ArrayList<PVector[]>();
		  for(int i = 0; i < organism.size() - 1; i++) {
			  PVector[] edge = {organism.get(i), organism.get(i+1)};
			  edges.add(edge);
		  }
		  PVector[] lastEdge = {organism.get(organism.size() - 1), organism.get(0)};
		  edges.add(lastEdge);
		  
		  ArrayList<Float> diameters = new ArrayList<Float>();
		  for(PVector[] l : edges){
		    float d = maxLinePointDist(l,organism);
		    diameters.add(d);
		    // println(d);
		  }
		  Collections.sort(diameters, null);
		  // println(diameters);
		  return diameters.get(0);

	}
	
	public static float meanDiameter(ArrayList<PVector> organism){
		  // println(" called");
		  ArrayList<PVector[]> edges = new ArrayList<PVector[]>();
		  for(int i = 0; i < organism.size() - 1; i++) {
			  PVector[] edge = {organism.get(i), organism.get(i+1)};
			  edges.add(edge);
		  }
		  PVector[] lastEdge = {organism.get(organism.size() - 1), organism.get(0)};
		  edges.add(lastEdge);
		  
		  ArrayList<Float> diameters = new ArrayList<Float>();
		  for(PVector[] l : edges){
		    float d = maxLinePointDist(l,organism);
		    diameters.add(d);
		    // println(d);
		  }
		  float totalDiameter = 0;
		  for(float d : diameters) {
			  totalDiameter+=d;
		  }
		  
		  return totalDiameter/(float)diameters.size();

	}
	
	public static PVector getIntersectionPoint(PVector[] a, PVector[] b) {
	    float x1 = b[0].x;
	    float y1 = b[0].y;
	    float x2 = b[1].x;
	    float y2 = b[1].y;

	    float x3 = a[0].x;
	    float y3 = a[0].y;
	    float x4 = a[1].x;
	    float y4 = a[1].y;

	    float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	    if (den == 0) return null;

	    float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
	    float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;
	    
	    if (t > 1f || t < 0f || u > 1f || u < 0f) return null;
	    
	    float x = x1 + (t*(x2-x1));
	    float y = y1 + (t*(y2-y1));

	    return new PVector(x,y);
	}
	
	public static Boolean isInside(PVector p, ArrayList<PVector> o) {
		ArrayList<PVector[]> edges = new ArrayList<PVector[]>();
		for(int i = 0; i < o.size(); i++) {
			int next;
			if(i == o.size()-1) {
				next = 0;
			}else {
				next = i+1;
			}
			PVector[] edge = {o.get(i), o.get(next)};
			edges.add(edge);
		}
		float x3 = p.x;
		float y3 = p.y;

		float x4 = 9999;	//this is a constant to draw a long horizontal line
		float y4 = p.y;

		int intersections = 0;
		for (PVector[] e : edges) {
			float x1 = e[0].x;
			float y1 = e[0].y;
			float x2 = e[1].x;
			float y2 = e[1].y;
			float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

			float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
			float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;

			if (t <= 1f && t >= 0f && u <= 1f && u >= 0f) intersections++;
		}
		return intersections%2 == 1;
	} // isInside
	
	public static PVector getHullLineIntersection(ArrayList<PVector> hull, PVector[] line) {
		// check that one point from line is inside hull
		if(!isInside(line[0], hull) && !isInside(line[1], hull)) return null;
		if(isInside(line[0], hull) && isInside(line[1], hull)) return null;
		
		PVector intersection = null;
		for(int i = 0; i < hull.size(); i++){
			int next;
			if(i == hull.size()-1) {
				next = 0;
			}else {
				next = i+1;
			}
			PVector[] edge = {hull.get(i), hull.get(next)};
			intersection = getIntersectionPoint(line,edge);
			
		}

		return intersection;

	}
	
	//=================COVERAGE===========
	public float coverage(int res) {
		float coverage = 0;
		int coveredTiles = 0;
		int w = this.e.getWidth();
		int h = this.e.getHeight();
		for(int x = 0; x < w/res; x++) {
			for(int y = 0; y < h/res; y++) {
				float sum = 0;
				float tileX = x * res + (res/2);
				float tileY = y * res + (res/2);
				PVector tileCentre = new PVector(tileX, tileY);
				for(Organism o : c.organisms){
					PVector centre = o.getHullCentroid();
					PVector[] dist = {centre, tileCentre};
					ArrayList<PVector> hull = o.getConvexHull();
					float r = meanDiameter(hull)/2;
					for(int i = 0; i < hull.size(); i++){
	                    int next;
	                    if(i == hull.size() - 1){
	                        next = 0;
	                    }else{
	                        next = i+1;
	                    }
	                    PVector[] edge = {hull.get(i), hull.get(next)};
	                    PVector intersection = Simulation.getIntersectionPoint(dist, edge);
	                    // intersection = Simulation.getIntersectionPoint(dist, edge);
	                    
	                    if(intersection != null){
	                        r = PVector.dist(centre, intersection);
	                        break;
	                    }else{
	                        r = Simulation.meanDiameter(hull)/2;
	                    }
	                    
	                }
					
					float d = PVector.dist(tileCentre, centre);
					sum+=(r/d);
				}
				if(sum > 0.85) coveredTiles++;
			}
		}
		int totalTiles = (w/res) * (h/res);
		coverage = (float)coveredTiles/(float)totalTiles;
		return coverage;
	}
	
	public float occupiedArea() {
		float occupiedArea = 0;
		for(Organism o:this.colony3d.getLayer(0).organisms) {
			occupiedArea+=o.getArea();
		}
		float totalArea = this.e.getWidth() * this.e.getHeight();
		
		
		
		float occupiedAreaRatio = occupiedArea/totalArea;
		
		return occupiedAreaRatio;
	}
	
	public float relativeCoverage() {
		float relCov = this.coverage(10) - this.occupiedArea();
		if(relCov < 0) {
			return 0;
		}else {
			return relCov;			
		}
	}
	
	
	//=================CONVEXITY===========
	//convexity ratio of single organism
	public float getConvexityRatio(Organism organism){
//		ArrayList<PVector> convexHull = organism.getConvexHull();
		float totalPerimeter = organism.getPerimeter();
		float convexPerimeter = organism.getConvexHullPerimeter();

		return convexPerimeter/totalPerimeter;
	} 

	//returns median of convexity ratio of all organisms
	public float convexity(){
//		ArrayList<Float> convexityRatios = new ArrayList<Float>();
		//for every layer
		float totalPerimeter = 0;
		float totalConvexPerimeter = 0;
		for(Colony layer : this.colony3d.layers){
			// for every organism
			for(Organism org : layer.organisms){
//				int n = layer.organisms.indexOf(org);
				float perimeter = org.getPerimeter();
				float convexPerimeter = org.getConvexHullPerimeter();
				totalPerimeter+=perimeter;
				totalConvexPerimeter+=convexPerimeter;
//				float cr = this.getConvexityRatio(org);
//				System.out.println("Org "+totalOrgs+" Conv: "+ cr);
				
//				convexityRatios.add(this.getConvexityRatio(org));
			}
		}
//	    Collections.sort(convexityRatios, null);
//	    return convexityRatios.get(convexityRatios.size()/2);
//		for (float r : convexityRatios) convexitySum+=r;
//		return convexitySum / (float) totalOrgs;
		float result = totalConvexPerimeter / totalPerimeter;
		if(Float.isNaN(result)) {
			return 1f;
		}else {
			return result;
		}
	}
	  
	//=================COMPACTNESS===========
	/**
	 * a measure of compactness using the Polsby-Popper test
	 * @return average compactness between 0 and 1, where 1 is a perfect circle (full compactness)
	 */
	public float compactness() {
		// Initialise array to store compactness scores
		ArrayList<Float> compactnessScores = new ArrayList<Float>();
		
		// for every layer
		for(Colony layer : this.colony3d.layers) {
			ArrayList<Float> layerCompactness = new ArrayList<Float>();
			// for every organism in layer
			for(Organism org : layer.organisms) {
				float compactness = (float)(4*Math.PI*org.getArea())/(org.getPerimeter()*org.getPerimeter());
				layerCompactness.add(compactness);
			}
			float layerCompactnessSum = 0;
			for(float score : layerCompactness) {
				layerCompactnessSum+=score;
			}
			compactnessScores.add(layerCompactnessSum/layerCompactness.size());
		}
		float totalSum = 0;
		for(float score : compactnessScores) {
			totalSum+=score;
		}
		return totalSum/compactnessScores.size();
	}
	
		
	//=======DISPERSION RATE OF CONSECUTIVE ANGLES========

	public ArrayList<Float> getAngles(){
		ArrayList<Float> angles = new ArrayList<Float>();
		//for every layer
		for(Colony layer : this.colony3d.layers){
			//ArrayList<ArrayList<PVector[]>> layer = this.colony3d.layers.get(i);
			//for every organism
			for (Organism org : layer.organisms){
//				ArrayList<PVector[]> org = layer.get(j);
				// for every spring
				for (int i = 0; i < org.springs.size()-1; i++){
					PVector[] spring = {org.springs.get(i).sp.loc,org.springs.get(i).ep.loc};
					PVector[] nextSpring = {org.springs.get(i+1).sp.loc, org.springs.get(i+1).ep.loc};
					float magOfSum = spring[0].dist(nextSpring[1]); //distance from sp1 to ep2
					float sumOfMags = spring[0].dist(spring[1]) + nextSpring[0].dist(nextSpring[1]); //distance from sp1 to ep2 via ep1

					// println("Layer "+i+" Edges "+k+"-"+(k+1)+" direct dist: "+magOfSum+" | length of edges: "+sumOfMags);

					// float change = 1 - (magOfSum/sumOfMags);
					angles.add(magOfSum/sumOfMags);

				}
			}
		}
		return angles;
	}
		
	public ArrayList<Float> getRealAngles(){
		ArrayList<Float> angles = new ArrayList<Float>();
		// a.b/mag(a)*mag(b)
		for(Colony layer:this.colony3d.layers) {
			for(Organism org : layer.organisms) {
				for(int i = 0; i < org.springs.size()-1;i++) {
					Spring s = org.springs.get(i);
					Spring nexts = org.springs.get(i+1);
					PVector v1 = PVector.sub(s.sp.loc, s.ep.loc);
					PVector v2 = PVector.sub(nexts.ep.loc, nexts.sp.loc);
					float cosA = PVector.dot(v1, v2)/(v1.mag() * v2.mag());
					float a = (float) Math.acos((double) cosA);
					angles.add(a);
				}
			}
		}
			
		return angles;
	}
			
	public float angleDispersionCoefficient(){
		ArrayList<Float> a = this.getRealAngles();
		if(a.size() < 1) return 0f;
		Collections.sort(a, null);
		float q1 = a.get((int)(a.size()/4));
		float q3 = a.get((int)(3 * (a.size()/4)));
			
		float result = (q3-q1)/(q3+q1);
		if(Float.isNaN(result)) {
			return 0f;
		}else {
			return result;
		}
	} //angle dispersion coefficient
	
	/**
	 * Splitting Score
	 * the more a colony splits, the higher its score goes
	 * splitting constant set at 0.1
	 * @return float
	 */
	public float splittingScore() {
		double s = (double)this.c.getSplits();
		double splitC = Math.pow(0.1, (1/(s+1)));
		return (float)splitC;
	}//splittingScore
		  
	public float complexity() {
		float complexity = 0;
	  
		// get convexity
		// low is good
		float convexity = this.convexity();
		complexity+=(1 - convexity);
			  
		// get angle dispersion coefficient
		// high is good
		float angleDispersion = this.angleDispersionCoefficient();
		complexity += angleDispersion;
		
		float splitScore = this.splittingScore();
		complexity+=splitScore;
		return complexity/3;
	}
	
	

	public void setCostOfLiving(float cost) {
		this.c.setCostOfLiving(cost);
	}
		  
	public void setSpringRLMult(float rl) {
		this.c.setSpringRLMult(rl);
	}
}




