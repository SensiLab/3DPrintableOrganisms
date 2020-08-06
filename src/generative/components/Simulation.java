package generative.components;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;

//import fr.inria.optimization.cmaes.IntDouble;
import processing.core.*;


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
	int colonySize;
	int colonyWarmup;
	int envRes = 20;
	ArrayList<ArrayList<ArrayList<OVector[]>>> pointCloud;
	Colony3D colony3d;
	Colony c;
	float printabilityThreshold = 0.2f;
	float minPrintDiameter = 7.5f;
	float outputSize = 120f;
	float outputLayerHeight = 0.2f;
	
	float angleDispersion;
	float printability;

	// MAIN CONSTRUCTOR	
	public Simulation(Chromosome _c, Environment _e, int _colonySize, int _colonyWarmup){
		this.e = _e;
		this.chromosome = _c;
		this.colonySize = _colonySize;
		this.colonyWarmup = _colonyWarmup;

		this.c = new Colony(this.chromosome, this.e);

		this.pointCloud = new ArrayList<ArrayList<ArrayList<OVector[]>>>();
		this.colony3d = new Colony3D();

		// println("Simulation initialised\n"+" Env: "+this.e+"\n");

	} //constructor
	


	//======================GENERATE==================
	// Method to generate the geometric data of an object
	public void generate() {
		
		//copy environment
		Environment env = new Environment(this.e);
		this.c.setEnvironment(env);

//		this.pointCloud.clear();
		this.colony3d.clear();

		for (int i = 0; i < this.colonySize + this.colonyWarmup; i++) {
			env.updateNutrientField(this.c);
		      
			this.c.update();

			if (i >= this.colonyWarmup) {
				
				Colony layer;
				// Copy colony
				if (this.c.organisms.size() != 0) {
					layer = new Colony(this.c);
					layer.setZPosition((float) (i - this.colonyWarmup));
				} else {
					layer = new Colony();
				}
				
				// Add layer to colony 3d
				this.colony3d.addLayer(layer);

			}
		}
	}
	
	public Colony3D getColony3D() {
		return this.colony3d;
	}
	
	float diameterFactor(float d, float threshold){
		if(d > 2*threshold) return 1;
		if(d >= threshold) return MathLib.map(d,threshold,(2*threshold),0.0f,1.0f);
		return 0;
	}
	
	/**
	 * Utility methods for the evaluation of 3d colonies
	 */
	
	/**
	 * Method to find the closest spring to a given cell 
	 * @param pt (the cell)
	 * @param layer (the colony in which the spring is)
	 * @return array with the index of the organism and the spring in the colony
	 */
	static int[] closestEdge(OVector pt, Colony layer){
		  float minDist = 99999999;
		  int closestOrg = -1;
		  int closestEdge = -1;
		  for(int i = 0; i < layer.organisms.size(); i++){
			  Organism organism = layer.organisms.get(i);
			  for(int j = 0; j < organism.springs.size(); j++){
				  OVector[] edge = {organism.springs.get(j).sp.loc, organism.springs.get(j).ep.loc};
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
	}
	
	
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
		
//		long startScale = System.nanoTime();
//		 scale 3d colony
		Colony3D scaledColony3d = this.colony3d.scale(this.outputSize, this.outputLayerHeight, this.e);
//		long endScale = System.nanoTime();
//		System.out.println("Scale colony: "+ ((endScale - startScale)/1000000));
		
		// Initialise array to store scores per vertex
		// The array is nested in order to keep track of scores of vertices in layers below
		ArrayList<ArrayList<ArrayList<Float>>> colonyScores = new ArrayList<ArrayList<ArrayList<Float>>>();
		
		// add scores for layer 0
		ArrayList<ArrayList<Float>> layer0Scores = new ArrayList<ArrayList<Float>>(); 
		// for every organism in layer 0
//		long startScoresLayer0 = System.nanoTime();
		for(Organism o : scaledColony3d.layers.get(0).organisms) {
			ArrayList<Float> orgLayer0Scores = new ArrayList<Float>();
			
			ArrayList<OVector> cHull = o.getConvexHull();
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
//		long endScoresLayer0 = System.nanoTime();
//		System.out.println("Scores Layer 0: "+ ((endScoresLayer0 - startScoresLayer0)/1000000));
		// scores for layer 0 added
		
		// add scores for all other layers
		for(int i = 1; i < scaledColony3d.layers.size(); i++) {
//			long startScoresLayer = System.nanoTime();
			ArrayList<ArrayList<Float>> layerScores = new ArrayList<ArrayList<Float>>();
			Colony layer = scaledColony3d.layers.get(i);
			Colony layerBelow = scaledColony3d.layers.get(i-1);
			
			//for every organism in layer
			for(int j = 0; j < layer.organisms.size(); j++) {
				Organism o = layer.organisms.get(j);
				
				// calculate the diameter factor here
				ArrayList<OVector> cHull = o.getConvexHull();
				float hullMinDiameter = minDiameter(cHull);
				float diameterFactor = diameterFactor(hullMinDiameter, this.minPrintDiameter);				
				
				ArrayList<Float> orgScores = new ArrayList<Float>();
				for(int k = 0; k < o.springs.size(); k++) {
					Spring s = o.springs.get(k);
					float edgeLen = s.getLen();

					// calculate the score for this edge here
					
					//check support for start point
					OVector sp = s.sp.loc;
					// Get the index of the organism and spring in the colony layerBelow
					int[] spClosestEdgeCoordinates = closestEdge(sp, layerBelow);
					// get the support score of the edge closest to sp on layer below
					float closestEdgeScoreSp = colonyScores.get(i-1).get(spClosestEdgeCoordinates[0]).get(spClosestEdgeCoordinates[1]);
					// Get the actual spring
					Spring closestSpringSp = layerBelow.organisms.get(spClosestEdgeCoordinates[0]).springs.get(spClosestEdgeCoordinates[1]);
					// make a reference to the positions of the closest spring's sp and ep
					OVector[] closestEdgeSp = {closestSpringSp.sp.loc, closestSpringSp.ep.loc};
					
					float distSP = distPointLine(sp, closestEdgeSp);
					float scoreSp = printabilityScore(distSP,this.printabilityThreshold) * closestEdgeScoreSp;
					
					//check support for mid point
					OVector mp = OVector.getMidPoint(s.sp.loc, s.ep.loc);
					int[] mpClosestEdgeCoordinates = closestEdge(mp, layerBelow);
					// get the support score of the edge closest to mp on layer below
					float closestEdgeScoreMp = colonyScores.get(i-1).get(mpClosestEdgeCoordinates[0]).get(mpClosestEdgeCoordinates[1]);
					Spring closestSpringMp = layerBelow.organisms.get(mpClosestEdgeCoordinates[0]).springs.get(mpClosestEdgeCoordinates[1]);
					OVector[] closestEdgeMp = {closestSpringMp.sp.loc, closestSpringMp.ep.loc};
					float distMP = distPointLine(mp, closestEdgeMp);
					float scoreMp = printabilityScore(distMP,this.printabilityThreshold) * closestEdgeScoreMp;
					
					//check support for end point
					OVector ep = s.ep.loc;
					int[] epClosestEdgeCoordinates = closestEdge(ep, layerBelow);
					// get the support score of the edge closest to ep on layer below
					float closestEdgeScoreEp = colonyScores.get(i-1).get(epClosestEdgeCoordinates[0]).get(epClosestEdgeCoordinates[1]);
					Spring closestSpringEp = layerBelow.organisms.get(epClosestEdgeCoordinates[0]).springs.get(epClosestEdgeCoordinates[1]);
					OVector[] closestEdgeEp = {closestSpringEp.sp.loc, closestSpringEp.ep.loc};
					float distEP = distPointLine(ep, closestEdgeEp);
					float scoreEp = printabilityScore(distEP,this.printabilityThreshold) * closestEdgeScoreEp;
					
					float edgeScore = (diameterFactor * ((scoreSp + (2 * scoreMp) + scoreEp )/4));
					
//					System.out.println("Edge len: "+edgeLen+" Edge score: "+edgeScore+" (diameter factor = "+diameterFactor+")");
					
					// add scores to overall counters
					overallPerimeter+=edgeLen;
					supportedPerimeter+=(edgeScore*edgeLen);
					orgScores.add(edgeScore);
				}
				layerScores.add(orgScores);
			}
			colonyScores.add(layerScores);
//			long endScoresLayer = System.nanoTime();
//			System.out.println("Scores Layer "+i+": "+ ((endScoresLayer - startScoresLayer)/1000000)+" layer len: "+layer.getTotalSprings());
		}
		return supportedPerimeter/overallPerimeter;
	}
	
	
	/**
	 * Calculates the distance between a point and a line on the xy plane
	 * @param pt
	 * @param ln
	 * @return distance as a float
	 */
	
	static float distPointLine(OVector pt, OVector[] ln) {
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
	 * @param line as an OVector[] with a start point and an end point
	 * @param array of points
	 * @return the distance between the line and the point farthest away from it
	 */
	static float maxLinePointDist(OVector[] ln, ArrayList<OVector> points){
		float maxDist = 0;
		for(OVector pt : points){
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
	public static float minDiameter(ArrayList<OVector> organism){
		  // println("minDiameter called");
		  ArrayList<OVector[]> edges = new ArrayList<OVector[]>();
		  for(int i = 0; i < organism.size() - 1; i++) {
			  OVector[] edge = {organism.get(i), organism.get(i+1)};
			  edges.add(edge);
		  }
		  OVector[] lastEdge = {organism.get(organism.size() - 1), organism.get(0)};
		  edges.add(lastEdge);
		  
		  ArrayList<Float> diameters = new ArrayList<Float>();
		  for(OVector[] l : edges){
		    float d = maxLinePointDist(l,organism);
		    diameters.add(d);
		    // println(d);
		  }
		  Collections.sort(diameters, null);
		  // println(diameters);
		  return diameters.get(0);

	}
	
	//=================CONVEXITY===========
	//convexity ratio of single organism
	public float getConvexityRatio(Organism organism){
//		ArrayList<OVector> convexHull = organism.getConvexHull();
		float totalPerimeter = organism.getPerimeter();
		float convexPerimeter = organism.getConvexHullPerimeter();

		return convexPerimeter/totalPerimeter;
	} 

	//returns median of convexity ratio of all organisms
	public float convexity(){
		ArrayList<Float> convexityRatios = new ArrayList<Float>();
		//for every layer
		for(Colony layer : this.colony3d.layers){
			// for every organism
			for(Organism org : layer.organisms){
				convexityRatios.add(this.getConvexityRatio(org));
			}
		}
//	    Collections.sort(convexityRatios, null);
//	    return convexityRatios.get(convexityRatios.size()/2);
		float convexitySum = 0;
		for (float r : convexityRatios) convexitySum+=r;
		return convexitySum / (float) convexityRatios.size();
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
//		  float dispersionRateAngles(){
//		    ArrayList<Float> angleDifference = new ArrayList<Float>();
//		    //for every layer
//		    for(int i = 0; i < this.pointCloud.size(); i++){
//		      ArrayList<ArrayList<OVector[]>> layer = this.pointCloud.get(i);
//		      //for every organism
//		      for (int j = 0; j < layer.size(); j++){
//		        ArrayList<OVector[]> org = layer.get(j);
//		        for (int k = 0; k < org.size()-1; k++){
//		          OVector[] spring = org.get(k);
//		          OVector[] nextSpring = org.get(k+1);
//		          float magOfSum = spring[0].dist(nextSpring[1]); //distance from sp1 to ep2
//		          float sumOfMags = spring[0].dist(spring[1]) + nextSpring[0].dist(nextSpring[1]); //distance from sp1 to ep2 via ep1
//
//		          // println("Layer "+i+" Edges "+k+"-"+(k+1)+" direct dist: "+magOfSum+" | length of edges: "+sumOfMags);
//
//		          float change = 1 - (magOfSum/sumOfMags);
//		          angleDifference.add(change);
//
//		        }
//		      }
//		    }
//		    // println("Angle Variance: " + Math.variance(angleDifference) + " Angle diff avg: "+Math.average(angleDifference));
//		    // println(millis());
//		    return Math.variance(angleDifference);
//		  } //dispersion rate angles

		public ArrayList<Float> getAngles(){
		    ArrayList<Float> angles = new ArrayList<Float>();
		    //for every layer
		    for(Colony layer : this.colony3d.layers){
//		      ArrayList<ArrayList<OVector[]>> layer = this.colony3d.layers.get(i);
		      //for every organism
		      for (Organism org : layer.organisms){
//		        ArrayList<OVector[]> org = layer.get(j);
		    	// for every spring
		        for (int i = 0; i < org.springs.size()-1; i++){
		          OVector[] spring = {org.springs.get(i).sp.loc,org.springs.get(i).ep.loc};
		          OVector[] nextSpring = {org.springs.get(i+1).sp.loc, org.springs.get(i+1).ep.loc};
		          float magOfSum = spring[0].dist(nextSpring[1]); //distance from sp1 to ep2
		          float sumOfMags = spring[0].dist(spring[1]) + nextSpring[0].dist(nextSpring[1]); //distance from sp1 to ep2 via ep1

		          // println("Layer "+i+" Edges "+k+"-"+(k+1)+" direct dist: "+magOfSum+" | length of edges: "+sumOfMags);

		          // float change = 1 - (magOfSum/sumOfMags);
		          angles.add(magOfSum/sumOfMags);

		        }
		      }
		    }
		    // println("Angle Variance: " + Math.variance(angleDifference) + " Angle diff avg: "+Math.average(angleDifference));
		    // println(millis());
		    return angles;
		  }

		  public float angleDispersionCoefficient(){
		    ArrayList<Float> a = this.getAngles();
		    Collections.sort(a, null);
		    float q1 = a.get((int)(a.size()/4));
		    float q3 = a.get((int)(3 * (a.size()/4)));
		    return (q3-q1)/(q3+q1);
		  }
//
//		  int[] angleFrqDistro(int buckets){
//		    int[] f = new int[buckets];
//		    ArrayList<Float> a = this.getAngles();
//		    for(int i = 0; i < a.size(); i++){
//		      int b = int(a.get(i)/(1.0/buckets));
//		      if(b==buckets) b-=1;
//		      f[b]++;
//		    }
//		    //initialise buckets
//		    return f;
//		  }


//		  
//
//		  //Solidity
//		  float getSolidtyRatio(ArrayList<OVector[]> organism){
//		    ArrayList<OVector[]> convexHull = getConvexHull(organism);
//		    float totalArea = organismArea(organism);
//		    float convexArea = organismArea(convexHull);
//
//		    return totalArea/convexArea;
//		  }
//
//		  float solidityMedian(){
//		      ArrayList<Float> solidityRatios = new ArrayList<Float>();
//		    //for every layer
//		    for(ArrayList<ArrayList<OVector[]>> layer : this.pointCloud){
//		      // for every organism
//		      for(ArrayList<OVector[]> org : layer){
//		        solidityRatios.add(getSolidtyRatio(org));
//		      }
//		    }
//		    Collections.sort(solidityRatios, null);
//		    return solidityRatios.get(solidityRatios.size()/2);
//		  }
//
//		  //===================CIRCULAR VARIANCE=============
//		  ArrayList<Float> getRadii(ArrayList<ArrayList<ArrayList<OVector[]>>> pointCloud){
//		  ArrayList<Float> radii = new ArrayList<Float>();
//		  for(ArrayList<ArrayList<OVector[]>> layer:pointCloud){
//		    for(ArrayList<OVector[]> org:layer){
//		      OVector centroid = getCentroid(org);
//		      for(OVector[] l : org){
//		        radii.add(l[0].dist(centroid));
//		      }
//		    }
//		  }
//		  return radii;
//		}
//
//		float getCircularVariance(){
//		  ArrayList<Float> radii = getRadii(this.pointCloud);
//		  return Math.variance(radii);
//		}
//
//		  //====================TOTAL FITNESS=================
//		  float getFitness() {
//		    //use sigma function to add printability score
//		    //The values used for the fuction are: m = 1.0, mp = 0.7, k = 15
//		    // if(this.pointCloud.size() < 1) return 0;
//		    if(this.lifeSpan() <= .1) return 0;
//		    float p = Math.sigmoid(this.printability, 1.0, 0.7, 15)/4;
//		    float c = this.angleDispersionCoefficient()/4;
//		    float l = this.edgeLenDispersionCoefficient()/4;
//		    float a = this.areaDispersionCoefficient()/4;
//		    return (p + c + l + a);
//		  }
}



