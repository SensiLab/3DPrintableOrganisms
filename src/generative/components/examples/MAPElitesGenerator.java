package generative.components.examples;

import java.io.File;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import processing.core.*;
import processing.data.JSONArray;
import processing.data.JSONObject;
import generative.components.*;
import peasy.*;

public class MAPElitesGenerator extends PApplet{
	
	PeasyCam camera;
	
	Random generator;
	
	String destinationFolder = "/Users/ccru0002/Documents/Work/01_3D-printable-organisms/CodeDev/MAP_Elites_Tests/temp/";
	
	double[] parent = {-1, -1, -1, -1, -1};
	Chromosome chr;
	// Environment size set in millimeters
	Environment e;
	Cell c;
	Cell c2;
	ArrayList<Cell> cells;
	Organism org;
	Colony col;
	Simulation sim;
	Colony3D col3d;
	
	int colonyLayers;
	int colonyWarmup;
	
	ArrayList<PVector> locations;
	
	JSONArray json;
	
	int lambda = 1;
	
	double mutationRate = 0.5;
	double mutationFactor = 0.01;
	
	int foodSources = 5;
	int envSeed = 25;
	
	public MAPElitesGenerator(String[] args) {
		for(String s : args) {
			System.out.println(s);
			String[] arg = s.split("=");
			switch(arg[0]) {
			case "-p":
				String parentStr = String.valueOf(arg[1]);
				String[] genesStr = parentStr.replace("(","").replace(")","").split(",");
				if(genesStr.length != parent.length) {
					System.out.println("Please provide parent data using -x= followed by 5 doubles separated by commas.\n"
							+ "Example: -p=(double,double,double,double,double)");
					System.exit(0);
				}
				for(int i = 0; i < genesStr.length; i++) {
					parent[i] = Double.parseDouble(genesStr[i]);
				}
				break;
			case "-mr":
				mutationRate = Double.valueOf(arg[1]);
				break;
			case "-mf":
				mutationFactor = Double.valueOf(arg[1]);
				break;
			case "-es":
				int es = Integer.valueOf(arg[1]);
				if(es < 0) {
					System.out.println("-es takes values between 0 and 49. The random seed for the environment  has been set to 0");
					System.exit(0);
				}else if(es > 49) {
					System.out.println("-es takes values between 0 and 49. The random seed for the environment  has been set to 49");
					System.exit(0);
				}else {
					envSeed = es;
				}
				break;
			case "-ef":
				foodSources = Integer.valueOf(arg[1]);
				break;
			case "-l":
				lambda = Integer.valueOf(arg[1]);
				break;
			case "-df":
				destinationFolder = String.valueOf(arg[1]);
				break;
			case "-ol":
				locations = new ArrayList<PVector>();
				String orgLocsString = String.valueOf(arg[1]);
				String[] orgLocsArray = orgLocsString.split(";");
				for(String locationStr : orgLocsArray) {
					String[] vectorStr = locationStr.replace("(", "").replace(")", "").split(",");
					if(vectorStr.length == 2) {
						float locX = Float.parseFloat(vectorStr[0]);
						float locY = Float.parseFloat(vectorStr[1]);
						PVector orgLoc = new PVector(locX, locY);
						locations.add(orgLoc);
					}
				}
				break;
			case "-cw":
				colonyWarmup = Integer.valueOf(arg[1]);
				break;
			case "-cl":
				colonyLayers = Integer.valueOf(arg[1]);
				
			default:
				break;
			}
			
		}

	}
	
	public void settings() {
		size(600,600,P3D);
		
	}
	
	public void setup() {
		
		background(220);
		// Camera setup
		camera = new PeasyCam(this, width/2, height/2, 0, 1.5*height);
	    camera.rotateZ(PI/4);
	    camera.rotateX(-PI/4);
	    camera.pan(0,-width/2.35);
	    
	    // Initialise JSON file
	    json = new JSONArray();
	    
	    // generate timestamp
	    DateTimeFormatter timestampFormat = DateTimeFormatter.ofPattern("yyyyMMdd-HH_mm_ss");
	    String timestamp = LocalDateTime.now().format(timestampFormat);
	    
	    // check that imgs folder is empty
	    File img_folder = new File(destinationFolder+"imgs/");
	    for(File f : img_folder.listFiles()) {
	    	if(!f.isDirectory()) f.delete();
//	    	System.out.println(f.getName());
	    	
	    }
	 // check that ptData folder is empty
	    File ptData_folder = new File(destinationFolder+"ptData/");
	    for(File f : ptData_folder.listFiles()) {
	    	if(!f.isDirectory()) f.delete();
//	    	System.out.println(f.getName());
	    	
	    }
	    
//	    for(File file : img_folder.listFiles()) {
//	    	println(file);
//	    }
	    
//	    System.out.println(Arrays.toString(locations));
	    
	    // generate individuals
	    for(int i=0; i < lambda; i++) {
//	    	System.out.println("Parent[0]: "+parent[0]);
	    	generator = new Random();
	    	double[] genes = new double[parent.length];
	    	//check parent: if parent[0] == -1, generate random
	    	for(int j=0; j< parent.length; j++) {
//	    		System.out.println("Gene "+j+ " before: "+parent[j]);
//	    		System.out.println("-------------------");
	    		if(parent[0] == -1) {
	    			
	    			genes[j] = generator.nextDouble();
	    		}else {
	    			
	    			// flip a coin
	    			double m = generator.nextDouble();
	    			// change allele j if the coin toss yields a number smaller than the mutation rate
	    			if(m <= mutationRate) {
	    				double mutation = generator.nextDouble() * 2 - 1;
//	    				System.out.println("Mutation: "+mutation);
	    				genes[j]+=(parent[j]+ (mutation * mutationFactor));
	    			}else {
	    				genes[j] = parent[j];
	    			}
	    		}
//	    		System.out.println("Gene "+j+ " after: "+genes[j]);
    		}
	    	System.out.println("Genes: "+Arrays.toString(genes));
	    	// generate chromosome
	    	Chromosome chr = new Chromosome(genes);
	    	
	    	//initialise environment
	    	e = new Environment();
			e.setRandomSeed(CommonConsts.getEnvSeed(envSeed));
			e.setFoodSources(foodSources);
			
			
			if(colonyLayers == 0) colonyLayers = 500;
			if(colonyWarmup == 0) colonyWarmup = 20;
			
			// start simulation
			// check for locations
			boolean isComplete = false;
			int attempts = 0;
			int maxAttempts = 50;
			while(!isComplete && attempts < maxAttempts) {
//				System.out.println("attempt "+(attempts+1));
				if(locations != null && locations.size() > 0) {
					sim = new Simulation(chr, e, colonyLayers, colonyWarmup, locations);
				}else {
					sim = new Simulation(chr, e, colonyLayers, colonyWarmup);
				}
				sim.generate();
//				System.out.println("Colony "+i+" is complete = "+sim.isComplete());
				if(sim.isComplete() < 1) {
					double[] oldGenes = chr.getGenes2();
					double[] newGenes = new double[oldGenes.length];
//					System.out.println("Start new genes");
					//change genes
					for(int j = 0; j < oldGenes.length; j++) {
						// flip a coin
						double m = generator.nextDouble();
						// change allele j if the coin toss yields a number smaller than the mutation rate
						if(m <= mutationRate) {
							double mutation = generator.nextDouble() * 2 - 1;
//							System.out.println("Mutation: "+mutation);
							newGenes[j] = oldGenes[j] + (mutation * mutationFactor);
						}else {
							newGenes[j] = oldGenes[j];
						}
						if(newGenes[j] < 0) newGenes[j] = 0;
						if(newGenes[j] > 0) newGenes[j] = 1;
					}
//					System.out.println("New genes generated");
//					System.out.println(Arrays.toString(newGenes));
					//apply to chromosome
					chr = new Chromosome(genes);
					
				}else {
					isComplete = true;
				}
				attempts++;
			}
			
			
			
			
			
			//calculate printability
			float p = sim.printability();
			float printability = MathLib.sigmoid(p, 1f, 0.8f, 30f);
			

			float cr = sim.coverageRatio();
//			float  = MathLib.sigmoid(p, 1f, 0.8f, 30f);
			
			//calculate convexity
			float convexity = sim.convexity();
			
			//calculate angle dispersion coefficient
			float angleDisp = sim.angleDispersionCoefficient();
			
			//calculate edge length dispersion coefficient
			float edgeLenDisp = sim.lengthDispersionCoefficient();
			
			//calculate total splits score
			float splits = sim.totalSplits();
			
			//calculate internal volume
			float internalVolume = sim.getVolume();
			
			//calculate area standard deviation
			float areaStDev = sim.areaStandardDev();
			
			//calculate first split layer
			float firstSplit = sim.firstSplit();
			
			//calculate average growth rate
			float avgGR = sim.growthRate();
			
			//calculate cell growth rate
			float cellGR = sim.cellGrowthRate();
			
			//calculate energy rate of change
			float energyChange = sim.energyChange();
			
			//calculate hamming distance
			float hammingDistance = sim.hammingDistance(60);
			
			float overlaps = sim.overlaps();
			
			//generate image name
			String img_name = timestamp + "_" +String.format("%03d", i);
			
			// ==========DRAW=============
		    // generate image
			stroke(0,100);
		    CommonFunctions.displayColony3D(this, sim.getColony3D(), 3);
		    
		    //save image
		    saveFrame(destinationFolder+"imgs/"+img_name+".png");
		    
		    //save point data
		    sim.saveColony3D(destinationFolder+"ptData/"+img_name);
		    
		    //clear background
		    background(220);
			
			println(img_name + ": " + printability + " (" + p +") " + timestamp);
			
			//initialise JSONObject
			JSONObject ind = new JSONObject();
			
			//identification of object
			ind.setString("id", img_name);
			ind.setString("genes", Arrays.toString(genes));
			ind.setInt("env_rand_seed", envSeed);
			
			//features
			ind.setFloat("printability", printability);
			ind.setFloat("coverage_ratio", cr);
			ind.setFloat("convexity", convexity);
			ind.setFloat("angle_dispersion", angleDisp);
			ind.setFloat("edge_length_dispersion", edgeLenDisp);
			ind.setFloat("splits", splits);
			ind.setFloat("internal_volume", internalVolume);
			ind.setFloat("area_stdev", areaStDev);
			ind.setFloat("first_split", firstSplit);
			ind.setFloat("growth_rate", avgGR);
			ind.setFloat("cell_growth_rate", cellGR);
			ind.setFloat("energy_change_rate", energyChange);
			ind.setFloat("hamming_distance", hammingDistance);
			ind.setFloat("overlaps", overlaps);
			
			
			
	    	
	    	json.setJSONObject(i, ind);
	    	
	    	
	    }
	    // save json file
	    saveJSONArray(json, destinationFolder+"data.json");
	    


	    //exit	    
	    exit();
		
	}
	
	
	public static void main(String[] args) {
		String[] processingArgs = {"VisTests"};
		MAPElitesGenerator vis2D = new MAPElitesGenerator(args);
		PApplet.runSketch(processingArgs, vis2D);
		
	}
}
