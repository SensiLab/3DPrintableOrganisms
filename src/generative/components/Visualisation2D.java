package generative.components;

import generative.components.CommonConsts;
import java.util.Random;
import processing.core.PApplet;
import processing.core.PVector;
import processing.data.Table;
import processing.data.TableRow;
import peasy.*;

public class Visualisation2D extends PApplet {
	
	PeasyCam camera;
	
	Random generator;
	
//	long envSeed = CommonConsts.getEnvSeed(25);
	
	double[] genes;
	Chromosome chr;
	// Environment size set in millimeters
	Environment e;
	Colony c;
	Simulation s;
	Colony3D col;
	int counter = 0;
	
	public void settings() {
		size(600,600, P3D);
	

		
		
//		s.generate();
//		col = s.colony3d;
//		noLoop();
	}
	
//	@Override
	public void setup() {
//		camera = new PeasyCam(this, 300, -300, 0, 675);
		camera = new PeasyCam(this, 300, 300, 0, 900);
	    camera.rotateZ(PI/4);
	    camera.rotateX(-PI/4);
	    camera.pan(0,-225);
	    
	    this.generator = new Random();
	    
//	    for(int i = 0; i < 520; i++) {
//	    	e.updateNutrientField(c);
//	    	c.update();
//	    	if(i >=20) {
//	    		col.addLayer(c);
//	    	}
//	    }
	    
	    String folderPath = "../out/evo_output/IEEE_Model-Alternatives_001/";
	    String alt = "02_No-Physics";
	    /*
	     * Alternatives consist in:
	     * 01 - Full model -> no modification to the original model
	     * 02 - No physics -> comment out Organism 877, 899, 916, 917 and Colony 196
	     * 02 - No energy consumption -> comment out Cell 137 and Organism 885
	     * 03 - No organism division -> comment out Colony 235-246, uncomment Colony 248-253
	     * 04 - No metabolic system -> comment out Colony 194
	     */
	    
	    // initialise table
	    Table dataOut = new Table();
        dataOut.addColumn("iteration");
        dataOut.addColumn("envrandseed");
        dataOut.addColumn("genes");
        dataOut.addColumn("printability");
        dataOut.addColumn("complexity");
        dataOut.addColumn("area_mean");
        dataOut.addColumn("area_dev");
	    
	    // repeat 100 times
	    for(int i = 0; i < 100; i++){
			e = new Environment();
			e.setRandomSeed(CommonConsts.getRandomSeed());
	    	
	    	
	    	// Generate random genes
	    	double[] genes = new double[5];
	    	for(int g = 0; g < genes.length; g++) {
	    		genes[g] = generator.nextDouble();
	    	}
	    	
	    	chr = new Chromosome(genes);
	    	s = new Simulation(chr, e, 500, 20);	    	
	    	
	    	
	    	// generate individual
	    	s.generate();
	    	
	    	// evaluate individual
	    	
	    	// complexity
	    	float complexity = s.complexity();
//	    	// printability
	    	float printability = s.printability();
	    	//area mean
	    	float areamean = s.areaMean();
//	    	// areadev
	    	float areadev = s.areaStandardDev();
//	    	
	    	println("Complexity:\t"+complexity+"\n"+
	    			"Printability:\t"+printability+"\n"+
	    			"Area mean:\t"+areamean+"\n"+
	    			"Area deviation:\t"+areadev+"\n");
//	    	
	    	// store data in table
	    	TableRow newRow = dataOut.addRow();
            newRow.setInt("iteration",(i+1));
            newRow.setLong("envrandseed", e.getRandomSeed());
            newRow.setString("genes",chr.chrToString());
            newRow.setFloat("printability", printability);
            newRow.setFloat("complexity", complexity);
	    	newRow.setFloat("area_mean", areamean);
            newRow.setFloat("area_dev", areadev);
	    	
	    	// store image
          //======RENDER=======
            stroke(30, 150);
            strokeWeight(0.5f);
            Colony3D col = s.getColony3D();
            Colony topC = col.getLayer(col.size() - 1);
            for(int j = 0; j < col.size(); j+=3){
                Colony c = col.getLayer(j);
                displayColony(c, j);

            }
            saveFrame(folderPath+"imgs/"+alt+"/iteration_"+String.format("%04d", i)+".png");

//             reset background
            background(220);
            
            println("Full Model, iteration: "+i);
//             break;
	    	
	    }
	    
	    //save table
	    saveTable(dataOut, folderPath+"data/"+alt+".csv");	    
		exit();
	}
	
//	public void draw() {
//		if(counter > 520) {
//			fill(0);
//			text("Counter:\t"+(counter-1),10,15);
//			text("Env Size:\t"+e.getWidth(),10,30); 
//			noLoop();
//		}
//		background(220);
//        stroke(20,50);
//        strokeWeight(1f);
//        for(int i = 0; i < col.getLayers().size();i++) {
//        	if(i%3==0) {
//        	Colony cc = col.getLayer(i);
//        	displayColony(cc, i);
//        	}
//        }
//        displayEnvironment(e);
//        e.updateNutrientField(this.c);
////        for (int i = 0; i < this.colonyLayers + this.colonyWarmupLayers; i++) {
//			
//		      
//		this.c.update();
//		counter++;

//	}
	
	public void displayColony(Colony _c){
		// Scale colony based on drawing size
		
	    for(Organism o : _c.organisms){
	        for(Spring s : o.getSprings()){
	        	float spX = s.sp.loc.x * width/_c.e.width;
	        	float spY = s.sp.loc.y * height/_c.e.height;
	        	float epX = s.ep.loc.x * width/_c.e.width;
	        	float epY = s.ep.loc.y * height/_c.e.height;
	            line(spX,spY,epX,epY);
//	            ellipse(spX, spY, width * .01f, height * .01f);
	        }
	        
	    }
	}
	
	public void displayColony(Colony _c, float layerZ){
		// Scale colony based on drawing size
		
	    for(Organism o : _c.organisms){
	        for(Spring s : o.getSprings()){
	        	float spX = s.sp.loc.x * width/_c.e.width;
	        	float spY = s.sp.loc.y * height/_c.e.height;
	        	float spZ = layerZ;
	        	float epX = s.ep.loc.x * width/_c.e.width;
	        	float epY = s.ep.loc.y * height/_c.e.height;
	        	float epZ = layerZ;
	            line(spX, spY, spZ, epX, epY, epZ);
//	            ellipse(spX, spY, width * .01f, height * .01f);
	        }
	        
	    }
	}
	
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
