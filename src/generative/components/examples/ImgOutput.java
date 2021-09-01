package generative.components.examples;

import generative.components.CommonConsts;
import processing.core.*;
import processing.core.PApplet;
import processing.data.*;
import peasy.*;

import generative.components.*;


public class ImgOutput extends PApplet{
	
	PeasyCam camera;
	
	Environment e;
	
	//=====VIZ COLOURS===========
	int bg = 220;
	int bg_data = 200;

	//=====Environment===========
	int envWidth = 600;
	int envHeight = 600;
	int res = 15;
	float drag = 0.4f;
	int fss = 5;
	int fsSizeMin = 2;
	int fsSizeMax = 5;
	float foodGrowthRate = 5;
	float foodDecayRate = 0.2f;
	long ranSeed;
	
	//======DATA STORAGE=========
	float[] adisp;
	float[] conv;
	float[] print;
	double[][] chromosomes;
	
	
	public void settings() {
		size(envWidth,envHeight,P3D);
		
	}
	
	public void setup() {
		camera = new PeasyCam(this, 300, -300, 0, 675);
	    camera.rotateZ(PI/4);
	    camera.rotateX(-PI/4);
	    camera.pan(0,-225);
		
		background(220);
		
		String folderPath = "../out/evo_output/IEEE_Refine_014/";
		int pathLen = folderPath.length();
	    String splitForPrefix = folderPath.substring(pathLen-4,pathLen);
	    println("Listing all filenames in a directory: ");
	    String[] filenames = listPaths(folderPath);
	    for(String f : filenames) {
	    	f = f.replace("\\","/");
	    	
	    	if(f.contains("GenesData")) {
//	    	if(f.contains("GenesData") && f.contains("025_2021")) {
	    		println(f);
	    		
	    		String[] prfxs = f.split("/");
	            int prfxsLen = prfxs.length;
	    		
	            String prefix = prfxs[prfxsLen-1].split("_")[0];
	            int prefixInt = Integer.parseInt(prefix);
	            Long rSeed = CommonConsts.getEnvSeed(prefixInt);
	            println("Prefix: "+prefixInt);
	            println("Env Seed: "+rSeed);
	            
	            Table t = loadTable(f);
	            
	            
	            Table dataOut = new Table();
	            dataOut.addColumn("iteration");
	            dataOut.addColumn("genes");
	            dataOut.addColumn("total_fitness");
	            dataOut.addColumn("printability");
	            dataOut.addColumn("complexity");
	            
	            for(int i = 1; i < t.getRowCount() ; i++){
	            	TableRow r = t.getRow(i);
	                int iteration = r.getInt(0);
	                float fitness = r.getFloat(4);
	                String chrString = r.getString(5).replace("(", "").replace(")","").replace(";",",");
	                String[] chrStringArray = chrString.split(",");
	                
	                double[] genes = new double[chrStringArray.length];
	                for(int k = 0; k < chrStringArray.length; k++){
	                        genes[k] = Double.parseDouble(chrStringArray[k]);
	                }

	                Chromosome chr = new Chromosome(genes);
	                Environment e = new Environment();
	                e.setRandomSeed(rSeed);
	                Simulation s = new Simulation(chr, e, 500, 20);
	                
	                s.generate();
	                
	              //=======RECORD DATA======
	                float complexity = s.complexity();
	                float printability = s.printability();

	                TableRow newRow = dataOut.addRow();
	                newRow.setInt("iteration",iteration);
	                newRow.setString("genes",chrString);
	                newRow.setFloat("total_fitness", ((printability+complexity)/2));
	                newRow.setFloat("printability", printability);
	                newRow.setFloat("complexity", complexity);
	                
	                
	              //======RENDER=======
	                stroke(30, 150);
	                strokeWeight(0.5f);
	                Colony3D col = s.getColony3D();
	                Colony topC = col.getLayer(col.size() - 1);
	                for(int j = 0; j < col.size(); j+=3){
	                    Colony c = col.getLayer(j);
	                    displayColony(c, j);

	                }
	                saveFrame(folderPath+"imgs/run_"+prefix+"/iteration_"+String.format("%04d", iteration)+".png");

	                // reset background
	                background(220);
	                
	                println(prefix+" "+iteration);
	                // break;
	            }
	            
	            saveTable(dataOut, folderPath+"data/run_"+prefix+"_(new_printability).csv");
	            
	    	}
	    }
		
		
	}
//	
//	public void draw() {
//		background(220);
//		pushMatrix();
//		scale(1,-1,1);
//		rect(20,20,100,100);
//		popMatrix();
//	}
	
	void displayColony(Colony _c, int _h){
	    pushMatrix();
	    scale(1, -1, 1);
	    for(Organism o : _c.organisms){
	        for(Spring s : o.getSprings()){
	            line(s.sp.loc.x, s.sp.loc.y, _h, s.ep.loc.x, s.ep.loc.y, _h);
	        }

	    }
	    popMatrix();
	}
//	
	public static void main(String args[]) {
		String[] appletArgs = new String[] {"generative.components.examples.ImgOutput"};
		PApplet.main(appletArgs);
		
		
	}
}
