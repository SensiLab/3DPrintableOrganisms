package generative.components.examples;

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

	long[] envSeeds = {-1632478226818474777L, 8558593120860338761L, 9063855533770827947L, 8141159607898132170L, -2228090770267545711L,
					5760679812470627183L, 8122425018693031776L, 7233387513096301004L, 8063979647430030046L, 1209585713320896596L, 
					724224084362998227L, -5071659569307755225L, 5823929498051207574L, -551995978363359132L, -5192310287763653706L,
					5108100186193579751L, 4676335984590537776L, 7191218674032234786L, 7152999658342319501L, 7062340219330815167L,
					4871260159946803014L, -5288381650226216036L, -1164060861243111837L, -7174953927911303402L, -9094118736910298012L,
					-3497492774682784886L, -1682527197870981019L, 476367449420437529L, -2508839713788940676L, -6172147982322942294L, 
					-8486584721678693775L, 7170492297316362611L, -7692137076181517887L, 3861397882188512846L, -3110163977026367975L,
					-8634546392007882087L, 7557350025982015301L, 4025513581029972815L, -1483865547251428910L, -309284057398871553L,
					-186292376937715704L, -4463713339575002910L, 7878392066209648151L, -4195456362027586512L, 1429266236894580111L,
					3045732330329506990L, 4597571254417314837L, 5927936720052420491L, -1593416348006868653L, 3027521254663839597L};
	
	
	public void settings() {
		size(envWidth,envHeight,P3D);
		
	}
	
	public void setup() {
		camera = new PeasyCam(this, 300, -300, 0, 675);
	    camera.rotateZ(PI/4);
	    camera.rotateX(-PI/4);
	    camera.pan(0,-225);
		
		background(220);
		
		String folderPath = "../out/evo_output/IEEE_Refine_013/";
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
	            Long rSeed = envSeeds[prefixInt];
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
	                Environment e = new Environment(envWidth, envHeight, res, drag, fss, fsSizeMin, fsSizeMax, foodGrowthRate, foodDecayRate, rSeed);
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
