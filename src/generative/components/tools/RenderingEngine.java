package generative.components.tools;

import generative.components.*;
import processing.core.*;
import processing.core.PVector;
import processing.data.*;
import java.io.*;
import java.util.Arrays;
import java.util.Scanner;

public class RenderingEngine {
	public static void main(String[] args) throws IOException {
		String folderPath = "../out/evo_output/Test_Full_Fitness_002";
		File folder = new File(folderPath);
		
		for (final File fileEntry : folder.listFiles()) {
			String fileName = fileEntry.getName();
			String run = fileName.substring(0, 3);
	        
			if (fileName.contains("Genes")) {
	        	System.out.println(run);
	        	Table t = new Table(fileEntry);
	        	String envData = "";
	        	
	        	for (final File fitnessFile : folder.listFiles()) {
	        		String fitFileName = fitnessFile.getName();
	        		String runF = fitFileName.substring(0,3);
	        		if(runF.equals(run) && !fitFileName.equals(fileName)) {
	        			
	        			Table ft = new Table(fitnessFile);
	        			TableRow envDataRow = ft.getRow(1);
	        			for(int i = 2; i < envDataRow.getColumnCount(); i++) {
	        				envData+=envDataRow.getString(i);
	        				
	        			}
//	        			ft.print();
	        			System.out.println(envDataRow.getColumnCount());
	        			System.out.println(envData);
	        			
	        		}
	        	}
	        	
	        	
//	        	System.out.println(t.getColumn(0));
//	        	t.print();
	        	TableRow tr = t.getRow(0);
	        	String[] envSeedCell = tr.getString(1).split("=");
	        	Long envSeed = Long.parseLong(envSeedCell[1]);
	        	System.out.println(envSeed);
//	        	int[] colTypes = tr.getColumnTypes();
//	        	System.out.println(Arrays.toString(colTypes));
//	            System.out.println(fileEntry.getName());
	        	
	        	Table dataOut = new Table();
	            dataOut.addColumn("iteration");
	            dataOut.addColumn("genes");
	            dataOut.addColumn("total_fitness");
	            dataOut.addColumn("printability");
	            dataOut.addColumn("complexity");
	            dataOut.addColumn("convexity");
	            dataOut.addColumn("angle_dispersion");
	            
	            for(int i = 2; i < 152 ; i++){
	            	TableRow r = t.getRow(i);
	                int iteration = r.getInt(0);
	                float fitness = r.getFloat(4);
	                String chrString = r.getString(5).replace("(", "").replace(")","").replace(";",",");
	                String[] chrStringArray = chrString.split(",");

	                double[] genes = new double[chrStringArray.length];
	                for(int k = 0; k < chrStringArray.length; k++){
	                	genes[k] = Double.parseDouble(chrStringArray[k]);
	                }

//	                Chromosome chr = new Chromosome(genes);
//	                Environment e = new Environment(envWidth, envHeight, res, drag, fss, fsSizeMin, fsSizeMax, foodGrowthRate, foodDecayRate, rSeed);
//	                Simulation s = new Simulation(chr, e, 500, 20);
//	                s.generate();
	            }
	        	
	        	break;
	        }
	    }
		
	}
}
