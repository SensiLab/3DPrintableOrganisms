package generative.components;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import processing.core.PVector;

public class Colony3D{

	
	ArrayList<Colony> layers;
	
	public Colony3D(){
		this.layers = new ArrayList<Colony>();
	}
	
	public void addLayer(Colony layer) {
		this.layers.add(layer);
	}
	
	public void clear() {
		this.layers.clear();
	}
	
	public int size() {
		return this.layers.size();
	}
	
	public Colony getLayer(int i) {
		return this.layers.get(i);
	}
	
	public ArrayList<Colony> getLayers(){
		return this.layers;
	}
	
	public Colony3D scale(float xyFactor, float zFactor, Environment e) {
		Colony3D scaledColony = new Colony3D();
		
		for(Colony layer : this.layers) {
			Colony scaledLayer = new Colony();
			for(Organism org : layer.organisms) {
				Organism scaledOrg = new Organism();
				for (Cell c : org.getCells()) {
					float x = c.loc.x * xyFactor / e.width;
					float y = c.loc.y * xyFactor / e.height;
					float z = c.loc.z * zFactor;
					PVector newLoc = new PVector(x,y,z);
					Cell scaledC = new Cell(newLoc, layer.getChromosome());
					scaledOrg.addCell(scaledC);
				}
				scaledOrg.connectCells();
				scaledLayer.organisms.add(scaledOrg);
			}
			scaledColony.addLayer(scaledLayer);
		}
		return scaledColony;
	}
	
	public void saveColony(String filePath, String fileName) {
		
		// Create file
		File output = new File(filePath+fileName+".txt");
		try {
			FileWriter outputWriter = new FileWriter(filePath+fileName+".txt");
			outputWriter.write("Total Layers: "+this.layers.size()+"\n");
			
			// check every layer
			for(Colony layer : this.layers) {
				// check every organism
				for(Organism o : layer.organisms) {
//					outputWriter.write("");
					for(Cell c : o.cells) {
						float x = c.loc.x;
						float y = c.loc.y;
						float z = c.loc.z;
						outputWriter.write("("+x+","+y+","+z+")");
						if(o.cells.indexOf(c) < o.cells.size() - 1) {
							outputWriter.write(";");
						}
					}
					
				}
				outputWriter.write("\n");
			}
//			outputWriter.write("All layers should follow here!!\n");
			outputWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		// make every organism a string
		ArrayList<String> allOrgs = new ArrayList<String>();
	    //for every layer
	    for(Colony layer : this.layers){
	      //for every organism
	      for(Organism org : layer.organisms){
	        ArrayList<String> organism = new ArrayList<String>();
	        StringBuilder sb = new StringBuilder();
	        //add the start point for every line of org
	        for(Cell c : org.cells){
	          PVector sp = c.loc;
	          String pt = "("+Float.toString(sp.x)+","+Float.toString(sp.y)+","+Float.toString(sp.z)+")";
	          sb.append(pt);
	          if(org.cells.indexOf(c) != org.cells.size()-1) sb.append(";");
	        }
	        allOrgs.add(sb.toString());
	        // println(sb.toString());
	      }
	    }
	    String[] stringsOut = new String[allOrgs.size()];
	    for(int i = 0; i < allOrgs.size(); i++){
	      stringsOut[i] = allOrgs.get(i);
	    }
	    		
	}
	
}
