package generative.components;

import java.util.ArrayList;
import java.util.Collections;

import processing.core.PVector;
import processing.core.*;

public class CommonFunctions {
	
	/**
	 * Calculates the distance between a point and a line on the xy plane
	 * @param pt
	 * @param ln
	 * @return distance as a float
	 */
	
	public static float distPointLine(PVector pt, PVector[] ln) {
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
	
//	public static 
	
	public static int[] closestEdge(PVector pt, Colony layer){
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
	
	public static double std(double[] data) {
		double sum = 0;
		double stDev = 0;
		int len = data.length;
		
		for(double n : data) {
			sum+=n;
		}
		
		double mean = sum/len;
		
		for(double n : data) {
			stDev += Math.pow(n - mean,  2);
		}
		
		return Math.sqrt(stDev)/len;
	}
	
	public static void displayColony(PApplet _pa, Colony _c, float layerZ){
		// Scale colony based on drawing size
		// Horizontal connections
	    for(Organism o : _c.organisms){
	        for(Spring s : o.getSprings()){
	        	float spX = s.sp.loc.x * _pa.width/_c.e.width;
	        	float spY = s.sp.loc.y * _pa.height/_c.e.height;
	        	float spZ = layerZ;
	        	float epX = s.ep.loc.x * _pa.width/_c.e.width;
	        	float epY = s.ep.loc.y * _pa.height/_c.e.height;
	        	float epZ = layerZ;
	            _pa.line(spX, spY, spZ, epX, epY, epZ);
//	            ellipse(spX, spY, width * .01f, height * .01f);
	        }
	        
//	        for(Cell c : o.getCells()) {
//	        	drawVector(_pa, c.vel, c.loc.x, c.loc.y, 20);
//	        }
	        
	        
	    }
	    
	   
	    
	}
	
	public static void displayColony3D(PApplet _pa, Colony3D _c, int _interval){
		// Scale colony based on drawing size
		
	    for(int i = 0; i < _c.size();i+=_interval) {
	    	Colony _c2d = _c.getLayer(i);
	    	float layerZ = _c2d.getHeight();
	    	for(Organism o : _c2d.organisms){
		        for(Spring s : o.getSprings()){
		        	float spX = s.sp.loc.x * _pa.width/_c2d.e.width;
		        	float spY = s.sp.loc.y * _pa.height/_c2d.e.height;
		        	float spZ = layerZ;
		        	float epX = s.ep.loc.x * _pa.width/_c2d.e.width;
		        	float epY = s.ep.loc.y * _pa.height/_c2d.e.height;
		        	float epZ = layerZ;
		            _pa.line(spX, spY, spZ, epX, epY, epZ);
//		            ellipse(spX, spY, width * .01f, height * .01f);
		        }
		        
		    }
	    }
		
	}
	
	public static void displayColony3DMesh(PApplet _pa, Colony3D _c) {
		// draw horizontal lines
		for(int i = 0; i < _c.size(); i++) {
			Colony _c2d = _c.getLayer(i);
			float layerZ = _c2d.getHeight();
			for(Organism o : _c2d.organisms){
		        for(Spring s : o.getSprings()){
		        	float spX = s.sp.loc.x * _pa.width/_c2d.e.width;
		        	float spY = s.sp.loc.y * _pa.height/_c2d.e.height;
		        	float spZ = layerZ;
		        	float epX = s.ep.loc.x * _pa.width/_c2d.e.width;
		        	float epY = s.ep.loc.y * _pa.height/_c2d.e.height;
		        	float epZ = layerZ;
		            _pa.line(spX, spY, spZ, epX, epY, epZ);
//		            ellipse(spX, spY, width * .01f, height * .01f);
		        }
			}
			if(i < _c.size()-1) {
				Colony next = _c.getLayer(i+1);
				float nextZ = next.getHeight();
				if(next.organisms.size() == _c2d.organisms.size()) {
					for(int j = 0; j < _c2d.organisms.size(); j++) {
						Organism current_o = _c2d.organisms.get(j);
						Organism next_o = next.organisms.get(j);
						if(current_o.getCells().size() == next_o.getCells().size()) {
							for(int k = 0; k < current_o.getCells().size();k++) {
								Cell bottom = current_o.getCells().get(k);
								Cell top = next_o.getCells().get(k);
								float bottomX = bottom.loc.x * _pa.width/_c2d.e.width;
					        	float bottomY = bottom.loc.y * _pa.height/_c2d.e.height;
					        	float bottomZ = layerZ;
					        	float topX = top.loc.x * _pa.width/_c2d.e.width;
					        	float topY = top.loc.y * _pa.height/_c2d.e.height;
					        	float topZ = nextZ;
								_pa.line(bottomX,  bottomY,  bottomZ,  topX, topY, topZ);
							}
						}
					}
				}
			}
			
		}
	}
	
	public static void displayEnvironment(PApplet _pa, Environment e) {
		float gridSize = _pa.width/e.cols;
		//draw grid
		_pa.stroke(100f);
		_pa.strokeWeight(0.5f);
		for(int i = 1; i < e.cols; i++) {
			float x = i * gridSize;
			_pa.line(x,0,x,_pa.height);
			_pa.line(0,x,_pa.width,x);
		}
		//draw food
		for(int i = 0; i < e.cols; i++) {
			for(int j = 0; j < e.rows; j++) {
				int[] tile = {i,j};
//				if(e.lookupNutrients(tile) > 0) {
//					
//					_pa.textSize(9);
//					String nutrients = String.format("%.2f", e.lookupNutrients(tile));
//					float r = e.lookupNutrients(tile) * gridSize;
//					float x = tile[0] * gridSize + gridSize/2;
//					float y = tile[1] * gridSize + gridSize/2;
//					_pa.pushMatrix();
//					_pa.translate(x, y);
//					_pa.fill(50,20);
////					_pa.noFill();
//					_pa.rectMode(2);
//					_pa.rect(0, 0, gridSize/2, gridSize/2);
//					_pa.fill(0);
//					_pa.textAlign(2,2);
//					_pa.text(nutrients,0,0);
//					_pa.popMatrix();
//				}
				
				
				if(e.lookupNutrients(tile) > 0) {
					float r = e.lookupNutrients(tile) * gridSize;
					float x = tile[0] * gridSize + gridSize/2;
					float y = tile[1] * gridSize + gridSize/2;
					_pa.noStroke();
					_pa.fill(255,200);
					_pa.ellipse(x,y,r,r);
				}
			}
		}
	}
	
	public static void drawVector(PApplet _pa, PVector v, float x, float y, float scale) {
		  _pa.pushMatrix();
		  float arrowsize = 4;
		  _pa.translate(x, y);
		  _pa.stroke(0, 100);
		  _pa.rotate(v.heading());
		  float len = v.mag() * scale;
		  _pa.line(0, 0, len, 0);
		  _pa.line(len, 0, len-arrowsize, arrowsize/2);
		  _pa.line(len, 0, len-arrowsize, -arrowsize/2);
		  _pa.popMatrix();
	}
	
	public static double hammingDistance(int[] a, int[] b) {
		if(a.length != b.length) {
			System.out.println("Hamming Distance requires both arrays to be of the same length");
			return -9999;
		}
		int[] c = new int[a.length];
		for(int i = 0; i < a.length; i++) {
			boolean first = a[i] != 0;
			boolean second = b[i] != 0;
			if(first ^ second) {
				c[i] = 1;
			}else {
				c[i] = 0;
			}
		}
		int totalDist = 0;
		for(int i : c) {
			totalDist+=i;
		}
		return (double)totalDist/c.length;
		
	}
	
	public static float dispersionCoefficient(ArrayList<Float> a){
//		ArrayList<Float> a = this.getRealAngles();
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
	
	// shift cells
	public static ArrayList<Cell> shiftCells(ArrayList<Cell> c, int s){
		ArrayList<Cell> shifted = new ArrayList<Cell>();
		for(int i = s; i < c.size(); i++) {
			shifted.add(c.get(i));
		}
		for(int i = 0; i < s; i++) {
			shifted.add(c.get(i));
			
		}
		return shifted;
		
	}//shift cells
	
	// bottomMost cell
	public static int bottomMostCell(ArrayList<Cell> _cells) {
		int winner_index = -1;
		PVector winner = new PVector();
		for (int i=0; i < _cells.size(); i++) {
			PVector p = _cells.get(i).loc;
			if (p.y > winner.y) {
				winner = p;
				winner_index = i;
			}else if(p.y == winner.y){
				if(p.x > winner.x){
					winner = p;
					winner_index = i;
				}
			}else{
				continue;
			}
		}
		return winner_index;
	}//bottomMost cell
	
	// polar sort cells
	public static ArrayList<Cell> polarSort(ArrayList<Cell> cells) {
		//make a copy of pts array
		ArrayList<Cell> cellsCopy = new ArrayList<Cell>();
		for (Cell c : cells) cellsCopy.add(new Cell(c));
		// Remove duplicates
		for (int i =0; i < cellsCopy.size()-1;i++) {
			float ix = cellsCopy.get(i).loc.x;
			float iy = cellsCopy.get(i).loc.y;
			for (int j = i+1; j < cellsCopy.size(); j++) {
				float jx = cellsCopy.get(j).loc.x;
				float jy = cellsCopy.get(j).loc.y;
				if(ix == jx && iy == jy) {
					cellsCopy.remove(j);
				}
			}
		}
		//initialise sorted array
		ArrayList<Cell> sortedCells = new ArrayList<Cell>();

		//get bottomMost
		int bottom_most_index = bottomMostCell(cellsCopy);
		

		// add bottomMost to sorted and remove it from cellsCopy		
		sortedCells.add(cellsCopy.get(bottom_most_index));
		cellsCopy.remove(bottom_most_index);
		
		//create bottom most vector
		PVector bm = sortedCells.get(0).loc;

		
		double prevMinCotan = -1;
		while (cellsCopy.size() > 0) {
			double minCotan = 999999999;
			PVector winner = new PVector();
			int winner_index = -1;
			for(int i = 0; i < cellsCopy.size(); i++) {
				PVector p = cellsCopy.get(i).loc;
				double c = (p.x - bm.x)/(p.y - bm.y);
				if(c < minCotan){
					minCotan = c;
					winner = p;
				}else if(c == minCotan){
					minCotan = c;
					if(bm.dist(p) <= bm.dist(winner)) winner = p;
				}else{
					continue; 
				}
			}

			PVector lastLoc = sortedCells.get(sortedCells.size() - 1).loc;
			if(winner.x != lastLoc.x && winner.y != lastLoc.y) {
				sortedCells.add(new Cell(cellsCopy.get(winner_index)));
				cellsCopy.remove(winner_index);
			}else {
				cellsCopy.remove(winner_index);
			}
			prevMinCotan = minCotan;
		}
		return sortedCells;
	}
	
	public static float[][] transpose(float[][] array) {
	    // empty or unset array, nothing do to here
	    if (array == null || array.length == 0)
	        return array;

	    int width = array.length;
	    int height = array[0].length;

	    float[][] array_new = new float[height][width];

	    for (int x = 0; x < width; x++) {
	        for (int y = 0; y < height; y++) {
	            array_new[y][x] = array[x][y];
	        }
	    }
	    return array_new;
	}
	
	public static void print2DArray(float[][] array) {
		float[][] transposed = transpose(array);
		for(float[] row : transposed) {
			printRow(row);
		}
	}
	 public static void printRow(float[] row) {
        for (float i : row) {
            System.out.print(String.format("%.2f", i));
            System.out.print("\t");
        }
        System.out.println();
    }
}
