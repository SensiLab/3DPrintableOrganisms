package generative.components;

import processing.core.PApplet;

public class Visualisation2D extends PApplet {
	
	Colony3D col;
	
	Visualisation2D(Colony3D _col){
		this.col = _col;
	}
	
	public void settings() {
		
	}
	
//	@Override
//	public void setup() {
//		col.setPApplet(this);
//		background(220);
//		
//	}
	
	public void draw() {
		background(220);
        stroke(20,50);
        strokeWeight(1f);
        for(Colony c : col.layers) {
        	displayColony(c);
        }
	}
	
	void displayColony(Colony _c){
	    for(Organism o : _c.organisms){
	        for(Spring s : o.getSprings()){
	            line(s.sp.loc.x,s.sp.loc.y,s.ep.loc.x,s.ep.loc.y);
	        }
	        
	    }
	}
	
	

}
