package generative.components;

import java.io.Serializable;

public class FoodSource implements Serializable{
	int col;
	int row;

	int totalEnergy;
	float freq;

	public FoodSource(int _col, int _row, int _totalEnergy) {
		this.col = _col;
		this.row = _row;

	    this.totalEnergy = _totalEnergy;
	    
//	    System.out.println("New food source:");
//	    System.out.println("Loc: ["+this.col+", "+this.row+"]");
//	    System.out.println("Total energy: "+this.totalEnergy);
	}
	
	public float getTotalEnergy() {
		return this.totalEnergy;
	}
	
	public int getRow() {
		return this.row;
	}
	
	public int getCol() {
		return this.col;
	}

	public Boolean inOrganism(Organism o, int resolution) {
		float x3 = (this.col * resolution) + resolution/2;
		float y3 = (this.row * resolution) + resolution/2;

		float x4 = 999999f;
		float y4 = (this.row * resolution) + resolution/2;

		int intersections = 0;
		for (Spring s : o.springs) {
			float x1 = s.sp.loc.x;
			float y1 = s.sp.loc.y;
			float x2 = s.ep.loc.x;
			float y2 = s.ep.loc.y;
			float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

			float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
			float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;

			if (t <= 1f && t >= 0f && u <= 1f && u >= 0f) intersections++;
		}
		return intersections%2 == 1;
	}
}
