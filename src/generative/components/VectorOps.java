package generative.components;


import processing.core.PVector;

public class VectorOps{

	public static PVector getMidPoint(PVector v1, PVector v2) {
		float dx = v2.x - v1.x;
		float dy = v2.y - v1.y;
		float dz = v2.z - v1.z;
		
		return new PVector(v1.x + (dx/2), v1.y + (dy/2), v1.z + (dz/2));
	}
	
	public static PVector getNormal(PVector v1, PVector v2) {
		float dx = v2.x - v1.x;
		float dy = v2.y - v1.y;
		
		PVector normal = new PVector(dy, -dx);
		
		return normal.normalize();
	}
	
	public static PVector getNormal(float x1, float y1, float x2, float y2) {
		float dx = x2 - x1;
		float dy = y2 - y1;
		
		PVector normal = new PVector(dy, -dx);
		
		return normal.normalize();
	}
	
}
