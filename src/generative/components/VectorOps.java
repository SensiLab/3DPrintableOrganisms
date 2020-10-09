package generative.components;

import java.io.Serializable;
import processing.core.PVector;

public class VectorOps{

	/**
	 * Constructors
	 */
	
	
//	public VectorOps(){
//		this.x = 0f;
//		this.y = 0f;
//		this.z = 0f;
//	}
//	
//	public VectorOps(float _x, float _y){
//		this.x = _x;
//		this.y = _y;
//		this.z = 0f;
//	}
//	
//	public VectorOps(float _x, float _y, float _z){
//		this.x = _x;
//		this.y = _y;
//		this.z = _z;		
//	}
	
//	==============METHODS================
	
//	public PVector add(PVector other){
//		this.x += other.x;
//		this.y += other.y;
//		this.z += other.z;
//		
//		return this;
//	}
	
//	public PVector sub(PVector other) {
//		this.x-=other.x;
//		this.y-=other.y;
//		this.z-=other.z;
//		
//		return this;
//	}
//	
//	public static PVector sub(PVector a, PVector b) {
//		PVector target = new PVector(a.x-b.x, a.y-b.y, a.z-b.z);
//		return target;
//	}
//	
//	public PVector mult(float n) {
//		this.x*=n;
//		this.y*=n;
//		this.z*=n;
//		return this;
//	}
//	
//	public PVector div(float n) {
//		this.x/=n;
//		this.y/=n;
//		this.z/=n;
//		return this;
//	}
//	
//	public float magSq() {
//		return (x*x + y*y + z*z);
//	}
//	
//	public float mag(){
//		return (float) Math.sqrt(x*x + y*y + z*z);
//	}
	
//	public PVector normalize() {
//		    float m = mag();
//		    if (m != 0 && m != 1) {
//		      div(m);
//		    }
//		    return this;
//		  }
//	
//	public PVector setMag(float n) {
//		this.normalize();
//		this.mult(n);
//		return this;
//	}
//	
//	
//	public PVector limit(float max) {
//	    if (magSq() > max*max) {
//	      normalize();
//	      mult(max);
//	    }
//	    return this;
//	  }
	
//	public float dist(PVector other) {
//		float dx = this.x - other.x;
//		float dy = this.y - other.y;
//		float dz = this.z - other.z;
//		
//		return (float) Math.sqrt(dx*dx+dy*dy+dz*dz);
//	}
//	
//	public PVector rotate(float angle) {
//		float newX = (float)((this.x * Math.cos((double)angle)) - (this.y*Math.sin((double)angle)));
//		float newY = (float)((this.y * Math.cos((double)angle)) + (this.x*Math.sin((double)angle)));
//		return new PVector(x,y);
//	}
//	
//	public PVector cross(PVector other) {
//		float crossX = this.y * other.z - other.y * z;
//	    float crossY = this.z * other.x - other.z * x;
//	    float crossZ = this.x * other.y - other.x * y;
//	    
//	    PVector target = new PVector(crossX, crossY, crossZ);
//	    
//	    return target;
//	}
	
//	/**
//	 * Calculate the dot product between two vectors
//	 * @param other (in PVector format)
//	 * @return dot product
//	 */
//	public float dot(PVector other) {
//		return (this.x*other.x+this.y*other.y+this.z*other.z);
//	}
//	
//	/**
//	 * Calculate dot product between two vectors
//	 * static function
//	 * @param Vector 1 and Vector 2
//	 * @return dot product between v1 and v2
//	 */
//	public static float dot(PVector v1, PVector v2) {
//		return(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
//	}
	
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
	
//	public float heading() {
//		float angle = (float) Math.atan2(y, x);
//		return angle;
//	}
	
}
