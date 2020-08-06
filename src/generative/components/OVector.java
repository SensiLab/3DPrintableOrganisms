package generative.components;

import java.io.Serializable;

public class OVector implements Serializable {
	public float x;
	public float y;
	public float z;
	
	/**
	 * Constructors
	 */
	
	
	public OVector(){
		this.x = 0f;
		this.y = 0f;
		this.z = 0f;
	}
	
	public OVector(float _x, float _y){
		this.x = _x;
		this.y = _y;
		this.z = 0f;
	}
	
	public OVector(float _x, float _y, float _z){
		this.x = _x;
		this.y = _y;
		this.z = _z;		
	}
	
//	==============METHODS================
	
	public OVector add(OVector other){
		this.x += other.x;
		this.y += other.y;
		this.z += other.z;
		
		return this;
	}
	
	public OVector sub(OVector other) {
		this.x-=other.x;
		this.y-=other.y;
		this.z-=other.z;
		
		return this;
	}
	
	public static OVector sub(OVector a, OVector b) {
		OVector target = new OVector(a.x-b.x, a.y-b.y, a.z-b.z);
		return target;
	}
	
	public OVector mult(float n) {
		this.x*=n;
		this.y*=n;
		this.z*=n;
		return this;
	}
	
	public OVector div(float n) {
		this.x/=n;
		this.y/=n;
		this.z/=n;
		return this;
	}
	
	public float magSq() {
		return (x*x + y*y + z*z);
	}
	
	public float mag(){
		return (float) Math.sqrt(x*x + y*y + z*z);
	}
	
	public OVector normalize() {
		    float m = mag();
		    if (m != 0 && m != 1) {
		      div(m);
		    }
		    return this;
		  }
	
	public OVector setMag(float n) {
		this.normalize();
		this.mult(n);
		return this;
	}
	
	
	public OVector limit(float max) {
	    if (magSq() > max*max) {
	      normalize();
	      mult(max);
	    }
	    return this;
	  }
	
	public float dist(OVector other) {
		float dx = this.x - other.x;
		float dy = this.y - other.y;
		float dz = this.z - other.z;
		
		return (float) Math.sqrt(dx*dx+dy*dy+dz*dz);
	}
	
	public OVector cross(OVector other) {
		float crossX = this.y * other.z - other.y * z;
	    float crossY = this.z * other.x - other.z * x;
	    float crossZ = this.x * other.y - other.x * y;
	    
	    OVector target = new OVector(crossX, crossY, crossZ);
	    
	    return target;
	}
	
	/**
	 * Calculate the dot product between two vectors
	 * @param other (in PVector format)
	 * @return dot product
	 */
	public float dot(OVector other) {
		return (this.x*other.x+this.y*other.y+this.z*other.z);
	}
	
	/**
	 * Calculate dot product between two vectors
	 * static function
	 * @param Vector 1 and Vector 2
	 * @return dot product between v1 and v2
	 */
	public static float dot(OVector v1, OVector v2) {
		return(v1.x*v2.x+v1.y*v2.y+v1.z*v2.z);
	}
	
	public static OVector getMidPoint(OVector v1, OVector v2) {
		float dx = v2.x - v1.x;
		float dy = v2.y - v1.y;
		float dz = v2.z - v1.z;
		
		return new OVector(v1.x + (dx/2), v1.y + (dy/2), v1.z + (dz/2));
	}
	
	public static OVector getNormal(OVector v1, OVector v2) {
		float dx = v2.x - v1.x;
		float dy = v2.y - v1.y;
		
		OVector normal = new OVector(dy, -dx);
		
		return normal.normalize();
	}
	
	public static OVector getNormal(float x1, float y1, float x2, float y2) {
		float dx = x2 - x1;
		float dy = y2 - y1;
		
		OVector normal = new OVector(dy, -dx);
		
		return normal.normalize();
	}
	
}
