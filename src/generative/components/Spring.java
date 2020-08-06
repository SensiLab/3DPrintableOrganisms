package generative.components;

import java.util.ArrayList;

public class Spring{
	public Cell sp;
	public Cell ep;
	
	private float restLen;
	public float lenCoef = 0.1f;
	private float sprCoef;
	public float minRestLen;
	public float maxLen;
	public float newSpringLenMult;
	
	float splitThreshold = 0.3f;
	
	public Spring(Cell _sp, Cell _ep) {
	    this.sp = _sp;
	    this.ep = _ep;

	    this.minRestLen = this.sp.repRad + this.ep.repRad;

	    this.restLen = 20;
	    this.sprCoef = 0.5f;

	    this.newSpringLenMult = 0.6f;	    
	}
	
	public Spring(Cell _sp, Cell _ep, float _rl, float _sc, float _nl) {
		this.sp = _sp;
		this.ep = _ep;

		this.restLen = _rl;
		this.sprCoef = _sc;
		this.newSpringLenMult = _nl;

		this.minRestLen = this.sp.repRad + this.ep.repRad;
	}
	
	// ================METHODS===========
	/**
	 * Getter method to retrieve the rest length of springs
	 * @return the rest length of the spring
	 */
	
	public float getRestLen() {
		return this.restLen;
	}
	
	/**
	 * Getter method to retrieve split threshold
	 * 
	 */
	public float getSplitThreshold() {
		return this.splitThreshold;
	}
	
	/**
	 * Setter method to set splitThreshold
	 */
	public void setSplitThreshold(float t) {
		this.splitThreshold = t;
	}

	/**
	 * Getter method to retrieve the spring coefficient of a
	 * spring (i.e. how hard or soft the spring is)
	 * @return the spring coefficient of the current spring
	 */
	public float getSpringCoef() {
		return this.sprCoef;
	}
	
	public void setSpringCoef(float sc) {
		this.sprCoef = sc;
	}
	
	/**
	 * Measures the actual length of a spring
	 * @return spring current length
	 */
	public float getLen() {
		return sp.loc.dist(ep.loc);
	}
	
	/**
	 * Setter for spring rest length
	 * @param l
	 */
	public void setRestLen(float l) {
		this.restLen = l;
	}

	/**
	 * Measures the energy that corresponds to a spring by
	 * taking half of the energies of the cells connected by it
	 * @return 
	 */
	public float getEnergy() {
		return (this.sp.energy/2) + (this.ep.energy/2);
	}
	
	/**
	 * Finds the mid point between the two ends of a spring
	 * @return midpoint in PVector format
	 */
	public OVector getMidpoint() {
		float dx = this.ep.loc.x - this.sp.loc.x;
		float dy = this.ep.loc.y - this.sp.loc.y;
		OVector mp = new OVector(sp.loc.x + (dx/2), sp.loc.y + (dy/2));
		return mp;
	}
	
	/**
	 * calculates the normal of a spring (perpendicular
	 * to the spring)
	 * @return normalised PVector perpendicular to spring
	 */
	public OVector getNormal() {
		float dx = this.ep.loc.x - this.sp.loc.x;
		float dy = this.ep.loc.y - this.sp.loc.y;
		
		OVector normal = new OVector(dy, -dx);
		
		return normal.normalize();
	}
	
	/**
	 * Calculates the dot product between the normal of a spring and the
	 * vector between the midpoint of the spring and a light source.
	 * @param ls (Location vector of light source)
	 * @return dot product
	 */
	public float getDot(OVector ls) {
		OVector mp = this.getMidpoint();
		OVector pointer = OVector.sub(ls, mp);
		pointer.normalize();
		
		OVector normal = this.getNormal();
		
		return OVector.dot(normal, pointer);
	}
	
	/**
	 * calculates the intersection point between two springs
	 * @param other (Spring)
	 * @return intersection point PVector
	 */
	public OVector getIntersectionPoint(Spring other) {
	    float x1 = other.sp.loc.x;
	    float y1 = other.sp.loc.y;
	    float x2 = other.ep.loc.x;
	    float y2 = other.ep.loc.y;

	    float x3 = this.sp.loc.x;
	    float y3 = this.sp.loc.y;
	    float x4 = this.ep.loc.x;
	    float y4 = this.ep.loc.y;

	    float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
	    if (den == 0) return null;

	    float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4))/den;
	    float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3))/den;
	    
	    if (t > 1f || t < 0f || u > 1f || u < 0f) return null;
	    
	    float x = x1 + (t*(x2-x1));
	    float y = y1 + (t*(y2-y1));

	    return new OVector(x,y);
	}

	public void update() {
		// set rest length to the sum of the energies contained in connected cells
	    float rl = (this.sp.energy + this.ep.energy);
	    this.restLen = rl;
	    
	    if (this.getLen() > this.restLen) {
			// contract spring
	    	// pull in ep
			OVector spToEp = OVector.sub(this.sp.loc, this.ep.loc);
			this.ep.applyForce(spToEp.mult(this.sprCoef));
			// pull in sp
			OVector epToSp = OVector.sub(this.ep.loc, this.sp.loc);
			this.sp.applyForce(epToSp.mult(this.sprCoef));
	    } else {
	    	// expand spring
	    	// push ep
			OVector spToEp = OVector.sub(this.ep.loc, this.sp.loc);
			this.ep.applyForce(spToEp.mult(this.sprCoef));
			
			//push sp
			OVector epToSp = OVector.sub(this.sp.loc, this.ep.loc);
			this.sp.applyForce(epToSp.mult(this.sprCoef));
		}
	}
}