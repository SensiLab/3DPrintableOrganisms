package generative.components;

public class MathLib {
	
	  /*
	  max is the max value of the curve, mp is the value of x at the curve's mid
	  point and k is the steepness of the curve
	  */
	  public static float sigmoid(float x, float max, float mp, float k) {
	    return  max/ (float)(1 + Math.exp(-k*(x - mp)));
	  }
	  
	  /*
	  method to map one value from one range to another
	  */
	  public static float map(float value, float valueRangeMin, float valueRangeMax, float newRangeMin, float newRangeMax) {
		  float r = ((value - valueRangeMin) / (valueRangeMax - valueRangeMin)) * (newRangeMax - newRangeMin) + newRangeMin;
		  return r;
	  }
	  
	  public static double mapDouble(double value, float valueRangeMin, float valueRangeMax, float newRangeMin, float newRangeMax) {
		  double r = ((value - valueRangeMin) / (valueRangeMax - valueRangeMin)) * (newRangeMax - newRangeMin) + newRangeMin;
		  return r;
	  }
	  
	  
	//returns the min distance between a point pt and a line ln
	  public static float distPointLine(OVector pt, OVector[] ln) {
		  float x = pt.x;
		  float y = pt.y;

		  float x1 = ln[0].x;
		  float y1 = ln[0].y;

		  float x2 = ln[1].x;
		  float y2 = ln[1].y;

		  float d = Math.abs(((y2 - y1)*x)-((x2-x1)*y) + (x2*y1) - (y2*x1))/ (float) Math.sqrt(Math.pow(y2 - y1, 2) + Math.pow(x2 - x1,2));
		  return d;
	  }
}
