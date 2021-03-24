package generative.components;

import java.util.ArrayList;
import java.util.Random;
import processing.core.PVector;

public class Environment{

	long randSeed;
	Random generator;
//	int noiseSeed;
	
	long randSeedOut;
//	int noiseSeedOut;
	
	long maxSeed = 2147483647;
	
//	float[][] dragField;
	float drag;
	float[][] nutrientField;
	int cols, rows;
	int width, height;
	int resolution; //the size of the tiles
	
//	Colony colony;
	
	float xoffIncrements = .13f;
	float yoffIncrements = .13f;
	
	float foodGrowthP = 0.0001f;
	
	ArrayList<FoodSource> foodSources;
	ArrayList<FoodSource> initFS;
	int numOfFoodSources = 5;
	int foodSourceSizeMin;
	int foodSourceSizeMax;
//	int totEnergyOfFoodSources = 50;
	float foodGrowthRate = 10;
	float foodDecayRate = 0.05f;
	
	int timestep;
//
	public Environment(int _width, int _height, int res, float _drag, int _numOfFs, int fsSizeMin, int fsSizeMax, float _foodGrowthRate, long _randSeed) {
		
		this.randSeed = _randSeed;

		this.width = _width;
		this.height = _height;
		this.resolution = res;
		this.cols = this.width/this.resolution;
		this.rows = this.height/this.resolution;
		
		this.drag = _drag;

		this.numOfFoodSources = _numOfFs;
		this.foodSourceSizeMin = fsSizeMin;
		this.foodSourceSizeMax = fsSizeMax;
		this.foodGrowthRate = _foodGrowthRate;

//		dragField = new float[this.cols][this.rows];
		this.nutrientField = new float[this.cols][this.rows];

		this.initFS = new ArrayList<FoodSource>();
		this.foodSources = new ArrayList<FoodSource>();

		this.init();
	}
	
public Environment(int _width, int _height, int res, float _drag, int _numOfFs, int fsSizeMin, int fsSizeMax, float _foodGrowthRate, float _foodDecayRate, long _randSeed) {
		
		this.randSeed = _randSeed;

		this.width = _width;
		this.height = _height;
		this.resolution = res;
		this.cols = this.width/this.resolution;
		this.rows = this.height/this.resolution;
		
		this.drag = _drag;

		this.numOfFoodSources = _numOfFs;
		this.foodSourceSizeMin = fsSizeMin;
		this.foodSourceSizeMax = fsSizeMax;
		this.foodGrowthRate = _foodGrowthRate;
		this.foodDecayRate = _foodDecayRate;

//		dragField = new float[this.cols][this.rows];
		this.nutrientField = new float[this.cols][this.rows];

		this.initFS = new ArrayList<FoodSource>();
		this.foodSources = new ArrayList<FoodSource>();

		this.init();
	}
		
	//========COPY CONSTRUCTOR========
	public Environment(Environment _e){
		this.randSeed = _e.randSeedOut;
		this.width = _e.width;
		this.height = _e.height;
		//    this.noiseSeed = _e.noiseSeed;

		this.resolution = _e.resolution;
		this.cols = this.width/this.resolution;
		this.rows = this.height/this.resolution;
    
		this.drag = _e.drag;

		this.numOfFoodSources = _e.numOfFoodSources;
		this.foodSourceSizeMin = _e.foodSourceSizeMin;
		this.foodSourceSizeMax = _e.foodSourceSizeMax;
		this.foodGrowthRate = _e.foodGrowthRate;


		//    dragField = new float[this.cols][this.rows];
		this.nutrientField = new float[this.cols][this.rows];

		this.initFS = new ArrayList<FoodSource>();
		this.foodSources = new ArrayList<FoodSource>();

		this.init();

	}

  
	// METHODS
	public int getWidth() {
		return this.width;
	}
	
	public int getHeight() {
		return this.height;
	}
	
	public int getRows() {
		return this.rows;
	}
	
	public int getCols() {
		return this.cols;
	}
	
	public float getEnergy(int _col, int _row) {
		return this.nutrientField[_col][_row];
	}
	
	public int getResolution() {
		return this.resolution;
	}
	
	public FoodSource[] getFoodSources() {
		FoodSource[] fs = this.foodSources.toArray(new FoodSource[this.foodSources.size()]);
		return fs;
	}
	
	public long getRandomSeed() {
		return this.randSeedOut;
	}
	
	public void setRandomSeed(long _rs) {
		this.randSeed = _rs;
		this.init();
	}
	
	public String getParams() {
		String size = "EnvSize=("+this.width+";"+this.height+")"+"\n";
		String res = "EnvRes="+this.resolution+"\n";
		String d = "EnvDrag="+this.drag+"\n";
		String fss = "EnvFoodSources="+this.numOfFoodSources+"\n";
		String fsSizeMin = "EnvFoodSourcesMin="+this.foodSourceSizeMin+"\n";
		String fsSizeMax = "EnvFoodSourcesMax="+this.foodSourceSizeMax+"\n";
		String fsGR = "EnvFoodSourceGR="+this.foodGrowthRate+"\n";
		String fsDR = "EnvFoodSourcesDR="+this.foodDecayRate+"\n";
		return size+res+d+fss+fsSizeMin+fsSizeMax+fsGR+fsDR;
	}
	
	void init() {
		if(this.width != this.height){
			throw new ArithmeticException("width and height MUST be the same");
		}
		this.generator = new Random();
		if(this.randSeed == -1){
			long seed = this.generator.nextLong();
			this.randSeedOut = seed;
			long r = this.randSeedOut;
      
		}else{
			this.randSeedOut = this.randSeed;
		}
		this.generator.setSeed(this.randSeedOut);
	
    //define position of food sources.
		if (this.initFS.size() == 0) {
			for (int i = 0; i< this.numOfFoodSources; i++) {
				FoodSource fs = new FoodSource(this.generator.nextInt(cols), this.generator.nextInt(rows), this.generator.nextInt(this.foodSourceSizeMax - this.foodSourceSizeMin) + this.foodSourceSizeMin);
				this.initFS.add(fs);
			}
		}
		this.foodSources = initFS;
		this.timestep = 0;
	}//init
	
	public void setfoodDecayRate(float rate) {
		this.foodDecayRate = rate;
	}

	public void updateNutrientField(Colony c){
		// initialise field for next time step
		float[][] newNutrientField = new float[this.cols][this.rows];
		for (int i = 0; i < this.cols; i++) {
			for (int j = 0; j < this.rows; j++) {
				newNutrientField[i][j] = this.nutrientField[i][j];
			}
		}
		
		//add food from food sources
		for (int i = 0; i < this.foodSources.size(); i++) {
			FoodSource fs = this.foodSources.get(i);
			int col = fs.col;
			int row = fs.row;

			if (fs.totalEnergy > 0) {
				newNutrientField[col][row] += 1f/this.foodGrowthRate;
				fs.totalEnergy-=(1f/this.foodGrowthRate);
			}
			
			if (fs.totalEnergy <=0) {
				// Remove the empty food source
				this.foodSources.remove(i);
				
				// Create a new food source
				FoodSource newFs = new FoodSource(this.generator.nextInt(cols), this.generator.nextInt(rows), this.generator.nextInt(this.foodSourceSizeMax - this.foodSourceSizeMin) + this.foodSourceSizeMin);

				this.foodSources.add(newFs);
			}
		}

		// update
		for (int i = 0; i< this.cols; i++) {
			for (int j = 0; j< this.rows; j++) {
				// 
				if (this.nutrientField[i][j] > 1f) {
					float energy = this.nutrientField[i][j] - 1;
					newNutrientField[i][j] = 1f;
					//	diffuse to neighbours
					ArrayList<int[]> neighbours = new ArrayList<int[]>();
					for (int x = -1; x<=1; x++) {
						for (int y = -1; y <= 1; y++) {
							int[] n = new int[2];
							n[0] = i + x;
							n[1] = j + y;
							if ((n[0] >= 0 && n[0] < this.cols) && x != 0 && y == 0) neighbours.add(n);
							if ((n[1] >= 0 && n[1] < this.rows) && y != 0 && x == 0) neighbours.add(n);
						}
					}
					for (int[] n : neighbours) {
						newNutrientField[n[0]][n[1]] += (energy / (float) neighbours.size());
					}
//  }

					if (this.nutrientField[i][j] > 0) {
						newNutrientField[i][j] -= this.foodDecayRate;
					} else {
						newNutrientField[i][j] = 0;
					}
				}
			}
		}
		
		for(int i = 0; i < this.cols; i++) {
			for(int j = 0; j < this.rows; j++) {
				this.nutrientField[i][j] = newNutrientField[i][j];
			}
		}
		
		this.timestep++;
	}

  public float lookupNutrients(PVector loc) {
    int col = (int) Math.min(Math.max(loc.x/resolution, 0), this.cols-1);
    int row = (int) Math.min(Math.max(loc.y/resolution, 0), this.rows-1);
    return this.nutrientField[col][row];
  }

  public float lookupNutrients(int[] tile) {
    int col = tile[0];
    int row = tile[1];
    return this.nutrientField[col][row];
  }
  
  public void depleteNutrients(int _col, int _row) {
	  this.nutrientField[_col][_row] = 0f;
  }

  int[] getCurrentTile(PVector loc) {
    

    int col = (int) Math.min(Math.max(loc.x/resolution, 0), this.cols-1);
    int row = (int) Math.min(Math.max(loc.y/resolution, 0), this.rows-1);

    int[] tile = {col, row};
    
    return tile;
  }
	
}
