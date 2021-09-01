package generative.components;

import java.util.ArrayList;
import java.util.Random;
import processing.core.PVector;

public class Environment{

	long randSeed;
	Random generator;
	
	int size = 600;
	int width;
	int height;
	int cols, rows;
	int resolution = 15; //the size of the tiles (use 15 for experiments)
	float drag = 0.4f;

	long randSeedOut;
	
//	long maxSeed = 2147483647;
	
	public float[][] nutrientField;

	
	float foodGrowthP = 0.0001f;
	
	ArrayList<FoodSource> foodSources;
	
//	ArrayList<FoodSource> initFS;
	int numOfFoodSources = 5;
	int foodSourceSizeMin = 100;
	int foodSourceSizeMax = 200;
//	int totEnergyOfFoodSources = 50;
	float foodGrowthRate = 5;
	float foodDecayRate = 0.001f;
	
	int timestep;

	public Environment() {
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

//		this.initFS = new ArrayList<FoodSource>();
		this.foodSources = new ArrayList<FoodSource>();

		this.init();

	}
	
	// INIT METHOD
	void init() {
//		System.out.println(this.width+" | "+this.height);
		//check if the environment is squared (not required in new version)
//		if(this.width != this.height){
//			throw new ArithmeticException("width and height MUST be the same");
//		}
		//set width and height
		this.width = this.size;
		this.height = this.size;
		
		//declare random seed
		this.generator = new Random();
		if(this.randSeed == 0L){
			long seed = this.generator.nextLong();
			this.randSeedOut = seed;
			long r = this.randSeedOut;
      
		}else{
			this.randSeedOut = this.randSeed;
		}
		this.generator.setSeed(this.randSeedOut);
		
		//assign cols and rows
		this.cols = this.width/this.resolution;
		this.rows = this.height/this.resolution;
		
		//define number of food sources
		
		// initialise nutrient field
		this.nutrientField = new float[this.cols][this.rows];

//		this.initFS = new ArrayList<FoodSource>();
		this.foodSources = new ArrayList<FoodSource>();
	
		//define position of food sources.
		for(int i =0; i < this.numOfFoodSources; i++) {
			int fsCol = this.generator.nextInt(this.cols);
			int fsRow = this.generator.nextInt(this.rows);
			int fsSize = this.generator.nextInt(this.foodSourceSizeMax - this.foodSourceSizeMin) + this.foodSourceSizeMin;
			FoodSource fs = new FoodSource(fsCol, fsCol, fsSize);
			this.foodSources.add(fs);
		}
//		if (this.initFS.size() == 0) {
//			for (int i = 0; i< this.numOfFoodSources; i++) {
//				FoodSource fs = new FoodSource(this.generator.nextInt(cols), this.generator.nextInt(rows), this.generator.nextInt(this.foodSourceSizeMax - this.foodSourceSizeMin) + this.foodSourceSizeMin);
//				this.initFS.add(fs);
//			}
//		}
		
//		this.foodSources = initFS;
		this.timestep = 0;
	}//init

  
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
	
	public void setFoodSources(int _n) {
		this.numOfFoodSources = _n;
		this.init();
	}
	
	public void addFoodSource(FoodSource _fs) {
		this.foodSources.add(_fs);
		this.numOfFoodSources = this.foodSources.size();
	}
	
	public String getParams() {
		String size = "EnvSize:\t\t("+this.width+", "+this.height+")"+"\n";
		String res = "EnvRes:\t\t\t"+this.resolution+"\n";
		String d = "EnvDrag:\t\t"+this.drag+"\n";
		String fss = "EnvFoodSources:\t\t"+this.numOfFoodSources+"\n";
		String fsSizeMin = "EnvFoodSourcesMin:\t"+this.foodSourceSizeMin+"\n";
		String fsSizeMax = "EnvFoodSourcesMax:\t"+this.foodSourceSizeMax+"\n";
		String fsGR = "EnvFoodSourceGR:\t"+this.foodGrowthRate+"\n";
		String fsDR = "EnvFoodSourcesDR:\t"+this.foodDecayRate+"\n";
		String rSeed = "Random Seed:\t\t"+this.randSeedOut+"\n";
		return size+res+d+fss+fsSizeMin+fsSizeMax+fsGR+fsDR+rSeed;
	}
	
	
	
	public void setfoodDecayRate(float rate) {
		this.foodDecayRate = rate;
	}
	
	public void updateNutrientField() {
		// initialise nutrient field for next time step
		float[][] newNutrientField = new float[this.cols][this.rows];
		
		// populate new nutrient field with decay rate
		for(int i=0;i<this.cols;i++) {
			for(int j=0;j<this.rows;j++) {
				newNutrientField[i][j] = this.foodDecayRate * -1;
			}
		}
		
		// add nutrients from food sources
		for(FoodSource fs: this.foodSources) {
			if(fs.totalEnergy > 0) {
				float eTransfer = Math.min(fs.totalEnergy, 1f/this.foodGrowthRate);
				int col = fs.col;
				int row = fs.row;
				
				if(this.nutrientField[col][row] == 0) {
					newNutrientField[fs.col][fs.row] = eTransfer;
				}else {
					newNutrientField[fs.col][fs.row] += eTransfer;
				}
				fs.totalEnergy-=eTransfer;
			}else {
				int col = this.generator.nextInt(cols);
				int row = this.generator.nextInt(rows);
				int newFsSize = (int) fs.getInitialEnergy();
				
				fs.col = col;
				fs.row = row;
				fs.totalEnergy = newFsSize;
			}
//						
		}
		
		
		//diffuse to neighbours
		for(int i = 0; i < this.cols; i++) {
			for(int j =0; j < this.rows; j++) {
				float totalEnergyInCell = this.nutrientField[i][j] + newNutrientField[i][j];
				
				if(totalEnergyInCell > 1) {
					float excessEnergy = totalEnergyInCell - 1.0f;
					
//					diffuse to neighbours
					// find neighbours
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
					
					//diffuse
					float diffuseEnergy = excessEnergy/neighbours.size();
					for(int k=0; k < neighbours.size(); k++) {
						int col_index = neighbours.get(k)[0];
						int row_index = neighbours.get(k)[1];
						newNutrientField[col_index][row_index]+=diffuseEnergy;
					}
				}
			}
		}
		
		// Add new field to existing one
		for(int i=0; i<this.cols; i++) {
			for(int j=0; j< this.rows; j++) {
				float newFood = this.nutrientField[i][j] + newNutrientField[i][j];
//				System.out.println("Cell ("+i+", "+j+") = "+newFood);
				if(newFood > 1) {
					newFood = 1;
				}else if(newFood <= 0) {
					newFood = 0;
				}
//				System.out.println("===========After==============");
//				System.out.println("Cell ("+i+", "+j+") = "+newFood);
//				System.out.println();
				this.nutrientField[i][j] = newFood;
			}
		}
	
		
	}

	public void updateNutrientField(Colony c){
		// initialise field for next time step
		float[][] newNutrientField = new float[this.cols][this.rows];
		
//		// copy current state of nutrient field
//		for (int i = 0; i < this.cols; i++) {
//			for (int j = 0; j < this.rows; j++) {
//				newNutrientField[i][j] = this.nutrientField[i][j];
//			}
//		}
		
		//add food from food sources
		for (int i = 0; i < this.foodSources.size(); i++) {
			FoodSource fs = this.foodSources.get(i);
			int col = fs.col;
			int row = fs.row;

			if (fs.totalEnergy > 0) {
				float energyTransfer = Math.min(fs.totalEnergy, 1f/this.foodGrowthRate);
				nutrientField[col][row] += energyTransfer;
				
				fs.totalEnergy-=energyTransfer;
			}
			
			if (fs.totalEnergy <=0) {
				// Remove the empty food source
				int newFsSize = (int) fs.getInitialEnergy();
				this.foodSources.remove(i);
				
				// Create a new food source
				int newFsCol = this.generator.nextInt(cols);
				int newFsRow = this.generator.nextInt(rows);
//				int newFsSize = this.generator.nextInt(this.foodSourceSizeMax - this.foodSourceSizeMin) + this.foodSourceSizeMin;
				FoodSource newFs = new FoodSource(newFsCol, newFsRow, newFsSize);

				this.foodSources.add(newFs);
			}
		}
		
		

//		if (this.nutrientField[i][j] > 0) {
//			newNutrientField[i][j] -= this.foodDecayRate;
//		} else {
//			newNutrientField[i][j] = 0;
//		}
		
		// update
		for (int i = 0; i< this.cols; i++) {
			for (int j = 0; j< this.rows; j++) {
				
				// 
				if (this.nutrientField[i][j] <= 1f) {
					newNutrientField[i][j] = nutrientField[i][j] - this.foodDecayRate;

				}else {
					float energy = this.nutrientField[i][j] - 1;
					
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
					
					if (neighbours.size() > 0) {
						for(int[] n : neighbours) {							
							System.out.println(n[0]+" "+n[1]);
						}
					}
					
					float spillEnergy = energy/(float)neighbours.size();
					for (int[] n : neighbours) {
						System.out.println("Neighbour "+n[0]+", "+n[1]+" will receive "+spillEnergy);
						newNutrientField[n[0]][n[1]] += spillEnergy;
						System.out.println(newNutrientField[n[0]][n[1]]);
					}
					
					newNutrientField[i][j] = 1f;
					
					CommonFunctions.print2DArray(newNutrientField);
					System.out.println("========================");
					
				}
				
				CommonFunctions.print2DArray(newNutrientField);
				System.out.println("+++++++++++++++++++++++++");

			}
		}
		
		for(int i = 0; i < this.cols; i++) {
			for(int j = 0; j < this.rows; j++){
				if(newNutrientField[i][j] <= 0) {
					newNutrientField[i][j] = 0;
				}
			}
		}
		
		
		CommonFunctions.print2DArray(newNutrientField);
		
		System.out.println();
		
//		for(int i = 0; i < this.cols; i++) {
//			for(int j = 0; j < this.rows; j++) {
//				this.nutrientField[i][j] = newNutrientField[i][j];
//			}
//		}
		this.nutrientField = newNutrientField;
		
		this.timestep++;
	}//updateNutrientField

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

	public int[] getCurrentTile(PVector loc) {
    

		int col = (int) Math.min(Math.max(loc.x/resolution, 0), this.cols-1);
		int row = (int) Math.min(Math.max(loc.y/resolution, 0), this.rows-1);

		int[] tile = {col, row};
    
		return tile;
	}
	
}
