package generative.components.examples;

import generative.components.*;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class fitnessChange {
	public static void main(String[] args) {
		//Instantiate the manager
		CMAESManager cma = new CMAESManager();
		
		//Environment (This environment will be used for all 
		Environment e = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, -732428508956609621L);	//Environment with randomisation
//		cma.setEnvironment(e);
		
		String destFolder = "../out/CMA_Testing/MultiFitnessFunc/";
		
		//Evolve with angle dispersion
//		cma.evolve(200, 1e-5, 3, e, destFolder, "Conv_");
		//Get the result
		double[] parent = cma.getCurrentBest();
		
		//use result to evolve with convexity.
//		cma.evolveFromParent(parent, 100, 3, e, destFolder, "Conv_");
		
	}
}
