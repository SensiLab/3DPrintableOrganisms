package generative.components.examples;

import java.util.ArrayList;

import generative.components.*;
//import processing.core.*;

public class Testing{
	public static void main(String[] args) {
//		PApplet.main("Testing");
//		
		
//		
//		
//		
//		long genStart = System.currentTimeMillis();
//		s.generate();
////		s.getColony3D().saveColony("C:\\Users\\camil\\01_Other Work\\Code Dev\\JavaTests\\", "outputTest001");
//		long genDone = System.currentTimeMillis() - genStart;
//		System.out.println("Fitness: generation done ("+genDone+" ms)");
//		System.out.println(s.isComplete());
//		System.out.println(s.convexity());
//		System.out.println(s.printability());
//		System.out.println(s.compactness());
//		
//		Colony3D col3d = s.getColony3D();
//		
////		for (Colony col : col3d.getLayers()) {
////			for
////		}
//		
//		
//	}
	
		//0.272272436f, 0.030558953f, 0.007127416f, 0.001067804f, 0.096573976f
		//0.382983121f, 0.343619075f, 0.455419871f, 0.102930382f, 0.662056974f
		 // //Printability 0.3
//		Chromosome c = new Chromosome(0.9948359404002058f, 0.486871528443838f, 0.900167237536065f, 0.984005146251034f, 0.11341161543211142f);
		
		// Printability 0.5
//		Chromosome c = new Chromosome(0.6769059731544329f, 0.64715821613544f, 0.831787344311771f, 0.95897969719842f, 0.10433459776762366f);
		
		// Printability 0.7
//		Chromosome c = new Chromosome(0.8688199862060835f, 0.137167605861153f, 0.807739016705708f, 0.563591815560236f, 0.4605745655195185f);
		
		// P = 0.30
//		Chromosome c = new Chromosome(0.8276204359321518f, 0.7561370192229118f, 0.22717792230091818f, 0.8641633772932387f, 0.49978244548286865f);
//		long randomSeed = -835719707812819656L;
		
		// P = 0.30 threshold = 0.4
//		Chromosome c = new Chromosome(0.8380956160822853f, 0.7452117316991451f, 0.21357897826961186f, 0.9274952228947192f, 0.401466123394427f);
//		long randomSeed = -7901120352229662547L;
		
		// P = 0.21
//		Chromosome c = new Chromosome(0.9244146461160524f, 0.816053005087995f, 0.29442386820485134f, 0.6706198376539838f, 0.3903223958507103f);
//		long randomSeed = 5790058062021088326L;
		
		// P = 0.3
//		Chromosome c = new Chromosome(0.310584834318495f, 0.8129281745353379f, 0.3817478288529372f, 0.9211542885462904f, 0.24952734664335402f);
//		int randomSeed = 2046917867;
		
		// P = 0.3
//		Chromosome c = new Chromosome(0.4321711261099972f, 0.18034756927965517f, 0.3054980995324926f, 0.8915447456741861f, 0.43713795136278966f);
//		long randomSeed = -8472561951181866436L;
		
		// P = 0.5
//		Chromosome c = new Chromosome(0.3048144759524228f, 0.7870718576746254f, 0.4893161145131731f, 0.9422068264268134f, 0.9367564373660789f);
//		int randomSeed = 95000209;
		
		// P = 0.67
//		Chromosome c = new Chromosome(0.9506315462433145f, 0.6002949067695207f, 0.6540177078355786f, 0.2614439105988119f, 0.36016599374998953f);
//		int randomSeed = -164133134;
		
		// P = 0.70 rs = 1592914707750605979L Printability threshold 0.4
//		Chromosome c = new Chromosome(0.21644576849090955f, 0.47409146416100667f, 0.21088299547328215f, 0.446006944809819f, 0.558186222322518f);
//		long randomSeed = 1592914707750605979L;
		
		// P = 0.79
//		Chromosome c = new Chromosome(0.8539204862877391f, 0.6216921788891419f, 0.08495788613655472f, 0.21614102547043051f, 0.8449273850077252f);
//		int randomSeed = -623149919;
		
//		// P = 0.006
//		Chromosome c = new Chromosome(0.12310398850101575f, 0.966879208054678f, 0.003251127591460579f, 0.004246237099604166f, 0.31588478142967263f);
//		Long randomSeed = 6943217681324779413L;
		
		// P = 0.013
//		Chromosome c = new Chromosome(0.21251506117575175f, 0.8959631415337008f, 2.838641732793625E-4f, 0.0014809155346076838f, 0.40820122051263774f);
//		Long randomSeed = -1691971534857694280L;
		
		// P = 0.6
		Chromosome c = new Chromosome(0.08898201295176912f, 0.554193187490498f, 0.9554839844980083f, 0.11216229720386414f, 0.6253394620436861f);
		Long randomSeed = -1691971534857694280L;
		
		Environment env = new Environment(600, 600, 20, 0.3f, 5, 30, 60, 20, randomSeed);
//		float energy = MathLib.map(c.get(2), 0f, 1f, 5f, 20f);
	//	System.out.println();
		Simulation s = new Simulation(c, env, 500,20);
		long startTime = System.nanoTime();
		s.generate();
		long generateTime = System.nanoTime();
		System.out.println("Generate: "+ ((float)(generateTime - startTime)/1000000));
//		float p = s.printability();
//		long printabilityTime = System.nanoTime();
//		System.out.println("Printability: "+ ((float)(printabilityTime - generateTime)/1000000));
//		float conv = s.convexity();
//		long convexityTime = System.nanoTime();
//		System.out.println("Convexity: "+ ((float)(convexityTime - printabilityTime)/1000000));
//		float comp = s.compactness();
//		long compactnessTime = System.nanoTime();
//		System.out.println("Compactness: "+ ((float)(compactnessTime - convexityTime)/1000000));
//		float angles = s.angleDispersionCoefficient();
//		System.out.println(angles);
		
		Visualisation2D vis = new Visualisation2D(s.getColony3D());
//		vis.init();
		vis.draw();
		
//		String out = "../out/gCode/20200827_tests/Printability_0013_thresh_040";
//		s.printColony(out);
//		System.out.println("Printability: "+s.printability());

	}
		

	
}
