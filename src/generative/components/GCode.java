package generative.components;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
//import java.math.*;
import java.util.*;

public class GCode {
	
	String gCode;

	//Printer parameters
	float filament_diameter = 1.75f;
	float extrusion_multiplier = 1.0f;

	float width_table = 250.00f;
	float depth_table = 250.00f;
	float height_printer = 20.00f;
	
	//print parameters
	float path_width = 0.4f;
	float layer_height = 0.2f;
	float extruded_path_section = path_width * layer_height;
	float filament_section = (float) Math.PI *  (float) Math.pow((double)filament_diameter/2.0,2.0);

	float flow = 1.25f;
	float min_flow = 0.85f;
	float flow_increments = 0.2f;
	float print_speed;
	float fan_speed;
	float extruder_temp = 230f; //215 for PLA, 230 for PET
	float bed_temp = 85f;  //60 for PLA, 85 for PET
	
	float print_output_x = 120;
	float print_output_y = 120;
	float print_output_z = 65;

	int total_layers = 0;
	int layer = 0;

	public GCode(){
		this.gCode = ";GCODE STARTS HERE\n";
		this.startPrint();
	}
	
	public void setMultiplier(float m) {
		this.extrusion_multiplier = m;
	}
	
	public void resetMultiplier() {
		this.extrusion_multiplier = 1.0f;
	}
	
	float extrude(OVector p1, OVector p2) {
		  float points_distance = p1.dist(p2);
		  //println("Distance: "+points_distance);
		  float vol_extruded_path = extruded_path_section * points_distance;
		  //println("Volume: "+vol_extruded_path);
		  float len_extruded_path = vol_extruded_path / filament_section;
		  return len_extruded_path * extrusion_multiplier;
	}
	
	void gCommand(String command) {
		  this.gCode+=(command+"\n");
	}
	
	void setSpeed(float speed) {
		  gCommand("G1 F" + speed + " ; Set speed");
	}
	
	void setFan(float _fan_speed) {
		if (_fan_speed > 0) {
			gCommand("M106 S"+ _fan_speed + " ; SET FAN SPEED");
		} else {
			gCommand("M107 ; FAN OFF");
		}
	}

	void enableFan(boolean enable) {
		if (enable) {
			gCommand("M106 S"+ fan_speed);
		} else {
			gCommand("M107");
		}
	}
	
	void startPrint() {
		gCommand(";START SIMULATION INFO");
		gCommand(";=========================================");

		gCommand(";=========================================");
		gCommand(";END MODEL INFO");
		gCommand(";=========================================");
		gCommand(";START PRINT INFO");
		gCommand(";=========================================");
		gCommand(";PRINT HEIGHT: "+print_output_z);
		gCommand(";PRINT LAYERS: "+total_layers);

		gCommand(";=========================================");
		gCommand(";END PRINT INFO");
		gCommand(";=========================================");
		gCommand("M73 P0 R0");
		gCommand("M73 Q0 S0");
		gCommand("M201 X1000 Y1000 Z1000 E5000 ; sets maximum accelerations, mm/sec^2");
		gCommand("M203 X200 Y200 Z12 E120 ; sets maximum feedrates, mm/sec");
		gCommand("M204 P1250 R1250 T1250 ; sets acceleration (P, T) and retract acceleration (R), mm/sec^2");
		gCommand("M205 X8.00 Y8.00 Z0.40 E1.50 ; sets the jerk limits, mm/sec");
		gCommand("M205 S0 T0 ; sets the minimum extruding and travel feed rate, mm/sec");
		gCommand("M107");
		gCommand("G90 ; use absolute coordinates");
		gCommand("M83 ; extruder relative mode");
		gCommand("M104 S"+extruder_temp+" ; set extruder temp");
		gCommand("M140 S"+bed_temp+" ; set bed temp");
		gCommand("M190 S"+bed_temp+" ; wait for bed temp");
		gCommand("M109 S"+extruder_temp+" ; wait for extruder temp");
		gCommand("G28 W ; home all without mesh bed level");
		gCommand("G80 ; mesh bed leveling");
		gCommand("G1 Y-3.0 F1000.0 ; go outside print area");
		gCommand("G92 E0.0");
		gCommand("G1 X60.0 E9.0 F1000.0 ; intro line");
		gCommand("G1 X100.0 E12.5 F1000.0 ; intro line");
		gCommand("G92 E0.0");
		gCommand("M221 S95");
		gCommand("M900 K30 ; Filament gcode");
		gCommand("G21 ; set units to millimeters");
		gCommand("G90 ; use absolute coordinates");
		gCommand("M83 ; use relative distances for extrusion");
		gCommand("");
		gCommand(";==========================================");
		gCommand(";END OF PREAMBLE");
		gCommand(";==========================================");
		gCommand("");

	} //startPrint()
	
	void endPrint(OVector _ep) {
		float x = _ep.x + 5;
		float y = _ep.y + 5;
		gCommand("");
		gCommand(";====================");
		gCommand(";START OF CODA");
		gCommand(";====================");
		gCommand("");
		gCommand("M204 S1000");
		gCommand("G1 X120.611 Y100.329 F10800.000");
		gCommand("G1 F8640;_WIPE");
		gCommand("G1 X"+x+" Y"+y+" E-0.76000");
		gCommand("G1 E-0.04000 F2100.00000");
		//gCommand("G1 Z1.000 F10800.000");
		gCommand("M107");
		gCommand("; Filament-specific end gcode");
		gCommand("G4 ; wait");
		gCommand("M221 S100");
		gCommand("M104 S0 ; turn off temperature");
		gCommand("M140 S0 ; turn off heatbed");
		gCommand("M107 ; turn off fan");
		gCommand("G91 ;USE RELATIVE COORDINATES");
		gCommand("G1 Z10 ; Move print head up");
		gCommand("G90 ;USE ABSOLUTE COORDINATES");
		gCommand("G1 X0 Y200 F3000 ; home X axis");
		gCommand("M84 ; disable motors");
		gCommand("M73 P100 R0");
		gCommand("M73 Q100 S0");

	}	// endPrint()
	
	void writeGcodeLayer(float layerNumber, Colony layer) {
		
	    gCommand(" ");
	    gCommand(";-----------------");
	    gCommand("; LAYER " + (layerNumber) + " STARTS HERE");
	    gCommand(";-----------------");
	    gCommand(" ");
	    
	    
		
		float width = layer.getEnvironment().width;
		float height = layer.getEnvironment().height;
		
		float z = this.layer_height * (layerNumber + 1);
		
		gCommand("G92 E0.0");

		gCommand("G1 E-0.80000 F2100.00000");
		gCommand("G1 Z" + (z + 0.4) + " F10800.000");

		setSpeed(2500);
		OVector origin = new OVector();
		//--------------traverse all organisms
		for (Organism org : layer.organisms) {
			//-------------traverse all edges
			for (Spring s : org.springs) {
				int j = org.springs.indexOf(s);
				//Get edge
				OVector[] l = {s.sp.loc, s.ep.loc};
				//get start point
				float spX = ((l[0].x) * print_output_x) / (width);
				float spY = ((l[0].y) * print_output_y) / (height);
		      
				OVector sp = new OVector(spX, spY);

				float epX = ((l[1].x) * print_output_x) / (width);
				float epY = ((l[1].y) * print_output_y) / (height);

				OVector ep = new OVector(epX, epY);


				float extrusion = extrude(sp, ep);
				if (j == 0) {
					origin = sp;


					// go to starting point
					gCommand("G1 X" + sp.x + " Y" + sp.y + " ;This is the initial point");

					// get layer height
					gCommand("G1 Z"+ z +" ; position nozzle at layer height");

					gCommand("G1 E0.80000 F2100.00000");

					setSpeed(print_speed);


					// go to next point
					gCommand("G1 X" +ep.x + " Y" + ep.y + " E" + extrusion+ " ;This is the second point");
				} else if (j < org.springs.size() -1) {
					gCommand("G1 X" + ep.x + " Y" + ep.y + " E" + extrusion);
				} else {
					//extrusion = extrude(sp,origin,_ext_mult);
					//println("Closing the cell");
					//println("Extrusion: ", extrusion);
					//println("From:  ", sp);
					//println("To:    ", ep);
					gCommand("G1 X" + origin.x + " Y" + origin.y + " E" + extrusion);
				}

				//write the layer transition stuff


				//final_point = ep;
			}
		}

	} // write layer
	
	public void export(String fileName, OVector _ep) {

		// End print
		endPrint(_ep);
		
		// Create file
		File output = new File(fileName+".gcode");
		
		
		// Write to file
		try {
			FileWriter outputWriter = new FileWriter(fileName+".gcode");
			outputWriter.write(this.gCode);
			outputWriter.close();
		}catch(IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	

	
}
