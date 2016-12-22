//this code calculate PSSM from given seed sequnces
// usage: java calculate_pssm seed_sequences
// example: java calculate_pssm SeedCcpA.txt
import java.io.*;
import java.util.*;
import java.util.Scanner;
import java.util.Random;


import ChoiPackage.*;

class calculate_pssm{
	public static void main(String[] args) {
	
	//Check usage
    try {
    	if (args.length != 1) {
        	System.err.println("usage: java calculate_pssm seed_sequences");
        	System.exit(0);
       	}
	}catch (Exception ex) {
			ex.printStackTrace();
	}//try
	
	double startscorenum= 1; //to avoid dividing 0
	int frequency = 0;
	
	//===================// 
	// Step 1: Read file //
	//===================// 
	String seedccpAFile = args[0];
	ReadFile seedccpA = new ReadFile(seedccpAFile,",");
	//seedccpA.printMatrix();
	
	//=============// 
	// Step 2:     //
	//=============// 
	MatrixG<String> Motif1 = new MatrixG<String>();
 	Motif1 = methods.FindRiMotif(seedccpA,startscorenum,frequency);	
 	Motif1.printMatrix();
	
	}
}