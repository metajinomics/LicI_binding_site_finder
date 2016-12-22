//this code find binding site from given binding seqeunces
//input : known binding sequences, target genome
// usage: java find_binding_site binding sequences target_genome output_file
// example: java find_binding_site SeedCcpA.txt ctResultSeq.txt output.txt
import java.io.*;
import java.util.*;
import java.util.Scanner;
import java.util.Random;


import ChoiPackage.*;

class find_binding_site{
	public static void main(String[] args) {
	
	//Check usage
    try {
    	if (args.length != 3) {
        	System.err.println("usage: java find_binding_site binding sequences target_genome output_file");
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
	
	//========================// 
	// Step 2: calculate PSSM //
	//========================// 
	MatrixG<String> Motif1 = new MatrixG<String>();
 	Motif1 = methods.FindRiMotif(seedccpA,startscorenum,frequency);	
 	//Motif1.printMatrix();
 	
 	//========================//     		
 	// 	Step 3: search bind   //
 	//========================//	
 	String fullSeq = args[1];
 	ReadFile Fullseq = new ReadFile(fullSeq,",");
 	int seqSkip = 1;//if > for first line on seq 1==yes
	int penalty = 0; 
	double scoreTol=-900000000;
	int palTol=-10;
 	MatrixG<String> searchResult;
	searchResult = new MatrixG<String>();
	searchResult = methods.searchBind(Motif1,Fullseq.get().get(0).get(0),scoreTol,palTol,seqSkip,penalty);
 		
 		
 	//========================//     		
 	// 	Step 5: write result  //
 	//========================//  	
 	String resultFile = args[2];

	try{
		PrintStream resultOut = new PrintStream(resultFile);
		resultOut.println(">>putative binding site,Score,palindrome,abs position,binding direction,Search Seq="
		+fullSeq+",Motif ="+seedccpAFile+",score Tol="+scoreTol
		+",palTol="+palTol+",startscorenum="+startscorenum);
		for (int i=0;i<searchResult.get().size();i++){
			for(int j=0;j<searchResult.get().get(0).size();j++){
			resultOut.print(searchResult.get().get(i).get(j)+",");
		}
		resultOut.println();
	}
	}catch (Exception ex) {
		ex.printStackTrace();
	}//try

	
	}
}