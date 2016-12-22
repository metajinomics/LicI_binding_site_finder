package ChoiPackage;

//***************************************************//
//  Class readFile : Read files                      //
//***************************************************//
//Last update: Sep./5/2014

import java.io.*;
import java.util.*;
public class ReadFile{  
	private ArrayList<ArrayList<String>> data; 

	public ReadFile(String fileName1,String Split) {
     	Scanner fileScanner = null;
		try {
			fileScanner = new Scanner(new File(fileName1));
		} catch(FileNotFoundException e) {
			System.err.println("Could not find file '" + fileName1 + "'.");
			System.exit(1);
		}//try
		String line = null;
		ArrayList <String> dataList = new ArrayList<String> ();
		data = new ArrayList<ArrayList<String>> ();
		while(fileScanner.hasNext()) {
			line = fileScanner.nextLine();
			if(line.substring(0,1).equals(">")){line = fileScanner.nextLine();}
			String[] entry = line.split(Split);
			dataList = new ArrayList<String> ();
 			for (int i=0; i<entry.length;i++){
					dataList.add(entry[i]);
				}
			data.add(dataList);
		}
       fileScanner.close();
     }//readFile()

//----------------//
//  Methods       //
//----------------//
	public ArrayList<ArrayList<String>> getreadFile () {
		return data;
	} 
	
	public ArrayList<ArrayList<String>> get () {
		return data;
	}
	
	public void printMatrix(){
		for (int i=0;i<data.size();i++){
			for(int j=0;j<data.get(i).size();j++){
				System.out.print(data.get(i).get(j)+" ");
			}
			System.out.println(" ");
		}	
  	}

}//class readFile