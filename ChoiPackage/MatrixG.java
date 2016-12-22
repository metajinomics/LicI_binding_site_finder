package ChoiPackage;

//***************************************************//
//  Class readFile : MatrixG                         //
//***************************************************//
//Last update: Sep./5/2014
import java.io.*;
import java.util.*;
public class MatrixG<T>{  
	private ArrayList<ArrayList<T>> data; 
	private ArrayList <T> dataList;
	
	public MatrixG() {

		dataList = new ArrayList<T> ();
		data = new ArrayList<ArrayList<T>> ();

     }//MatrixG

//----------------//
//  Methods       //
//----------------//
	
	public ArrayList<ArrayList<T>> get () {
		return data;
	}
	public ArrayList<T> getcol () {
		return dataList;
	} 
	
	public void addcol(T data){
		dataList.add(data);
  	}
  	
  	public void addrow(){
		data.add(dataList);
		dataList = new ArrayList<T> ();
  	}
  	
  	public void clear(){
  		dataList.clear();
  		data.clear();
		dataList = new ArrayList<T> ();
		data = new ArrayList<ArrayList<T>> ();
  	}
	
	public void printMatrix(){
		for (int i=0;i<data.size();i++){
			for(int j=0;j<data.get(i).size();j++){
				System.out.print(data.get(i).get(j)+" ");
			}
			System.out.println(" ");
		}	
  	}
  	
  	public void copy(ArrayList<ArrayList<T>> ori){
		for (int i=0;i<ori.size();i++){
			for(int j=0;j<ori.get(i).size();j++){
				dataList.add(ori.get(i).get(j));
			}
			data.add(dataList);
			dataList = new ArrayList<T> ();
		}	
  	}
  	
  	public void move(ArrayList<ArrayList<T>> ori){
		copy(ori);	
		ori.clear();
  	}

}//class readFile