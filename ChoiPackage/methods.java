package ChoiPackage;
import java.io.*;
import java.util.*;
import java.util.Scanner;
import java.util.Random;

public class methods
{

//----------------//
//  Methods       //
//----------------// 

	public static MatrixG<String> operonComp (MatrixG<String> searchResult,
	MatrixG<String>ExpressionGeneFile,MatrixG<String> operonList,int co){
		MatrixG<String> compResult = new MatrixG<String>();
		int flag=0;
		for	(int i=0; i<searchResult.get().size();i++){
 			for (int j=0;j<ExpressionGeneFile.get().size();j++){
 				if (searchResult.get().get(i).get(5).equals(
 				ExpressionGeneFile.get().get(j).get(0))){
 					flag = 1;
 					if(co == 1){//celC
 						if(Double.parseDouble(ExpressionGeneFile.get().get(j).get(6))<0){
 						compResult.addcol(searchResult.get().get(i).get(1));
 						compResult.addrow();
 						}
 					}else if (co == 2){//manB
 						if(Double.parseDouble(ExpressionGeneFile.get().get(j).get(6))>0){
 						compResult.addcol(searchResult.get().get(i).get(1));
 						compResult.addrow();
 						}
 					}else if (co == 3){//ccpA
 						compResult.addcol(searchResult.get().get(i).get(1));
 						compResult.addrow();
 					}// if co 	
 				}//if search result
 			}// for j
 			if (co == 4 && flag == 0){//ccpA
 				compResult.addcol(searchResult.get().get(i).get(1));
 				compResult.addrow();
 			}// if co
 			flag = 0; 
   		}//for i
		return compResult;
	}//operonComp
	
	public static MatrixG<String> operonSearch (MatrixG<String> Motifmat,
	MatrixG<String>OperonSeq,double scoreTol,int palTol,int seqSkip,int penalty,
	MatrixG<String>ExpressionGeneFile,MatrixG<String> operonList){
		MatrixG<String> searchResult = new MatrixG<String>();
		for	(int i=0; i<OperonSeq.get().size();i++){
 			//search 
 			MatrixG<String> celCsearchResult = new MatrixG<String>();
 			celCsearchResult = searchBind(Motifmat,OperonSeq.get().get(i).get(0)
 			,scoreTol,palTol,seqSkip,penalty);
 			
 			double tempscore=0.0;
 			int rownum = 0;
 			//get highest
 			for (int j=0;j<celCsearchResult.get().size();j++){
 				if (tempscore<Double.valueOf(celCsearchResult.get().get(j).get(1))){
 					tempscore=Double.valueOf(celCsearchResult.get().get(j).get(1));
 					rownum=j;
 				}
 			}
 			//add to result
 			for (int j=0;j<celCsearchResult.get().get(rownum).size();j++){
 				searchResult.addcol(celCsearchResult.get().get(rownum).get(j));
 			}
 			
 			
 			
 			searchResult.addcol(operonList.get().get(i).get(2));
 			
 			int temp = Integer.valueOf(celCsearchResult.get().get(rownum).get(3))
 			+Integer.valueOf(operonList.get().get(i).get(0));
 			searchResult.addcol(String.valueOf(temp));
 			for (int j=0;j<ExpressionGeneFile.get().size();j++){
 				if (operonList.get().get(i).get(2).equals(
 				ExpressionGeneFile.get().get(j).get(0))){
 					searchResult.addcol("expressed");
 					searchResult.addcol(ExpressionGeneFile.get().get(j).get(6));
 				}
 			}
 
 			searchResult.addrow();
   		}//for
   		
		return searchResult;
	}//operonSearch
	
	
	public static MatrixG<String> getSeqForOperon(ArrayList<ArrayList<String>> full,
	MatrixG<String> regionMat){
		
		MatrixG<String> resultSeq = new MatrixG<String>();
		
		String oneRegion = "";
		
		for (int i=0; i<regionMat.get().size();i++){
			oneRegion = full.get(0).get(0).substring(
			Integer.decode(regionMat.get().get(i).get(0)),
			Integer.decode(regionMat.get().get(i).get(1)));
			//System.out.println(regionMat.get().get(i).get(0)+" "+oneRegion);
			resultSeq.addcol(oneRegion);
			resultSeq.addrow();
		}//for i
		
		return resultSeq;
	}
	
	public static MatrixG<String> getOperonListNewAllgene (MatrixG<String> operon, 
	MatrixG<String> Allgene,int frontCut,double endCut,MatrixG<String> operonSource){
		MatrixG<String> operonList = new MatrixG<String>();
		int start = 0;
   		int end = 0;
   		int startregion = 0;
   		int endregion = 0;
   		int tempcomp = 0;

   		for (int k=0;k<operon.get().size();k++){
   			String temp = operon.get().get(k).get(4);
   			int templocus = Integer.decode(temp);
   			if (operon.get().get(k).get(5).equals("+")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = start-frontCut;
   				endregion = end;//this is difference

   				for(int j=1;j<operonSource.get().size();j++){
   					String tempAllgeneLocus = operonSource.get().get(j).get(9);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus) {
   						if (operonSource.get().get(j-1).get(5).equals("+")){
   							tempcomp = Integer.decode(operonSource.get().get(j-1).get(4));
   						}else if(operonSource.get().get(j-1).get(5).equals("-")){
   							tempcomp = Integer.decode(operonSource.get().get(j-1).get(4));
   						}//if
   						if (tempcomp>startregion){
   							if(tempcomp<start){startregion = tempcomp;}
   						}//if
   						operonList.addcol(String.valueOf(startregion));
   						operonList.addcol(String.valueOf(endregion));
   						operonList.addcol(operon.get().get(k).get(1));
   						operonList.addrow();
   					}//if
   				}//for j
  
   			}else if(operon.get().get(k).get(5).equals("-")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = start;//this is difference 
   				endregion = end+frontCut;
   	
   				for(int j=1;j<operonSource.get().size()-1;j++){
   					String tempAllgeneLocus = operonSource.get().get(j).get(9);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus) {
   						if (operonSource.get().get(j+1).get(5).equals("+")){
   							tempcomp = Integer.decode(operonSource.get().get(j+1).get(3));
   						}else if(operonSource.get().get(j+1).get(5).equals("-")){
   							tempcomp = Integer.decode(operonSource.get().get(j+1).get(3));
   					}//if
   					if (tempcomp<endregion){
   						if(tempcomp>end){endregion = tempcomp;}
   					}//if
   					operonList.addcol(String.valueOf(startregion));
   					operonList.addcol(String.valueOf(endregion));
   					operonList.addcol(operon.get().get(k).get(1));
   					operonList.addcol(operon.get().get(k).get(5));
   					operonList.addrow();
   				}//if
   			}//for j
   			
   			
   			}//if
   		}//for k
   		return operonList;
	}	
	
	
	public static MatrixG<String> getOperonListNew (MatrixG<String> operon, 
	MatrixG<String> Allgene,int frontCut,double endCut,MatrixG<String> operonSource){
		MatrixG<String> operonList = new MatrixG<String>();
		int start = 0;
   		int end = 0;
   		int startregion = 0;
   		int endregion = 0;
   		int tempcomp = 0;

   		for (int k=0;k<operon.get().size();k++){
   			String temp = operon.get().get(k).get(4);
   			int templocus = Integer.decode(temp);
   			if (operon.get().get(k).get(5).equals("+")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = start-frontCut; 
   				endregion = start+geneSize;//this is difference

   				for(int j=1;j<operonSource.get().size();j++){
   					String tempAllgeneLocus = operonSource.get().get(j).get(9);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus) {
   						if (operonSource.get().get(j-1).get(5).equals("+")){
   							tempcomp = Integer.decode(operonSource.get().get(j-1).get(4));
   						}else if(operonSource.get().get(j-1).get(5).equals("-")){
   							tempcomp = Integer.decode(operonSource.get().get(j-1).get(4));
   						}//if
   						if (tempcomp>startregion){
   							if(tempcomp<start){startregion = tempcomp;}
   						}//if
   						operonList.addcol(String.valueOf(startregion));
   						operonList.addcol(String.valueOf(endregion));
   						operonList.addcol(operon.get().get(k).get(1));
   						operonList.addrow();
   					}//if
   				}//for j
  
   			}else if(operon.get().get(k).get(5).equals("-")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = end-geneSize;//this is difference
   				endregion = end+frontCut;
   	
   				for(int j=1;j<operonSource.get().size()-1;j++){
   					String tempAllgeneLocus = operonSource.get().get(j).get(9);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus) {
   						if (operonSource.get().get(j+1).get(5).equals("+")){
   							tempcomp = Integer.decode(operonSource.get().get(j+1).get(3));
   						}else if(operonSource.get().get(j+1).get(5).equals("-")){
   							tempcomp = Integer.decode(operonSource.get().get(j+1).get(3));
   					}//if
   					if (tempcomp<endregion){
   						if(tempcomp>end){endregion = tempcomp;}
   					}//if
   					operonList.addcol(String.valueOf(startregion));
   					operonList.addcol(String.valueOf(endregion));
   					operonList.addcol(operon.get().get(k).get(1));
   					operonList.addcol(operon.get().get(k).get(5));
   					operonList.addrow();
   				}//if
   			}//for j
   			
   			
   			}//if
   		}//for k
   		return operonList;
	}

	public static MatrixG<String> getOperonList (MatrixG<String> operon, 
	MatrixG<String> Allgene,int frontCut,double endCut,MatrixG<String> operonSource){
		MatrixG<String> operonList = new MatrixG<String>();
		int start = 0;
   		int end = 0;
   		int startregion = 0;
   		int endregion = 0;
   		int tempcomp = 0;

   		for (int k=0;k<operon.get().size();k++){
   			String temp = operon.get().get(k).get(4);
   			int templocus = Integer.decode(temp);
   			if (operon.get().get(k).get(5).equals("+")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = start-frontCut;
   				endregion = start+geneSize;
   				for(int j=1;j<Allgene.get().size();j++){
   					String tempAllgeneLocus = Allgene.get().get(j).get(0);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus+1) {
   						if (Allgene.get().get(j).get(4).equals("f")){
   							tempcomp = Integer.decode(Allgene.get().get(j).get(3));
   						}else if(Allgene.get().get(j).get(4).equals("r")){
   							tempcomp = Integer.decode(Allgene.get().get(j).get(2));
   						}//if
   						if (tempcomp>startregion){
   						startregion = tempcomp;
   						}//if
   						operonList.addcol(String.valueOf(startregion));
   						operonList.addcol(String.valueOf(endregion));
   						operonList.addcol(operon.get().get(k).get(1));
   						operonList.addrow();
   					}//if
   				}//for j
   			}else if(operon.get().get(k).get(5).equals("-")){
   				start = Integer.decode(operon.get().get(k).get(2));
   				end = Integer.decode(operon.get().get(k).get(3));
   				int geneSize = end - start;
   				double tempdouble = (double) geneSize;
   				tempdouble = tempdouble * endCut;
   				geneSize = (int) tempdouble;
   				startregion = end-geneSize;
   				endregion = end+frontCut;
   				for(int j=1;j<Allgene.get().size();j++){
   					String tempAllgeneLocus = Allgene.get().get(j).get(0);
   					int tempAllLocus = Integer.decode(tempAllgeneLocus);
   					if (templocus == tempAllLocus-1) {
   						if (Allgene.get().get(j).get(4).equals("f")){
   							tempcomp = Integer.decode(Allgene.get().get(j).get(2));
   						}else if(Allgene.get().get(j).get(4).equals("r")){
   							tempcomp = Integer.decode(Allgene.get().get(j).get(3));
   					}//if
   					if (tempcomp<startregion){
   						startregion = tempcomp;
   					}//if
   					operonList.addcol(String.valueOf(startregion));
   					operonList.addcol(String.valueOf(endregion));
   					operonList.addcol(operon.get().get(k).get(1));
   					operonList.addcol(operon.get().get(k).get(5));
   					operonList.addrow();
   				}//if
   			}//for j
   			}//if
   		}//for k
   		return operonList;
	}
	
	public static MatrixG<String> MatrixAddColAllgene (MatrixG<String> temp, MatrixG<String> operon){
					if (temp.get().get(0).get(5).equals("+")){
						int last = temp.get().size();
						operon.addcol(temp.get().get(0).get(0));//operon number
						operon.addcol(temp.get().get(0).get(2));//locus
						operon.addcol(temp.get().get(0).get(3));//start
						operon.addcol(temp.get().get(last-1).get(4));//end
						operon.addcol(temp.get().get(0).get(9));//order number
				
						if (temp.get().get(0).get(5).equals("+") || temp.get().get(0).get(5).equals("+1")){
							operon.addcol("+");
						}else{operon.addcol("-");}//if
						
						for(int j=0;j<temp.get().size();j++){
							operon.addcol(temp.get().get(j).get(2));
						}//for j
						operon.addrow();
					}else{
						int last = temp.get().size();
						//System.out.println(last);
						operon.addcol(temp.get().get(last-1).get(0));//operon number
						operon.addcol(temp.get().get(last-1).get(2));//locus
						operon.addcol(temp.get().get(0).get(3));//start
						operon.addcol(temp.get().get(last-1).get(4));//end
						operon.addcol(temp.get().get(last-1).get(9));//order number
				
						if (temp.get().get(last-1).get(5).equals("+") || temp.get().get(last-1).get(5).equals("+1")){
							operon.addcol("+");
						}else{operon.addcol("-");}//if
						
						for(int j=last;j>0;j=j-1){
							operon.addcol(temp.get().get(j-1).get(2));
						}//for j
						operon.addrow();
					}//if
		return operon;
	}

	public static MatrixG<String> GetOperonNewAllgene (MatrixG<String> operonSource){
		MatrixG<String> operon = new MatrixG<String>();
		MatrixG<String> temp = new MatrixG<String>();
		String operonTemp = "";
		for (int i=1; i<operonSource.get().size();i++){
		
			if (operonTemp.equals(operonSource.get().get(i).get(0))){//in case more gene in the operon
				for(int j =0; j<operonSource.get().get(i).size();j++){
					temp.addcol(operonSource.get().get(i).get(j));
				}//for j
				temp.addrow();
				operonTemp = operonSource.get().get(i).get(0);
			}else if (operonTemp != operonSource.get().get(i).get(0)){//in case no more gene in the operon
				if(operonTemp != ""){
					operon =  MatrixAddColAllgene(temp,operon);
				}
				temp = new MatrixG<String>();
				for(int j =0; j<operonSource.get().get(i).size();j++){
					temp.addcol(operonSource.get().get(i).get(j));
				}
				temp.addrow();
				operonTemp = operonSource.get().get(i).get(0);
			}//if
		}//for i
		operon =  MatrixAddColAllgene(temp,operon);		
		return operon;
	}	
	
	public static MatrixG<String> MatrixAddCol (MatrixG<String> temp, MatrixG<String> operon){
					if (temp.get().get(0).get(5).equals("+")){
						operon.addcol(temp.get().get(0).get(0));//operon number
						operon.addcol(temp.get().get(0).get(2));//locus
						operon.addcol(temp.get().get(0).get(3));//start
						operon.addcol(temp.get().get(0).get(4));//end
						operon.addcol(temp.get().get(0).get(9));//order number
				
						if (temp.get().get(0).get(5).equals("+") || temp.get().get(0).get(5).equals("+1")){
							operon.addcol("+");
						}else{operon.addcol("-");}//if
						
						for(int j=0;j<temp.get().size();j++){
							operon.addcol(temp.get().get(j).get(2));
						}//for j
						operon.addrow();
					}else{
						int last = temp.get().size();
						//System.out.println(last);
						operon.addcol(temp.get().get(last-1).get(0));//operon number
						operon.addcol(temp.get().get(last-1).get(2));//locus
						operon.addcol(temp.get().get(last-1).get(3));//start
						operon.addcol(temp.get().get(last-1).get(4));//end
						operon.addcol(temp.get().get(last-1).get(9));//order number
				
						if (temp.get().get(last-1).get(5).equals("+") || temp.get().get(last-1).get(5).equals("+1")){
							operon.addcol("+");
						}else{operon.addcol("-");}//if
						
						for(int j=last;j>0;j=j-1){
							operon.addcol(temp.get().get(j-1).get(2));
						}//for j
						operon.addrow();
					}//if
		return operon;
	}
			
	public static MatrixG<String> GetOperonNew (MatrixG<String> operonSource){
		MatrixG<String> operon = new MatrixG<String>();
		MatrixG<String> temp = new MatrixG<String>();
		String operonTemp = "";
		for (int i=1; i<operonSource.get().size();i++){
		
			if (operonTemp.equals(operonSource.get().get(i).get(0))){
				for(int j =0; j<operonSource.get().get(i).size();j++){
					temp.addcol(operonSource.get().get(i).get(j));
				}//for j
				temp.addrow();
				operonTemp = operonSource.get().get(i).get(0);
			}else if (operonTemp != operonSource.get().get(i).get(0)){
				if(operonTemp != ""){
				
					operon =  MatrixAddCol(temp,operon);
				}
				temp = new MatrixG<String>();
				for(int j =0; j<operonSource.get().get(i).size();j++){
					temp.addcol(operonSource.get().get(i).get(j));
				}
				temp.addrow();
				operonTemp = operonSource.get().get(i).get(0);
			}//if
		}//for i
		operon =  MatrixAddCol(temp,operon);		
		return operon;
	}


	public static MatrixG<String> GetOperon (MatrixG<String> operonSource){
		MatrixG<String> operon = new MatrixG<String>();
		String operonTemp = "";
		for (int i=1; i<operonSource.get().size();i++){
			//System.out.println(operonTemp+" "+operonSource.get().get(i).get(0));
			if (operonTemp.equals(operonSource.get().get(i).get(0))){
				//operon.addcol(operonSource.get().get(i).get(0));
				operon.addcol(operonSource.get().get(i).get(2));
				operonTemp = operonSource.get().get(i).get(0);
			}else if (operonTemp != operonSource.get().get(i).get(0)){
				if(operonTemp != ""){operon.addrow();}
				operon.addcol(operonSource.get().get(i).get(0));
				operon.addcol(operonSource.get().get(i).get(2));
				operon.addcol(operonSource.get().get(i).get(3));
				operon.addcol(operonSource.get().get(i).get(4));
				operon.addcol(operonSource.get().get(i).get(9));//locus
				
				if (operonSource.get().get(i).get(5).equals("+") || operonSource.get().get(i).get(5).equals("+1")){
					operon.addcol("+");
				}else{
					operon.addcol("-");
				}
				
				operonTemp = operonSource.get().get(i).get(0);
			}//if
		}//for
		operon.addrow();
		return operon;
	}

	public static MatrixG<String> randomGenerator (MatrixG<String> motifFreq){
		String[] data={"A","G","C","T"};
		MatrixG<String> randomMatrix = new MatrixG<String>();
		MatrixG<String> randomMatrixResult = new MatrixG<String>();
		for (int i=0;i<motifFreq.get().size();i++){
			int[] weight={Integer.decode(motifFreq.get().get(i).get(0)),
					Integer.decode(motifFreq.get().get(i).get(1)),
					Integer.decode(motifFreq.get().get(i).get(2)),
					Integer.decode(motifFreq.get().get(i).get(3))};   
					
			Random random=new Random();
  			String randomtemp=null;
  			for (int j=0;j<100;j++){
       				randomtemp=nextRandom(data,weight,random);
      				randomMatrix.addcol(randomtemp);
      		}//for j
      		randomMatrix.addrow();
		}// for i
		
		for (int i=0;i<randomMatrix.get().get(0).size();i++){
			String TempString = "";
			for(int j=0;j<randomMatrix.get().size();j++){
				TempString = TempString + randomMatrix.get().get(j).get(i);
			}
			randomMatrixResult.addcol(TempString);
			randomMatrixResult.addrow();
		}//for	i
		
		return randomMatrixResult;
	}
	
	public static String nextRandom(String[] data,int[] weight,Random random){
		int totalWeight = sum(weight);
    	int n=random.nextInt(totalWeight);
		int runningTotal=0;
    	for (int i=0;i<weight.length;i++){
      		runningTotal+=weight[i];
			if (n<runningTotal) return data[i];
		}
    	return null; // impossible to get here
	}

	public static int sum(int[] a){
    	int s=0;
    	for (int i=0;i<a.length;i++) s+=a[i];
    	return s;
	}

	public static String[] compareRatio(MatrixG<String> MotifcelC, MatrixG<String> MotifmanB,MatrixG<String> celCsearchResult, 
	MatrixG<String> manBsearchResult){
		String[] forReturn = new String[13];
		double possiblecelC = possiblehigh(MotifcelC);
 		double possiblemanB = possiblehigh(MotifmanB);
 		double celChigh=0;
 		double manBhigh=0;
 		for(int i=0;i<celCsearchResult.get().size();i++){
			if (celChigh<Double.parseDouble(celCsearchResult.get().get(i).get(1))){
				celChigh = Double.parseDouble(celCsearchResult.get().get(i).get(1));
				for(int j=0;j<celCsearchResult.get().get(i).size();j++){
					forReturn[j]=celCsearchResult.get().get(i).get(j);
				}
			}//if
		}
		for(int i=0;i<manBsearchResult.get().size();i++){
			if (manBhigh<Double.parseDouble(manBsearchResult.get().get(i).get(1))){
				manBhigh = Double.parseDouble(manBsearchResult.get().get(i).get(1));
				for(int j=0;j<manBsearchResult.get().get(i).size();j++){
					forReturn[j+5]=manBsearchResult.get().get(i).get(j);
				}
			}//if
		}
		
		double celCRatio = celChigh/possiblecelC;
		double manBRatio = manBhigh/possiblemanB;
		if (celCRatio > manBRatio){
			forReturn[10]="1";
		}else{
			forReturn[10]="2";
		} 
		forReturn[11]=String.valueOf(celCRatio);
		forReturn[12]=String.valueOf(manBRatio);
		return forReturn;
	}

	public static double possiblehigh(MatrixG<String> motifRegion) {
		
      //Possible highest score
      double possibleScore = 0.0;
		for (int i=0;i<motifRegion.get().size();i++){
       		if(motifRegion.get().get(i).get(0).equals("A")){
       			possibleScore=possibleScore+Double.parseDouble(motifRegion.get().get(i).get(1));
       		}else if(motifRegion.get().get(i).get(0).equals("G")){
       			possibleScore=possibleScore+Double.parseDouble(motifRegion.get().get(i).get(2));
       		}else if(motifRegion.get().get(i).get(0).equals("C")){
       			possibleScore=possibleScore+Double.parseDouble(motifRegion.get().get(i).get(3));
       		}else if(motifRegion.get().get(i).get(0).equals("T")){
       			possibleScore=possibleScore+Double.parseDouble(motifRegion.get().get(i).get(4));
       		}else {
       		}//if
		}//for 
		return possibleScore;
	}

	public static MatrixG<String> getSeq(ArrayList<ArrayList<String>> full,
	ArrayList<ArrayList<String>> expression,int frontCut, double endCut){
		MatrixG<Integer> regionMat = new MatrixG<Integer>();
		MatrixG<String> resultSeq = new MatrixG<String>();
		regionMat = ExRegion(expression,frontCut,endCut);
		String oneRegion = "";
		
		for (int i=0; i<regionMat.get().size();i++){
			oneRegion = full.get(1).get(0).substring(regionMat.get().get(i).get(0),
			regionMat.get().get(i).get(1));
			
			resultSeq.addcol(oneRegion);
			resultSeq.addrow();
		}//for i
		
		return resultSeq;
	}
	public static MatrixG<Integer> ExRegion(ArrayList<ArrayList<String>> expression,
	int frontCut, double endCut){
		int start = 0;
 		int end = 0;
 		int geneSize = 0;
 		int startregion = 0;
 		int endregion = 0;
 		int[] region = new int[2];
 		MatrixG<Integer> regionMat = new MatrixG<Integer>();
 		for (int i=0; i<expression.size();i++){
 			start = Integer.decode(expression.get(i).get(2));
 			end = Integer.decode(expression.get(i).get(3));
 			
 			if (expression.get(i).get(1).equals("f")){
 				geneSize = end - start;
 				double tempdouble = (double) geneSize;
 				tempdouble = tempdouble * endCut;
 				geneSize = (int) tempdouble;
 				startregion = start-frontCut;
 				if (startregion <0){startregion =0;}
 				endregion = start+geneSize;
 				regionMat.addcol(startregion);
 				regionMat.addcol(endregion);
 				region[0] = startregion;
 				region[1] = endregion;
 			} else if (expression.get(i).get(1).equals("r")){
 				geneSize = start - end;
				double tempdouble = (double) geneSize;
				tempdouble = tempdouble * endCut;
				geneSize = (int) tempdouble;
				startregion = start-geneSize;
				if (startregion <0){startregion =0;}
				endregion = start+frontCut;
				regionMat.addcol(startregion);
 				regionMat.addcol(endregion);
 				region[0] = startregion;
 				region[1] = endregion;
 			}//if
 			regionMat.addrow();
 		}//for i
 		
		return regionMat;
	}
	


	public static MatrixG<String> ExPlSc(ArrayList<ArrayList<String>> expression,MatrixG<String>compResult){
		MatrixG<String> ExpressionPlusScore = new MatrixG<String>();
			int highestJ=0;
			double temp=0.0;
			for (int i=0; i<expression.size();i++){
			highestJ=0;
			temp=0.0;
				for(int j=0;j<compResult.get().size();j++){
					
					if(expression.get(i).get(0).equals(compResult.get().get(j).get(0))){
						if(temp<Double.parseDouble(compResult.get().get(j).get(8))){
							temp=Double.parseDouble(compResult.get().get(j).get(8));
							highestJ=j;
						}//if
					}//if
				}//for j
				for(int k=0;k<expression.get(i).size();k++){
 						ExpressionPlusScore.addcol(expression.get(i).get(k));
 				}//for k
 				ExpressionPlusScore.addcol(compResult.get().get(highestJ).get(7));
 				ExpressionPlusScore.addcol(compResult.get().get(highestJ).get(8));
 				ExpressionPlusScore.addrow();
			}//for i
		return ExpressionPlusScore;
	}
	
	
	public static MatrixG<String> compEx(ArrayList<ArrayList<String>> expression,
	MatrixG<String> searchResult,int frontCut, double endCut){
		MatrixG<Integer> regionMat = new MatrixG<Integer>();
		regionMat = ExRegion(expression,frontCut,endCut);
	
 		int startregion = 0;
 		int endregion = 0;
 		MatrixG<String> resultMat = new MatrixG<String>();
 		for (int i=0; i<expression.size();i++){
 			startregion=regionMat.get().get(i).get(0);
 			endregion=regionMat.get().get(i).get(1);

 			for (int j=0;j<searchResult.get().size();j++){
 				if (Integer.decode(searchResult.get().get(j).get(3))>startregion && 
 				Integer.decode(searchResult.get().get(j).get(3)) <endregion){
 					for(int k=0;k<expression.get(i).size();k++){
 						resultMat.addcol(expression.get(i).get(k));
 					}
 					for(int k=0;k<searchResult.get().get(j).size();k++){
 						resultMat.addcol(searchResult.get().get(j).get(k));
 					}
 					resultMat.addrow();
 				}
 			}//for
 		}//for i
		return resultMat;
	
	}
	
	public static void writeFile(MatrixG<String> result,PrintStream resultOut){
		
		for (int i=0;i<result.get().size();i++){
			for(int j=0;j<result.get().get(i).size();j++){
				resultOut.print(result.get().get(i).get(j)+",");
			}
			resultOut.println();
		}
	
	
	}
	
	public static void writeFileNew(MatrixG<String> result,PrintStream resultOut,String seperator){
		
		for (int i=0;i<result.get().size();i++){
			for(int j=0;j<result.get().get(i).size();j++){
				resultOut.print(result.get().get(i).get(j)+seperator);
			}
			resultOut.println();
		}
	
	
	}
	
	public static MatrixG<String> FindRiMotif(ReadFile seed,double startscorenum,int frequency){	
		MatrixG<String> Motif = new MatrixG<String>();
		double[][] scoreMat;
			
			scoreMat = FindPmat(seed,startscorenum);
			
		if(frequency==0){
			//log2 +1.9647995
			for (int i=0; i<scoreMat.length; i++) {
                for (int j=0; j<scoreMat[i].length; j++) {
                	scoreMat[i][j]=Math.log(scoreMat[i][j])/Math.log(2);
                    scoreMat[i][j]=scoreMat[i][j]+1.9647995;
                }
            }
			
			Motif=finalMat(scoreMat);
		}else if(frequency==1){
			for (int i=0; i<scoreMat[0].length-1; i++) {
                for (int j=0; j<scoreMat.length-1; j++) {
                    //System.out.print(Math.round(scoreMat[j][i] * 100)+" " );
                    long Temp=Math.round(scoreMat[j][i] * 100);
                    //System.out.print(Temp+" " );
                    Motif.addcol(String.valueOf(Temp));
                }
                //System.out.println();
                Motif.addrow();
        }//for i
		
		}//if ferequency
			
		return Motif;
	}
	
	public static double[][] FindPmat (ReadFile seed,double startscorenum){
		String[][] tri;
			String seedtemp=seed.get().get(0).get(0);;
			tri=new String[seed.get().size()][seedtemp.length()];
			int roofcount=0;
			for(int i=0; i<seed.get().size(); i++){
				seedtemp= seed.get().get(i).get(0);
				for(int j=0; j<seedtemp.length();j++){
					tri[i][j]=seedtemp.substring(j,j+1);
				}
			}
			//count nucleic acid
			double[][] scoreMat;
			scoreMat=new double[5][tri[1].length+1];
			
			double scoreA=startscorenum;
			double scoreG=startscorenum;
			double scoreC=startscorenum;
			double scoreT=startscorenum;
			
			
			for (int i=0; i<tri[0].length;i++){
				scoreA=startscorenum;
				scoreG=startscorenum;
				scoreC=startscorenum;
				scoreT=startscorenum;
				for (int j=0;j<tri.length;j++){
					if (tri[j][i].equals("A")) {
						scoreA=scoreA+1;
					}else if (tri[j][i].equals("G")) {
						scoreG=scoreG+1;
					}else if (tri[j][i].equals("C")) {
						scoreC=scoreC+1;
					} else if (tri[j][i].equals("T")) {
						scoreT=scoreT+1;
					} 
				}
				scoreMat[0][i]=scoreA;
				scoreMat[1][i]=scoreG;
				scoreMat[2][i]=scoreC;
				scoreMat[3][i]=scoreT;
				scoreMat[4][i]=scoreMat[0][i]+scoreMat[1][i]+scoreMat[2][i]+scoreMat[3][i];
			}
			
			//calculate score hori sum, calculate total number of nucleic acid
            for (int i=0; i<scoreMat.length; i++) {
                for (int j=0; j<scoreMat[i].length-1; j++) {
                    scoreMat[i][scoreMat[0].length-1]=scoreMat[i][scoreMat[0].length-1]+scoreMat[i][j];
                }
            }
		
            
    		//cal  divide, calculate frequency, p
    		for (int i=0; i<scoreMat.length; i++) {
                for (int j=0; j<scoreMat[i].length; j++) {
                    scoreMat[i][j]=scoreMat[i][j]/scoreMat[4][j];
                }
            }
            
            return scoreMat;
	}
	
	public static MatrixG<String> finalMat(double[][] scoreMat){
		MatrixG<String> Motif = new MatrixG<String>();
		String lar="A";
			//print out
			for (int i=0; i<scoreMat[0].length-1; i++) {
				if (scoreMat[0][i]>scoreMat[1][i] && scoreMat[0][i]>scoreMat[2][i] && scoreMat[0][i]>scoreMat[3][i]) {
					lar="A";
				}else if (scoreMat[1][i]>scoreMat[2][i] && scoreMat[1][i]>scoreMat[3][i]) {
					lar="G";
				}else if (scoreMat[2][i]>scoreMat[3][i]) {
					lar="C";
				}else {
					lar="T";
				}
				//motifOut.println(lar+","+scoreMat[0][i]+","+scoreMat[1][i]+","+scoreMat[2][i]+","+scoreMat[3][i]);
				String Temp[]={lar,String.valueOf(scoreMat[0][i]),String.valueOf(scoreMat[1][i]),String.valueOf(scoreMat[2][i]),String.valueOf(scoreMat[3][i])};
				for (int j=0;j<Temp.length;j++){
				Motif.addcol(Temp[j]);
				}
				Motif.addrow();
			}//for
		return Motif;
	}
	
	public static MatrixG<String> FindMotif(ReadFile seed,double startscorenum,int frequency){	
		MatrixG<String> Motif = new MatrixG<String>();
		try{
			//make matrix, read sequence
			double[][] scoreMat;
			
			scoreMat = FindPmat(seed,startscorenum);
			
		if(frequency==0){            
            //cal divide and log2, p/f and log2(p/f)
            for (int i=0; i<scoreMat.length; i++) {
                for (int j=0; j<scoreMat[i].length-1; j++) {
                    scoreMat[i][j]=scoreMat[i][j]/scoreMat[i][scoreMat[0].length-1];
					scoreMat[i][j]=Math.log(scoreMat[i][j])/Math.log(2);
                }
            }
			Motif=finalMat(scoreMat);
		}else if(frequency==1){
			for (int i=0; i<scoreMat[0].length-1; i++) {
                for (int j=0; j<scoreMat.length-1; j++) {
                    //System.out.print(Math.round(scoreMat[j][i] * 100)+" " );
                    long Temp=Math.round(scoreMat[j][i] * 100);
                    //System.out.print(Temp+" " );
                    Motif.addcol(String.valueOf(Temp));
                }
                //System.out.println();
                Motif.addrow();
            }//for i
		
		}//if ferequency
		
		
	
		}catch (Exception ex) {
		ex.printStackTrace();
		}//try
		return Motif;
	}//public
	
	public static void MatrixPrintout (double[][] scoreMat){
			//MatrixG<String> printout
            for (int i=0; i<scoreMat.length; i++) {
                for (int j=0; j<scoreMat[i].length; j++) {
                    System.out.print(" "+scoreMat[i][j]);
                }
                System.out.println("");
            }
	}      
	
	public static MatrixG<String> searchBind (MatrixG<String> Motif,String sequence,
	double scoreTol,int palTol,int seqSkip,int penalty){
		MatrixG<String> searchResult = new MatrixG<String>();
		int palscore=0;
		int start=0;
		if (seqSkip==1){start=1;}
		double score=0;
		double minusScore=0;
		double totScore=0;
		double totMinusScore=0;
		System.out.println("Progress");
		//for (int i=start; i<sequence.size();i++){//if sequence has many lines
		String seq = sequence;

		String tempSeq="";	
		String revTempSeq="";	
		
			for (int j=0;j<seq.length()-Motif.get().size()+1;j++){//until end of seq 
				score=0;
				totScore=0;
				minusScore=0;
				totMinusScore=0;
				tempSeq = seq.substring(j+0,j+Motif.get().size());
				//System.out.println(tempSeq);
				//forward score
				score = calScore(tempSeq,Motif);
				//System.out.println(score);
				
				//reverse score
				revTempSeq=revComp(tempSeq);
				//System.out.println(revTempSeq);
				minusScore = calScore(revTempSeq,Motif);
				//System.out.println(minusScore);
				//System.out.println(minusScore);
				//System.out.println(minusScore);
				
				//check palindrom
				palscore=pal(tempSeq);
				if(penalty==1){
					double pen=penalty(tempSeq,palscore);
					score=score+pen;
					minusScore=minusScore+pen;
				//System.out.println(pen);	
				}
			
				if(score>=minusScore && score>scoreTol && palscore >palTol){
					searchResult.addcol(tempSeq);
					searchResult.addcol(String.valueOf(score));
					searchResult.addcol(String.valueOf(palscore));
					searchResult.addcol(String.valueOf(j+1));
					searchResult.addcol("bind+");
					
					System.out.print("\33[1A\33[2K");
					System.out.println((j+1)+"/"+seq.length()+":"+(j+1)*100/seq.length()+"%");
					
				}else if(score<minusScore && minusScore>scoreTol && palscore >palTol){
					searchResult.addcol(revTempSeq);
					searchResult.addcol(String.valueOf(minusScore));
					searchResult.addcol(String.valueOf(palscore));
					searchResult.addcol(String.valueOf(j+1));
					searchResult.addcol("bind-");	
					
					System.out.print("\33[1A\33[2K");
					System.out.println((j+1)+"/"+seq.length()+":"+(j+1)*100/seq.length()+"%");
					
				}//if
				//System.out.println(score);
				if(searchResult.getcol().size() > 0){
					searchResult.addrow();
				}
			}//for j
			System.out.print("\33[1A\33[2K");
		//}//for i
		return searchResult;
	}  
	
	public static double penalty (String seq,int palscore){
		MatrixG<Double> penaltyMat = new MatrixG<Double>();
		double p=0;
		double temptotal=0;
		double average=0;
		double pen=0;
		if ( seq.length() % 2 == 0 ) { //when seq even mumber
			penaltyMat=calPenalty(seq.length()/2+1);
		} else { //when seq odd mumber
			penaltyMat=calPenalty(seq.length()/2);
		}//if
		pen=penaltyMat.get().get(3).get(palscore);
		return pen;
	} 

// This is mathematics for calculating penalty	
	public static MatrixG<Double> calPenalty (int length){
		MatrixG<Double> penaltyMat = new MatrixG<Double>();
		int q=0;
		double average=0;
		double temptotal=0;
		//System.out.println(length);
		for (int i=0;i<length;i++){
				penaltyMat.addcol(nCr(length-1,i)*Math.pow(0.25,i)*Math.pow(0.75,length-1-i));	
			} //for i
			penaltyMat.addrow();
		for (int i=0;i<length;i++){
				temptotal=temptotal+penaltyMat.get().get(0).get(i);
			}
			average=temptotal/(length);
			//generalize
		for (int i=0;i<length;i++){
				penaltyMat.addcol(penaltyMat.get().get(0).get(i)/average);
			}
			penaltyMat.addrow();
			//inverse
		for (int i=0;i<length;i++){
				penaltyMat.addcol(1/penaltyMat.get().get(1).get(i));
			}
			penaltyMat.addrow();
			//log2
		for (int i=0;i<length;i++){
				penaltyMat.addcol(Math.log(penaltyMat.get().get(2).get(i))/Math.log(2));
			}
			penaltyMat.addrow();
			//penaltyMat.printMatrix();
			//System.out.println();
			//System.out.println();
			//System.out.println(nCr(7,6));
			//System.out.println();
			//System.out.println();
		return penaltyMat;
	} 
	
	public static int nCr(int n,int r){
		return factorial(n)/(factorial(r)*factorial(n-r));
	}
	
	public static int factorial(int n){
		return n==0?1:n*factorial(n-1);
	}
	  
	public static int pal (String seq){
		String revTempSeq=revComp(seq);
		int palscore=0;
		
		if ( seq.length() % 2 == 0 ) { //when seq even mumber
			for (int i=0;i<seq.length()/2;i++){
				if(seq.substring(i+0,i+1).equals(revTempSeq.substring(i+0,i+1))){
					palscore++;
				}//if
			} 
		} else { //when seq odd mumber
			for (int i=0;i<(seq.length()-1)/2;i++){
				if(seq.substring(i+0,i+1).equals(revTempSeq.substring(i+0,i+1))){
					palscore++;
				}//if
			} 
		}
		return palscore;
	} 
		
	public static double calScore (String tempSeq,MatrixG<String> Motif){
		double score=0;	
		for (int k=0;k<Motif.get().size();k++){		
			if (tempSeq.substring(0+k,1+k).equals("A")){
				score=score+Double.parseDouble(Motif.get().get(k).get(1));
			} else if(tempSeq.substring(0+k,1+k).equals("G")){
				score=score+Double.parseDouble(Motif.get().get(k).get(2));
			} else if(tempSeq.substring(0+k,1+k).equals("C")){
				score=score+Double.parseDouble(Motif.get().get(k).get(3));
			} else if(tempSeq.substring(0+k,1+k).equals("T")){
				score=score+Double.parseDouble(Motif.get().get(k).get(4));
			}	
		}//for k
		return score;
	}      
	
	public static String revComp (String tempSeq){
		String tempRev="";
		for (int i=0;i<tempSeq.length();i++){
			if(tempSeq.substring(tempSeq.length()-1-i,tempSeq.length()-i).equals("A")){
				tempRev=tempRev+"T";
			}else if (tempSeq.substring(tempSeq.length()-1-i,tempSeq.length()-i).equals("G")){
				tempRev=tempRev+"C";
			}else if (tempSeq.substring(tempSeq.length()-1-i,tempSeq.length()-i).equals("C")){
				tempRev=tempRev+"G";
			}else if (tempSeq.substring(tempSeq.length()-1-i,tempSeq.length()-i).equals("T")){
				tempRev=tempRev+"A";
			}//if		
		}
		return tempRev;
	} 
            



}