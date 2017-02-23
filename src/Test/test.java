package Test;

import java.util.ArrayList;
import java.util.List;

public class test {
	public static void main(String[] args){
//		
//		String s ="-118.313527107239 33.7552446479486  0.01 0.01 0.01 0.26 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.51 0.01 0.01 0.01 0.01";
//		String[] ss = s.split(" ");
//		System.out.println(ss.length);
//		System.out.println(ss[2]);
//		
		List<Integer> inputList = new ArrayList<Integer>();
		ArrayList<ArrayList<Integer>> resultset = new ArrayList<>();
		
		inputList.add(1);
		inputList.add(2);
		inputList.add(3);
		inputList.add(4);
		
		int k=4;
		
		resultset = perm(inputList,k);
		
		for(int i=0; i<resultset.size(); i++){
			System.out.println(resultset.get(i));
		}
		
		System.out.println(resultset.size());
		
	}
	
	private static ArrayList<ArrayList<Integer>> perm(List<Integer> inputList,int k ){
		ArrayList<ArrayList<Integer>> resList = new ArrayList<>();
		
	    ArrayList<Integer> arr = new ArrayList<>();
		permInt(inputList, resList, 0,0,arr,k );//new char[inputList.size()] -> new ArrayList<String>() 
	    return resList;
	}

	private static void permInt(List<Integer> inputList, List<ArrayList<Integer>> resList, int ind, int gnd, ArrayList<Integer> arr,int k) {
	    if(ind == k){
	    	ArrayList<Integer> clone = (ArrayList<Integer>) arr.clone();
			resList.add(clone);
	        return;
	    }
	    for(int i = gnd;i<inputList.size();++i){
	    	arr.add(inputList.get(i));
	    	permInt(inputList,resList,ind+1,i+1,arr,k);
	    	arr.remove(arr.size()-1);
	    }
	}
}	
