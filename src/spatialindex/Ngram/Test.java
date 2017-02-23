package spatialindex.Ngram;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

public class Test {
	public static void main(String[] args) throws IOException {
		//---------------------------------------------------读取文件内的内容-----------------------------------------
		FileReader fr = new FileReader("E:/实验室论文cj/Venues/NYC/NYC-Venues.txt");
		BufferedReader br = new BufferedReader(fr);
		String s;
		
		Map<String, String> idString = new HashMap<String, String>();
		Map<String,Vector<String>> invertedList = new HashMap<String,Vector<String>>();
		int edThreshold = 2;//threshold门槛、入�?
		int K = 3;
		String query = "coffee shop";//查询关键�?
	
		while((s=br.readLine())!=null)
		{
			String[] temp = s.split("	");
			String id = temp[0].substring(1, temp[0].length()-1);
			String string = temp[1].substring(1, temp[1].length()-1);
			boolean flag = true;
			for(int j=0;j<string.length();j++){
				if(string.charAt(j)<0 || string.charAt(j)>=128)
					flag = false;
			}
			if(flag){
				idString.put(id, string);
			}
		}
		//--------------------------------------------------------------------------------------------------------
		
	
		for(Map.Entry<String, String> entry: idString.entrySet()){//表示访问�?有的Entry
			String[] temp = entry.getValue().split(" ");//将Entry中的数�?�读取到temp数组当中
			for(int i = 0;i < temp.length;i ++)
			{	 
				List<String> grams = new ArrayList<String>();//创建�?个新的数�?
				grams = Ngram.getKGramsList(temp[i], K);//依次取出k个�?�到数组�?
				if(grams!=null){
					for(String str: grams)//依次读取grams数组当中的�??
					{
						Vector<String> ids = new Vector<String>();
						if(!invertedList.containsKey(str)){//如果倒排文件中包含了字段
							ids.add(entry.getKey());//则从文件中读取完整的关键�?
							invertedList.put(str, ids);//将此关键字放到�?�排文件当中
						} else {//如果倒排文件中不包括此字�?
							ids = invertedList.get(str);
							ids.add(entry.getKey());
							invertedList.put(str, ids);
						}
					}
				}
			}
		}
		
		System.out.println("倒排文件的大小为�?" + invertedList.size());
		
		//-------------------------------------------------------------------------------------
		long startTime1=System.currentTimeMillis(); //�?始计�?
		int j = 0;
		int idnum = 0;
		int num1 = 0;
		for(Map.Entry<String, Vector<String>> entry: invertedList.entrySet()){
		if(query.contains(entry.getKey())){
			j++;
			Vector<String> idSet = new Vector<String>();
			idSet = entry.getValue();
			idnum += idSet.size();
			for(String id: idSet){
				String[] words = idString.get(id).split(" ");
				for(int i = 0;i < words.length;i ++){
					if(EditDistance.ld(words[i], query) <= edThreshold){
						break;
					}
				}
			}
		}
	}
		System.out.println("数据j的�?�为�?" + j);
		System.out.println("数据idnum的�?�为�?" + idnum);
		long endTime1=System.currentTimeMillis();
		System.out.println(num1+"利用N-gram+InvertedList�?�?要的时间�?"+(endTime1-startTime1)+"ms");
		//----------------------------------------------------------------------------------------------------------
		
		
		//-----------------------------------------------------------------------------------------------------------
		int num2=0;
		long startTime2=System.currentTimeMillis();
		for(Map.Entry<String,String> entry: idString.entrySet()){
			String[] strings = entry.getValue().split(" ");
			for(int i = 0;i < strings.length;i ++){
				if(EditDistance.ld(strings[i],query) <= edThreshold){
					num2++;
					break;
				}
			}
		}
		long endTime2=System.currentTimeMillis();
		System.out.println(num2+"利用线�?�扫描所�?要的时间�?"+(endTime2-startTime2)+"ms");
		//-------------------------------------------------------------------------------------------------------------
	}
}
