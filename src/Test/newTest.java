package Test;

import spatialindex.rtree.RTree;
import spatialindex.spatialindex.*;
import spatialindex.storagemanager.*;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class newTest {
	final int sumPoints = 140000;
	
	static int objective_number = 2;
	private IRPoint[] point = new IRPoint[sumPoints];//???????????
	private int countPoint = 0;//?????
	private int querycount = 0;//????????????????

	//------------------------------------------------------lsh????????----------------------------------------------------------------------------
	private ArrayList<ArrayList<double[]>> a = new ArrayList<ArrayList<double[]>>();  //LSH??a
	private double b;   //LSH??b???0~w?????????
	private ArrayList<HashMap<String,ArrayList<IRPoint>>> lshIndex = new ArrayList<HashMap<String,ArrayList<IRPoint>>>(); 
	private static double w = 0.5;   //??????????????LSH??w
	public static HashMap<Integer,IRPoint> map = new HashMap<Integer,IRPoint>();
	
	public static HashMap<Integer,HashMap<String,ArrayList> > inverted = new HashMap<>();

	public static double[] coord = new double[2];//???????????????
	
	//===============================================================================
	private QueryPoint[] querypoint = new QueryPoint[objective_number];
	public static int dimention = 20;
	private static int L = 3;  //hash???????,L
	private static int M = 8;  //???hash??????hash??????????, M
	
	double beta = 0.5;
	double alpha = 0.5;
	
//	newTest(int i)
//	{
//		querypoint = new QueryPoint[i];
//	}
	
	public static void main(String[] args) throws IOException
	{
		String QueryFile = "F:\\LSHquery\\queryset\\query"+objective_number+"_20_NYC.txt"
				//"E:\\LSHquery\\queryset\\query4_20_NYC.txt"
				//"E:\\LSHquery\\queryset\\query6_20_NYC.txt"
				//"E:\\LSHquery\\queryset\\query8_20_NYC.txt"
				/*"E:\\LSHquery\\queryset\\query10_20_NYC.txt",
				"E:\\LSHquery\\queryset\\query12_20_NYC.txt",*/
				//"E:\\LSHquery\\queryset\\query2_50_NYC.txt"
				//"E:\\LSHquery\\queryset\\query4_50_NYC.txt"
				//"E:\\LSHquery\\queryset\\query6_50_NYC.txt"
				//"E:\\LSHquery\\queryset\\query8_50_NYC.txt"
				/*"E:\\LSHquery\\queryset\\query10_50_NYC.txt",
				"E:\\LSHquery\\queryset\\query12_50_NYC.txt",*/
				//"E:\\LSHquery\\queryset\\query2_100_NYC.txt"
				//"E:\\LSHquery\\queryset\\query4_100_NYC.txt"
				//"E:\\LSHquery\\queryset\\query6_100_NYC.txt"
				//"E:\\LSHquery\\queryset\\query8_100_NYC.txt"
				/*"E:\\LSHquery\\queryset\\query10_100_NYC.txt",
				"E:\\LSHquery\\queryset\\query12_100_NYC.txt"*/
				//"E:\\LSHquery\\queryset\\query2_20_LA.txt"
				//"E:\\LSHquery\\queryset\\query4_20_LA.txt"
				//"E:\\LSHquery\\queryset\\query6_20_LA.txt"
				//"E:\\LSHquery\\queryset\\query8_20_LA.txt"
				/*"E:\\LSHquery\\queryset\\query10_20_LA.txt",
				"E:\\LSHquery\\queryset\\query12_20_LA.txt",*/
				//"E:\\LSHquery\\queryset\\query2_50_LA.txt"
				//"E:\\LSHquery\\queryset\\query4_50_LA.txt"
				//"E:\\LSHquery\\queryset\\query6_50_LA.txt"
				//"E:\\LSHquery\\queryset\\query8_50_LA.txt"
				/*"E:\\LSHquery\\queryset\\query10_50_LA.txt",
				"E:\\LSHquery\\queryset\\query12_50_LA.txt",*/
				//"E:\\LSHquery\\queryset\\query2_100_LA.txt"
				//"E:\\LSHquery\\queryset\\query4_100_LA.txt"
				//"E:\\LSHquery\\queryset\\query6_100_LA.txt"
				//"E:\\LSHquery\\queryset\\query8_100_LA.txt"
				/*"E:\\LSHquery\\queryset\\query10_100_LA.txt",
				"E:\\LSHquery\\queryset\\query12_100_LA.txt",*/;
//		String Dataset = "F:\\LSHquery\\dataset\\test.txt"
		String Dataset = "F:\\LSHquery\\dataset\\objects_NYC_20Topics.txt"
				//"E:\\LSHquery\\dataset\\objects_NYC_50Topics.txt"
				//"E:\\LSHquery\\dataset\\objects_NYC_100Topics.txt"
				//"E:\\LSHquery\\dataset\\objects_LA_20Topics.txt"
				//"E:\\LSHquery\\dataset\\objects_LA_50Topics.txt"
				//"E:\\LSHquery\\dataset\\objects_LA_100Topics.txt"
				;
		
		//for(int k=0; k<6; ++k){
			//for(int i=k*1; i<(k+1)*1; ++i){
				
				//for(int j=k; j<(k+1); ++j){
			
	//				if((i<5&&j<3) || ((i>4 && i<10)&&(j>2 && j<6)))
					//{
						String method = "IMOSKQ";//Top-K MoSKQ IMOSKQ
						newTest irtree = new newTest(/*(i%6+1)*2*/);
						if(method.equals("Top-K")){
//							System.out.println("read");
							irtree.ReadAll(Dataset);
						}else{
//							newTest irtree = new newTest(/*(i%6+1)*2*/);
							
							//System.out.println("????????LSH?????????");
						
//							Pattern p = Pattern.compile("E:\\\\LSHquery\\\\dataset\\\\objects_(.*)_(\\d*)Topics.txt");
//							Matcher m = p.matcher(Dataset[j]) ;
//							m.find();
//							newTest.dimention = new Integer(m.group(2));
//							System.out.println(newTest.dimention);
							irtree.CreateLSH(L,M,w,Dataset);
						}
						irtree.query(QueryFile,"F:\\LSHquery\\rtreedata.dat",method);
						
					//}
				
			
		
		System.out.println("?????????");
	}

	
	
	//------------------------------------createlsh---------------------------------------------------------------------------
	public void CreateLSH(int L, int M, double w, String filename){
		
		ReadAll(filename);
//		System.out.println("read points successfully!");
		
		setVarA(L,M);
//		System.out.println("set properties successfully!");
		
		Random ran = new Random(); //?????b
	    b = ran.nextDouble()*w;
		
        //---------------------LSH??????-----------------------
	    lshIndex = new ArrayList<HashMap<String,ArrayList<IRPoint>>>();
        
        for(int l = 0; l < L; l++) {
        	lshIndex.add(new HashMap<String,ArrayList<IRPoint>>());
		}

		for(int tmpnumber = 0;tmpnumber <countPoint;tmpnumber++) {//???????lshIndex??
			putLSH(point[tmpnumber],L,M,w);
		}
//		System.out.println("LSH??????????");
	}
	
	public void ReadAll(String filename)//???????????????????????????????????????keyword??topic???
	{
		//String filename = new String("E:\\LSHquery\\dataset\\objects_LA_20Topics.txt");
		
		LineNumberReader lr = null;
		try {
			lr = new LineNumberReader(new FileReader(filename));
		
		} catch (FileNotFoundException e) {
			System.err.println("Cannot open trajectory file " + filename + ".");
			System.exit(-1);
		}

		try {
			String line = lr.readLine();
			
			while((line != null) && (countPoint<sumPoints)) {
				//float[] tmp_topics = new float[dimention];//??????topic?????????????????topic????????100
				//int countTopic = 0;//????????topic??????
				//int tempnumber = 4;
					
//				StringTokenizer st = new StringTokenizer(line);
				String[] temp = line.split("	| ");

				point[countPoint] = new IRPoint();
				
				point[countPoint].m_pointID = Integer.parseInt(temp[0]); //?????????
				//point[countPoint].Keyword = temp[1];

				point[countPoint].m_pCoordinate[0] = Double.parseDouble(temp[4]);  //?????????????
				point[countPoint].m_pCoordinate[1] = Double.parseDouble(temp[5]);  //??????????????
				
				//String topicdistribution = temp[2];
				//String[] temp3 = topicdistribution.split(" ");
				for(int i=0; i<20; i++)
				{
					//while(countTopic < dimention) {
						//tmp_topics[i] = Float.parseFloat(temp3[i]);
						point[countPoint].topics[i] = Float.parseFloat(temp[i+7]);
					//}
				}
				//point[countPoint].topics = tmp_topics;

				newTest.map.put(point[countPoint].m_pointID, point[countPoint]);
				
				countPoint++;
				line = lr.readLine();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}	
	}

	private void setVarA(int L, int M) {
		Random normalRandom = new Random();
		double tmp = 0;
					
		for(int l = 0; l < L; l++) {
			a.add(new ArrayList<double[]>());
			
			for(int m = 0; m < M; m++) {//M??hash????	
				a.get(l).add(new double[dimention]);//Config.dimention???????topic????????????????????function?????????
				//a.get(l)??????hash table???hash function

				for(int k = 0; k < dimention; k++) {//??? bucket
								
					tmp = normalRandom.nextGaussian();//??????????????double?
					while(tmp < 0) {//????????????????
						tmp = normalRandom.nextGaussian();   
					}
					a.get(l).get(m)[k] = tmp;//????????????????????
				}
			}
		}
//		System.out.println("a:");
//		for(int i=0;i<a.size();++i){
//			for(int j=0;j<a.get(i).size();++j){
//				for(int k=0;k<a.get(i).get(j).length;++k)
//				System.out.println(a.get(i).get(j)[k]);
//			}
//		}
//		System.exit(0);
	}
	
	private void putLSH(IRPoint tmp, int L, int M, double W) {
//		ArrayList<Integer> bucketnumber = new ArrayList<Integer>();
//        System.out.println("tmp topics are:" + tmp.topics);
//        System.out.println("tmp ID are:" + tmp.m_pointID);
//        System.out.println();
		String key = "";

        for(int i = 0; i < L; i++) {
        	
			double b = Math.random();
			
			key = getKey(tmp.topics,i,a,b,M,W);
//			System.out.println("key" + key);
			
			tmp.bucketnumber.add(key);//??????hash table???????????????????????????????
				//????????????????????String?????
			
			if(lshIndex.get(i).get(key) == null) {//???hash?????????key?????????????key??value????
				ArrayList<IRPoint> tmpArr = new ArrayList<IRPoint>();
				tmpArr.add(tmp);
				lshIndex.get(i).put(key, tmpArr);
			}
			else {//???hash???????hash????key??????????????????
				lshIndex.get(i).get(key).add(tmp);
			}
		}
	}

	private String getKey(float[] topics, int L, ArrayList<ArrayList<double[]>> a, double b, int M, double W) {

		String result = "";//?????l??hash table??M??hash functions?????M??hash???
		
		for(int i = 0; i < M; i++) {
			int hashResult = hashFamily(topics,a.get(L).get(i),b,W,M);//a.get(l).get(i)??????l??hash table????i?????
			result += hashResult;
		}
		return result; //????hash???
		
	}

	private int hashFamily(float[] topics, double[] a, double b, double W, int M) { 
//		System.out.println("topics:");
//		for(int i=0; i< topics.length; ++i){
//			System.out.print(topics[i]);
//		}
//		System.out.println();
//		
		int h = 0;
		double tmp = b;
		
		for(int i =  0; i < dimention; i++) {
			tmp += topics[i]*a[i];
		}
		tmp = tmp/W;
		h = (int)tmp;
//		System.out.println("h"+h);
//		System.exit(0);
		return h;
	}
	
	private double guiyi(double x){
		return 2.0/(1+Math.pow(Math.E, -x)) -1;
	}
	
	private double ed(Point queryPoint,IRPoint o){//??????
		double d =0;
		for(int i = 0;i<o.m_pCoordinate.length;++i){
//			System.out.println("queryPoint:"+queryPoint.getCoord(i));
//			System.out.println("o:"+o.m_pCoordinate[i]);
			d+= Math.pow(queryPoint.getCoord(i)-o.m_pCoordinate[i], 2);
		}
		
		d = Math.sqrt(d);
		return d;
		
	}
	
	
	private List<IRPoint> allPointsInRadius(Point queryPoint2,int radius){
		List<IRPoint> l = new ArrayList<IRPoint>();
		for(int i = 0;i<point.length;++i){
			if(guiyi(ed(queryPoint2,point[i]))<=guiyi(radius)){
				l.add(point[i]);
			}
		}
		return l;
		
	}
	
	
	private List<List<IRPoint>> allCombines(List<IRPoint> allPointsInRadius,int queryCount, Point query, double alpha,double beta){
		List<List<IRPoint>> l = new ArrayList<List<IRPoint>>();
		for( int i =1;i<=queryCount;++i){
			l.addAll(allCombinesWithK(allPointsInRadius,i,query,alpha,beta));
		}
		return l;
		
	}
	
	private Collection<? extends List<IRPoint>> allCombinesWithK(
			List<IRPoint> allPointsInRadius, int k, Point query, double alpha,double beta) {
		List<List<IRPoint>> all = new ArrayList<List<IRPoint>>();
		permutationIRPoint(allPointsInRadius,all,0,new  ArrayList<IRPoint>(),k,query,alpha,beta);
		return all;
	}
	
	private void permutationIRPoint(List<IRPoint> inputList, List<List<IRPoint>> resList, int ind, ArrayList<IRPoint> arr,int k, Point query, double alpha,double beta) {
		 if(k == 0){
//			   resList.add((ArrayList<IRPoint>) arr.clone());
			 	double ds = ds(query,arr,alpha,beta);
//			 	System.out.println(ds);
			 	if(ds<cost)
			 		cost = ds;
			   return;
		   }
	   if(ind >= inputList.size())
		   return;
	  
	   for(int i = ind;i<inputList.size();++i){
		   arr.add(inputList.get(ind));
		   permutationIRPoint(inputList,resList,i+1,arr,k-1,query,alpha,beta);
		   arr.remove(arr.size()-1);
	   }
	   
	}
	
	private Collection<? extends List<IRPoint>> allCombinesWithK2(
			List<IRPoint> allPointsInRadius, int k, Point query, double alpha,double beta) {
		List<List<IRPoint>> all = new ArrayList<List<IRPoint>>();
		permutationIRPoint2(allPointsInRadius,all,0,new  ArrayList<IRPoint>(),k,query,alpha,beta);
		return all;
	}
	
	private void permutationIRPoint2(List<IRPoint> inputList, List<List<IRPoint>> resList, int ind, ArrayList<IRPoint> arr,int k, Point query, double alpha,double beta ) {
		 if(k == 0){
//			   resList.add((ArrayList<IRPoint>) arr.clone());
			 	double ds = ds2(query,arr,alpha,beta);
//			 	System.out.println(ds);
			 	if(ds<cost)
			 		cost = ds;
			   return;
		   }
	   if(ind >= inputList.size())
		   return;
	  
	   for(int i = ind;i<inputList.size();++i){
		   arr.add(inputList.get(ind));
		   permutationIRPoint2(inputList,resList,i+1,arr,k-1,query,alpha,beta);
		   arr.remove(arr.size()-1);
	   }
	   
	}
	
	
	private double ds(Point query,List<IRPoint> points,double alpha,double beta){
		double d = 0;
		
		double maxA=0;
		double maxB=0;
//		System.out.println("size:"+points.size());
		
		for(int i = 0;i<points.size();++i){ // ??query??point????????
			double d1 = guiyi(ed(query,points.get(i)));
			if(d1>maxA)
				maxA = d1;
		}
		
		//point??????????
		for(int i = 0;i<points.size();++i){
			for(int j = i+1;j<points.size();++j){
				double d1 = guiyi(ed(points.get(i),points.get(j)));
				
				if(d1>maxB)
					maxB = d1;
			}
		}
		
//		System.out.println(maxA+","+maxB);
		
		return beta*(alpha*maxA+(1-alpha)*maxB) + (1-beta)*semantics(query, points);
		
	}
	
	private double ds2(Point query,List<IRPoint> points,double alpha,double beta){
		double d = 0;
		
		double maxA=0;
		double maxB=0;
//		System.out.println("size:"+points.size());
		
		for(int i = 0;i<points.size();++i){ // ??query??point????????
			double d1 = guiyi(ed(query,points.get(i)));
			if(d1>maxA)
				maxA = d1;
		}
		
		//point??????????
		for(int i = 0;i<points.size();++i){
			for(int j = i+1;j<points.size();++j){
				double d1 = guiyi(ed(points.get(i),points.get(j)));
				
				if(d1>maxB)
					maxB = d1;
			}
		}
		
//		System.out.println(maxA+","+maxB);
		
		return beta*(alpha*maxA+(1-alpha)*maxB);
		
	}
	
	private double semantics(Point query,List<IRPoint> points){
		double d = 0;
		for(int i=0;i<querypoint.length;++i){
			double d1 = Double.MAX_VALUE;
			float[] f = querypoint[i].topics;
			for(int j=0;j<points.size();++j){
				double d2 = 0;
				for(int k =0;k<f.length;k++){
					d2+=(f[k]-points.get(j).topics[k])*(f[k]-points.get(j).topics[k]);
				}
				d2 = Math.sqrt(d2);
				d2 = guiyi(d2);
				if(d2<d1)
					d1 = d2;
			}
			d+=d1;
		}
		return d;
		
	}

	private double ed(IRPoint irPoint, IRPoint o) {
		double d =0;
		for(int i = 0;i<o.m_pCoordinate.length;++i){
			d+= Math.pow(irPoint.m_pCoordinate[i]-o.m_pCoordinate[i], 2);
		}
		
		d = Math.sqrt(d);
		return d;
	}



	private List<IRPoint> minCostCombine(Point query,List<IRPoint> allPointsInRadius,int queryCount,double alpha,double beta){
		List<List<IRPoint>> allCombines = allCombines(allPointsInRadius, querycount,query,alpha,beta);
		List<IRPoint> r = new ArrayList<IRPoint>();
//		double min = Double.MAX_VALUE;
//		int index = 0;
//		for(int i = 0;i<allCombines.size();++i){
//			double ds = ds(query,allCombines.get(i),alpha);
////			System.out.println(ds);
//			if(ds < min){
//				min = ds;
//				index = i;
//			}
//		}
////		System.out.println(min);
//		r.addAll(allCombines.get(index));
		return r;
	}
	
	//------------------------------------------------------------------------------------------------------------------------
	
	private double cost = Double.MAX_VALUE;
	private List<IRPoint> po = new ArrayList<IRPoint>();
	public void query(String queryfile,String rtreefile,String querytype)
		//????IRTree??????????
	{
		
		//--------------------------------search-----------------------------------------
		BufferedReader lr = null;
		try {
			lr = new BufferedReader(new FileReader(queryfile));
		} catch (FileNotFoundException e) {
			System.err.println("Cannot open file " + queryfile + ".");
			System.exit(-1);
		}
		
		int indexIO = 0;
		int leafIO = 0;
		
		double[] qf1 = new double[2];
		
		String line = null;
		try {
			line = lr.readLine();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		while(line != null)
		{
			float[] topic = new float[dimention];
			float temp = 0;
			int countnumber = 0;
			int number = 3;
			
			String[] querytemp2 = line.split(" ");
			
			querypoint[querycount] = new QueryPoint();
			
			//???????????
    		querypoint[querycount].m_pCoordinate[0] = Double.parseDouble(querytemp2[0]);
    		querypoint[querycount].m_pCoordinate[1] = Double.parseDouble(querytemp2[1]);
    		qf1[0] = Double.parseDouble(querytemp2[0]);
    		qf1[1] = Double.parseDouble(querytemp2[1]);
    		newTest.coord[0] = querypoint[querycount].m_pCoordinate[0];
    		newTest.coord[1] = querypoint[querycount].m_pCoordinate[1];
			
    		//???topics???
			while(countnumber < dimention) {
				temp = Float.parseFloat(querytemp2[number++]);//???????????
				topic[countnumber++] = temp;
			}
			querypoint[querycount].topics = topic;
			
			querycount++;
			
			try {
				line = lr.readLine();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		Point queryPoint = new Point(qf1);
		
		
		if(querytype.equals("Top-K")){
			long start = System.currentTimeMillis();
			System.out.println("in topk");
			
			List<IRPoint> po = new ArrayList<IRPoint>();
			List<IRPoint> optimalpo = new ArrayList<IRPoint>();
//			double cost = Double.MAX_VALUE;
			int radius = 0;
			
			
//			System.out.println(queryPoint.getCoord(0)+","+queryPoint.getCoord(1));
			
			List<IRPoint> minCostCombine = null;
			//upperbound = cost; lowerbound = 0.5*0.5*radius
			double last = 0;
			while(cost >= alpha*beta*guiyi(radius) )
			{
				radius += 5;
				List<IRPoint> allPointsInRadius = allPointsInRadius(queryPoint, radius);
				minCostCombine = minCostCombine(queryPoint,allPointsInRadius,querycount,alpha,beta);
//				cost = ds(queryPoint,minCostCombine,alpha);
				System.out.println("min:"+cost);
				if(cost == last){
					break;
				}
				last = cost;
				
			}
			
//			double min = ds(queryPoint,minCostCombine,alpha);
			System.out.println(cost);
			System.out.println("time:"+(System.currentTimeMillis()-start)/1000.0);
			
		}else{
		
		try
		{
			PropertySet ps1 = new PropertySet();

			Boolean b2 = new Boolean(true);
			ps1.setProperty("Overwrite", b2);

			ps1.setProperty("FileName", rtreefile);

			Integer i = new Integer(4096);
			ps1.setProperty("PageSize", i);
			IStorageManager diskfile = new DiskStorageManager(ps1);
			IBuffer file = new RandomEvictionsBuffer(diskfile, 10, false);

			PropertySet ps2 = new PropertySet();

			Double f = new Double(0.7);
			ps2.setProperty("FillFactor", f);

			i = new Integer(100);
			ps2.setProperty("IndexCapacity", i);   //capacity??????
			ps2.setProperty("LeafCapacity", i);

			i = new Integer(2);
			ps2.setProperty("Dimension", i);
			
			//---------------------------Create IRtree Index------------------------------
			
			ISpatialIndex tree = new RTree(ps2, file);
			
//			System.out.println("Build RTree Successfully!");
			
			int id;
			String ss;
			
			//point??????
			double[] d = new double[2];
			
			for(int j = 0;j < countPoint;j++)
			{
//				System.out.println("j??????"+j);
				id = point[j].m_pointID;
				
				d[0] = point[j].m_pCoordinate[0];//x
				d[1] = point[j].m_pCoordinate[1];//y
				Point p = new Point(d);
				
				//ss = point[j].Keyword;
				ArrayList<String> ids = point[j].bucketnumber;
				for(int k = 0;k<ids.size();++k){
					tree.insertData(ids.get(k).getBytes(), p, id);//
				}
//				tree.insertData(ss.getBytes(), p, id);//ss.getBytes()???String??????byte[]				
			}
			
//			System.err.println(tree);
			
			Integer indexID = (Integer) ps2.getProperty("IndexIdentifier");
			//System.err.println("Index ID: " + indexID);
			boolean ret = tree.isIndexValid();
			if (ret == false) System.err.println("Structure is INVALID!");
			tree.flush();
						
			tree.BuildInvertedList();
//			System.out.println("????IR?????!");
			//-------------------------------------------------------------------------------
			
			
			
			//---------------------------???query?????------------------------------------
			ArrayList<ArrayList<String>> querytemp = new ArrayList<ArrayList<String>>();
			
			for(int m = 0; m < querypoint.length; m++){//???query???????topic?????
				ArrayList<String> keytemp = new ArrayList<String>();
				String key = "";
				for(int n = 0; n<L; n++){//?????hash table??????????
					key = getKey(querypoint[m].topics, n, a, b, M, w);
					keytemp.add(key);
				}
				querytemp.add(keytemp);
			}

			//??????????????????
			ArrayList<ArrayList<String>> resultlist = new ArrayList<ArrayList<String>>();
			
			resultlist = permutation(querytemp);
//			System.out.println("?????????????????????MoSKQ!");
			//---------------------------------------------------------------------------
			long start = System.currentTimeMillis();//??????
			MyVisitor vis = new MyVisitor();
			
			if (querytype.equals("containment")) {/*??????????*/
				
			}
			else if (querytype.equals("intersection")) {
//					Region r = new Region(qf1, qf2);
//					tree.intersectionQuery(r, vis);
					// this will find all data that intersect with the query range.
			}
			else if(querytype.equals("pointLocation")) {/*??????*/
				
			}
			else if(querytype.equals("10NN")) {
					//Point p = new Point(qf1);
					//tree.nearestNeighborQuery(10, p, queryword, vis);
						// this will find the 10 nearest neighbors.
			}
			else if(querytype.equals("CSKQ_Type1")) {
//				List<IShape> p = new ArrayList<IShape>();//???????????
//				
//				for(int j=0;j<count;j++){
//					Point pp= new Point(qf1[j]);
//					p.add(pp);
//				}
//				tree.CollectiveQuery(p, queryword, vis);
			}
			else if(querytype.equals("CSKQ_Type2")) {
//				List<IShape> p = new ArrayList<IShape>();
//				
//				for(int j=0;j<count;j++){
//					Point pp= new Point(qf1[j]);
//					p.add(pp);
//				}						
//				tree.CollectiveQuery_Type2(p, queryword, vis);
				//this will find collective objects
			}
			else if(querytype.equals("CSKQ_Type3")) {
//				List<IShape> p = new ArrayList<IShape>();//???????????
//				
//				for(int j=0;j<count;j++){
//					Point pp= new Point(qf1[j]);
//					p.add(pp);
//				}
//				tree.CollectiveQuery_Type3(p, queryword, vis);
			}			
			else if(querytype.equals("Top-K")) {
				List<IRPoint> po = new ArrayList<IRPoint>();
				List<IRPoint> optimalpo = new ArrayList<IRPoint>();
				cost = Double.MAX_VALUE;
				int radius = 0;


				for(int k = 0;k<point.length;++k){
					System.out.println(point[k].m_pointID);
				}
				
				//upperbound = cost; lowerbound = 0.5*0.5*radius
//				while(cost >= 0.5*0.5*radius)
				{
//					radius += 10;
//					po = tree.getSubsetInRadius(radius,queryPoint,querycount);
//					cost = tree.getCost(po);
				}
				
				System.out.println("The final distance cost is:" + cost);
			}
			else if(querytype.equals("MoSKQ")) {
				/*
				 * ?????MoSKQ???????????CSKQ??????CollectiveQuery_temp
				 * */
				List<List<IRPoint>> po = new ArrayList<>();
				List<List<IRPoint>> poe = new ArrayList<>();
				cost = Double.MAX_VALUE;
				double costtemp = 0.0;
				Point p = new Point(qf1);
//				System.out.println("resultlist size:"+resultlist.size());
				for(int c=0; c<resultlist.size(); c++){
					po = tree.CollectiveQuery_temp(p, resultlist.get(c), vis);
//					for(int n=0;n<po.size();++n){
//						System.out.println(po.get(n).size());
//					}
//					System.out.println(po.size());
					permutationIRPoint(queryPoint,alpha,beta,po);
				}
				
				cost = ds(queryPoint,this.po,alpha,beta);
				System.out.println("The final distance cost is:" + cost);
			}
			else if(querytype.equals("IMOSKQ")) {
				int N = 3;
				
				Random r =new Random();
				
				
				List<List<IRPoint>> po = new ArrayList<>();
				List<List<IRPoint>> poe = new ArrayList<>();
				cost = Double.MAX_VALUE;
				double costtemp = 0.0;
				Point p = new Point(qf1);
//				System.out.println("resultlist size:"+resultlist.size());
				for(int c=0; c<N; c++){
					System.out.println("the "+ c +" start");
					po = tree.CollectiveQuery_temp(p, resultlist.get(r.nextInt(resultlist.size())), vis);
					printPo(po);
//					for(int n=0;n<po.size();++n){
//						System.out.println(po.get(n).size());
//					}
//					System.out.println(po.size());
					permutationIRPoint(queryPoint,alpha,beta,po);
				}
				cost = ds(queryPoint,this.po,alpha,beta);
				System.out.println("After random "+ N +" ,distance cost is:" + cost);
				//??
				List<IRPoint> po_temp = new ArrayList<IRPoint>();
				double min = 0;
				for(int k=0;k<this.po.size();++k){
					System.out.println("replace " + k);
					min = 0;
					po_temp.clear();
//					for(int v=0;v<this.po.size();v++){
//						po_temp.add(this.po.get(v));
//					}

					po_temp.addAll(this.po);
					for(int m=0;m<this.point.length;++m){
						if(!this.po.contains(point[m])) {
							po_temp.set(k, point[m]);
						}
						double newcost = ds(queryPoint,po_temp,alpha,beta);
						if(newcost - cost < min){

							min = newcost - cost;
							System.out.println("new min " + min);
							this.po.clear();
							this.po.addAll(po_temp);
//							for(int v=0;v<po_temp.size();v++){
//								this.po.add(po_temp.get(v));
//							}
						}
					}
					
					if(min >= 0){
						break;
					}else{
						if(k == this.po.size()-1)
							k=-1;
					}
					
				}
				cost = ds(queryPoint,this.po,alpha,beta);
				System.out.println("The final distance cost is:" + cost);
			}
			else {
				System.err.println("Unknown query type.");
				System.exit(-1);
			}

			indexIO += vis.m_indexIO;
			leafIO += vis.m_leafIO;
				// example of the Visitor pattern usage, for calculating how many nodes
				// were visited.
			System.out.println();
			line = lr.readLine();
		    
			//---------------------------------------------------------------------------------------------
		    long end = System.currentTimeMillis();
		    Thread.sleep(2);
		    System.err.println("Seconds: " + ((end - start) / 1000.0f));
			
		    diskfile.close();
		    
		}catch(Exception e)
		{
			e.printStackTrace();
		}
		}
	}
	
	private void printPo(List<List<IRPoint>> po){
		System.out.println("po size:"+po.size());
		for(int i = 0;i < po.size(); ++i){
			System.out.println(i + " size:"+po.get(i).size());
		}
		
	}
	
	private ArrayList<ArrayList<IRPoint>> permutationIRPoint(Point query,double alpha,double beta,List<List<IRPoint>> po){
		
		ArrayList<ArrayList<IRPoint>> resList = new ArrayList<ArrayList<IRPoint>>();
		ArrayList<IRPoint> arr = new ArrayList<IRPoint>();
		permutationIRPoint(query,alpha,beta,po, resList, 0, arr);//new char[inputList.size()] -> new ArrayList<String>() 
	    return resList;
	}

	private void permutationIRPoint(Point query,double alpha,double beta,List<List<IRPoint>> po, ArrayList<ArrayList<IRPoint>> resList, int ind, ArrayList<IRPoint> arr) {
	    if(ind == po.size()){
	    	double ds = ds2(query,arr,alpha,beta);
//		 	System.out.println(ds);
		 	if(ds<cost){
		 		cost = ds;
		 		this.po.clear();
		 		this.po.addAll(arr);
		 	}
//	        Object clone = arr.clone();
//			resList.add((ArrayList<IRPoint>) clone);
	        return;
	    }

	    for(IRPoint c: po.get(ind)){
	        arr.add(c);
	        permutationIRPoint(query,alpha,beta,po, resList, ind + 1, arr);
	        arr.remove(arr.size()-1);
	    }
	}

//	private ArrayList<String> GetOptimalCBS(ArrayList<ArrayList<String>> querytemp) {
//		ArrayList<String> ocbs = new ArrayList<String>();
//		ArrayList<String> tempocbs = new ArrayList<String>();
//		List<IRPoint> po = new ArrayList<IRPoint>();
//				
//		for(int i=0; i<querytemp.size(); ++i){
//			Random x = new Random();
//			
//			int temp = x.nextInt(3)+1;
//			
//			tempocbs.add(querytemp.get(i).get(temp));
//		}
//		
//		po = CollectiveQuery_temp(p, tempocbs, vis);
//		
//		return ocbs;
//	}

	public int getQuerycount() {
		return querycount;
	}

	public void setQuerycount(int querycount) {
		this.querycount = querycount;
	}

	private ArrayList<String> Replace(ArrayList<String> tempocbs,ArrayList<ArrayList<String>> querytemp) {
		Random x = new Random();
		ArrayList<String> tempocbs2 = new ArrayList<String>();
		
		int temp = x.nextInt(2);
		int temp2 = x.nextInt(3);
		
		tempocbs.set(temp, querytemp.get(temp).get(temp2));
		
		return tempocbs;
	}

	private ArrayList<ArrayList<IRPoint>> SubSet(List<IRPoint> po, int k) {

		ArrayList<ArrayList<IRPoint>> allset = new ArrayList<>();
		
		for(int i=1; i<=k; ++i){
			ArrayList<ArrayList<IRPoint>> resList = new ArrayList<>();
			resList = perm(po, k);
			allset.addAll(resList);
		}
		
		return allset;
	}

	private ArrayList<ArrayList<String>> permutation(ArrayList<ArrayList<String>> inputList){
		ArrayList<ArrayList<String>> resList = new ArrayList<ArrayList<String>>();
		ArrayList<String> arr = new ArrayList<String>();
		for(int i=0; i<inputList.size(); i++){
			arr.add("");
		}
	    permutationInt(inputList, resList, 0, arr);//new char[inputList.size()] -> new ArrayList<String>() 
	    return resList;
	}

	private void permutationInt(ArrayList<ArrayList<String>> inputList, ArrayList<ArrayList<String>> resList, int ind, ArrayList<String> arr) {
	    if(ind == inputList.size()){
	        Object clone = arr.clone();
			resList.add((ArrayList<String>) clone);
	        return;
	    }

	    for(String c: inputList.get(ind)){
	        arr.set(ind, c);
	        permutationInt(inputList, resList, ind + 1, arr);
	    }
	}
	
	private static ArrayList<ArrayList<IRPoint>> perm(List<IRPoint> inputList,int k ){
		ArrayList<ArrayList<IRPoint>> resList = new ArrayList<>();
		
	    ArrayList<IRPoint> arr = new ArrayList<>();
		permInt(inputList, resList, 0,0,arr,k );//new char[inputList.size()] -> new ArrayList<String>() 
	    return resList;
	}

	private static void permInt(List<IRPoint> inputList, List<ArrayList<IRPoint>> resList, int ind, int gnd, ArrayList<IRPoint> arr,int k) {
	    if(ind == k){
	    	ArrayList<IRPoint> clone = (ArrayList<IRPoint>) arr.clone();
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

class MyVisitor implements IVisitor {
	public int m_indexIO = 0;
  	public int m_leafIO = 0;
  
  	public void visitNode(final INode n)
  	{
  		if (n.isLeaf()) m_leafIO++;
  		else m_indexIO++;
  	}
  
  	public String visitData(final IData d)
  	{
  		String Keyword = new String(d.getData());
  		//System.out.println("??????id???" + d.getIdentifier());  		
  			// the ID of this data entry is an answer to the query. I will just print it to stdout.
  		return Keyword;
  	}
 }
