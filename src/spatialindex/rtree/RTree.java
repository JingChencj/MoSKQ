// Spatial Index Library
//
// Copyright (C) 2002  Navel Ltd.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License aint with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// Contact information:
//  Mailing address:
//    Marios Hadjieleftheriou
//    University of California, Riverside
//    Department of Computer Science
//    Surge Building, Room 310
//    Riverside, CA 92521
//
//  Email:
//    marioh@cs.ucr.edu

package spatialindex.rtree;

import java.util.*;
import java.util.function.Consumer;
import java.io.*;

import Test.IRPoint;
import Test.newTest;
import spatialindex.spatialindex.*;
import spatialindex.storagemanager.*;

public class RTree implements ISpatialIndex
{
	RWLock m_rwLock;

	IStorageManager m_pStorageManager;   //存储器管理

	int m_rootID;
	int m_headerID;

	int m_treeVariant;

	double m_fillFactor;

	int m_indexCapacity;

	int m_leafCapacity;

	int m_nearMinimumOverlapFactor;   //factor：因素     overlap：重叠
		// The R*-Tree 'p' constant（常量）, for calculating nearly minimum overlap cost.
		// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
		// for Points and Rectangles, Section 4.1]

	double m_splitDistributionFactor;   //split：分裂，分解
		// The R*-Tree 'm' constant, for calculating spliting distributions.
		// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method 
		// for Points and Rectangles, Section 4.2]

	double m_reinsertFactor;
		// The R*-Tree 'p' constant, for removing entries at reinserts.
		// [Beckmann, Kriegel, Schneider, Seeger 'The R*-tree: An efficient and Robust Access Method
		//  for Points and Rectangles, Section 4.3]

	int m_dimension;

	Region m_infiniteRegion;   //infinite：[数]无穷大

	Statistics m_stats;
	
	Node n;
	
	ArrayList m_writeNodeCommands = new ArrayList();
	ArrayList m_readNodeCommands = new ArrayList();
	ArrayList m_deleteNodeCommands = new ArrayList();

	public RTree(PropertySet ps, IStorageManager sm)
	{
		m_rwLock = new RWLock();
		m_pStorageManager = sm;
		m_rootID = IStorageManager.NewPage;//IStorageManager的静态属性，默认为-1
		m_headerID = IStorageManager.NewPage;
		m_treeVariant = SpatialIndex.RtreeVariantRstar;
		m_fillFactor = 0.7f;
		m_indexCapacity = 2;
		m_leafCapacity = 2;
		m_nearMinimumOverlapFactor = 32;
		m_splitDistributionFactor = 0.4f;
		m_reinsertFactor = 0.3f;
		m_dimension = 2;

		m_infiniteRegion = new Region();
		m_stats = new Statistics();

		Object var = ps.getProperty("IndexIdentifier");
		if (var != null)
		{
			if (! (var instanceof Integer)) throw new IllegalArgumentException("Property IndexIdentifier must an Integer");
			m_headerID = ((Integer) var).intValue();
			try
			{
				initOld(ps);
			}
			catch (IOException e)
			{
				System.err.println(e);
				throw new IllegalStateException("initOld failed with IOException");
			}
		}
		else
		{
			try
			{
				initNew(ps);
			}
			catch (IOException e)
			{
				System.err.println(e);
				throw new IllegalStateException("initNew failed with IOException");
			}
			Integer i = new Integer(m_headerID);
			ps.setProperty("IndexIdentifier", i);
		}
	}
	
	public void BuildInvertedList()
	{
		n = readNode(m_rootID);
		n.BuildInvertedFile();
	}
	
	//
	// ISpatialIndex interface
	//

	public void insertData(final byte[] data, final IShape shape, int id)
	{
		if (shape.getDimension() != m_dimension) throw new IllegalArgumentException("insertData: IShape dimensionality is incorrect.");

		m_rwLock.write_lock();

		try
		{
			Region mbr = shape.getMBR();

			byte[] buffer = null;

			if (data != null && data.length > 0)
			{
				buffer = new byte[data.length];
				System.arraycopy(data, 0, buffer, 0, data.length);
			}

			insertData_impl(buffer, mbr, id);
				// the buffer is stored in the tree. Do not delete here.		
		}
		finally
		{
			m_rwLock.write_unlock();
		}
	}

	public boolean deleteData(final IShape shape, int id)
	{
		if (shape.getDimension() != m_dimension) throw new IllegalArgumentException("deleteData: IShape dimensionality is incorrect.");

		m_rwLock.write_lock();

		try
		{
			Region mbr = shape.getMBR();
			return deleteData_impl(mbr, id);
		}
		finally
		{
			m_rwLock.write_unlock();
		}
	}

	public void containmentQuery(final IShape query, final IVisitor v)
	{
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("containmentQuery: IShape dimensionality is incorrect.");
		rangeQuery(SpatialIndex.ContainmentQuery, query, v);
	}

	public void intersectionQuery(final IShape query, final IVisitor v)
	{
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("intersectionQuery: IShape dimensionality is incorrect.");
		rangeQuery(SpatialIndex.IntersectionQuery, query, v);
	}

	public void pointLocationQuery(final IShape query, final IVisitor v)
	{
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("pointLocationQuery: IShape dimensionality is incorrect.");
		
		Region r = null;
		if (query instanceof Point)
		{
			r = new Region((Point) query, (Point) query);
		}
		else if (query instanceof Region)
		{
			r = (Region) query;
		}
		else
		{
			throw new IllegalArgumentException("pointLocationQuery: IShape can be Point or Region only.");
		}

		rangeQuery(SpatialIndex.IntersectionQuery, r, v);
	}
	
	public List<IRPoint> nearestNeighborQuery(int k, final IShape query, String keyword, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("nearestNeighborQuery: IShape dimensionality is incorrect.");
		NNComparator nnc = new NNComparator();
		po = nearestNeighborQuery(k, query, keyword, v, nnc);
		return po;
	}

	public List<IRPoint> nearestNeighborQuery(int k, final IShape query, String keyword, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();

		try
		{
			// I need a priority queue here. It turns out that TreeSet sorts unique keys only and since I am
			// sorting according to distances, it is not assured that all distances will be unique. TreeMap
			// also sorts unique keys. Thus, I am simulating a priority queue using an ArrayList and binarySearch.
			ArrayList queue = new ArrayList();
			ArrayList Result_queue = new ArrayList();

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));//加入根节点，距离为0

			int count = 0;
			double knearest = 0.0;

			while (queue.size() != 0)
			{
				NNEntry first = (NNEntry) queue.remove(0);

				if (first.m_pEntry instanceof Node)//如果first.m_pEntry是Node的一个实例
				{
					n = (Node) first.m_pEntry;
					v.visitNode((INode) n);//访问节点n

					for (int cChild = 0; cChild < n.m_children; cChild++)
					{
						IEntry e;

						if (n.m_level == 0)
						{
							e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
							//如果此节点是叶子节点，则依次读出叶子节点当中的对象的值
							String word = new String(n.m_pData[cChild]);//如果此节点是叶节点，则找出此实例当中的关键字						
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(query, e));
							//如果e在这里代表着对象，则表示查询点与此对象之间的最小距离；如果e在这里代表着节点，则表示查询点与此节点的MBR之间的距离
							
							// Why don't I use a TreeSet here? See comment above...
							if(word.equals(keyword)){
								Result_queue.add(nnc.getMinimumDistance(query, e));//将关键字相同的点之间的距离存到数组当中
								
								int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
								if (loc >= 0) queue.add(loc, e2);//增加条件，使得查询到的符合条件的点的关键字与查询的关键字相同
								else queue.add((-loc - 1), e2);									
							}
						}
						else
						{
							e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(query, e));
							
							int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
							if (loc >= 0) queue.add(loc, e2);
							else queue.add((-loc - 1), e2);
						}
					}
				}
				else
				{
					// report all nearest neighbors with equal furthest distances.
					// (neighbors can be more than k, if many happen to have the same
					//  furthest distance).
					if (count >= k && first.m_minDist > knearest) break;

					v.visitData((IData) first.m_pEntry);

					IRPoint temp = new IRPoint();
					temp.m_pointID = first.m_pEntry.getIdentifier();
					temp.Keyword = keyword;
					temp.distance = first.m_minDist;

					po.add(temp);
					m_stats.m_queryResults++;
					count++;
					knearest = first.m_minDist;
				}
			}
			
			Collections.sort(Result_queue);
			System.out.println("Key Distances:" + Result_queue);
		}
		finally
		{
			m_rwLock.read_unlock();
		}
		return po;
	}
	
	//-----------------------------------------------find top-k----------------------------------------------------------------------------
	
	public List<IRPoint> TOPK(int k, final IShape query, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("TOPK: IShape dimensionality is incorrect.");
		NNComparator nnc = new NNComparator();
		po = TOPK(k, query, v, nnc);
		return po;
	}

	public List<IRPoint> TOPK(int k, final IShape query, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();

		try
		{
			ArrayList queue = new ArrayList();
			ArrayList Result_queue = new ArrayList();

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));//加入根节点，距离为0

			int count = 0;
			double knearest = 0.0;

			while (queue.size() != 0)
			{
				NNEntry first = (NNEntry) queue.remove(0);

				if (first.m_pEntry instanceof Node)//如果first.m_pEntry是Node的一个实例
				{
					n = (Node) first.m_pEntry;
					v.visitNode((INode) n);//访问节点n

					for (int cChild = 0; cChild < n.m_children; cChild++)
					{
						IEntry e;

						if (n.m_level == 0)
						{
							e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
							//如果此节点是叶子节点，则依次读出叶子节点当中的对象的值
							String word = new String(n.m_pData[cChild]);//如果此节点是叶节点，则找出此实例当中的关键字						
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(query, e));
							//如果e在这里代表着对象，则表示查询点与此对象之间的最小距离；如果e在这里代表着节点，则表示查询点与此节点的MBR之间的距离
							
							// Why don't I use a TreeSet here? See comment above...
							Result_queue.add(nnc.getMinimumDistance(query, e));//将关键字相同的点之间的距离存到数组当中
							
							int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
							if (loc >= 0) queue.add(loc, e2);//增加条件，使得查询到的符合条件的点的关键字与查询的关键字相同
							else queue.add((-loc - 1), e2);									
						}
						else
						{
							e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(query, e));
							
							int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
							if (loc >= 0) queue.add(loc, e2);
							else queue.add((-loc - 1), e2);
						}
					}
				}
				else
				{
					// report all nearest neighbors with equal furthest distances.
					// (neighbors can be more than k, if many happen to have the same
					//  furthest distance).
					if (count >= k && first.m_minDist > knearest) break;

					v.visitData((IData) first.m_pEntry);

					IRPoint temp = new IRPoint();
					temp.m_pointID = first.m_pEntry.getIdentifier();
					//temp.Keyword = keyword;
					temp.distance = first.m_minDist;

					po.add(temp);
					m_stats.m_queryResults++;
					count++;
					knearest = first.m_minDist;
				}
			}
			
			Collections.sort(Result_queue);
			System.out.println("Key Distances:" + Result_queue);
		}
		finally
		{
			m_rwLock.read_unlock();
		}
		return po;
	}

	public double getCost(List<IRPoint> po){
		double a=0.5;
		//找出结果集中距离查询最远的一个对象，并找出此对象的关键字
		double maxDistance = 0.0;
		double temp = 0.0;

		//po[i].distance即节点与查询点之间的距离
		
		double qx = newTest.coord[0];
		double qy = newTest.coord[1];
		
		for(int i=0;i<po.size();i++){

			double x = po.get(i).m_pCoordinate[0];//double x = po.get(i).getLatitude();
			double y = po.get(i).m_pCoordinate[1];
			
			temp = Math.sqrt(Math.abs((x-qx)*(x-qx)-(y-qy)*(y-qy)));
			if(temp > maxDistance)
			{
				maxDistance = temp;
			}
		}

		List<Double> Result_queue = new ArrayList<Double>();
		
		double maxDisObj = 0.0;
		double temp2 = 0.0;
		if(po.size() > 1)
		{
			for(int i=0;i<po.size()-1;i++){
				for(int j=i+1;j<po.size();j++){
					double x1 = po.get(i).getLatitude();
					double y1 = po.get(i).getLongitude();
					
					double x2 = po.get(j).getLatitude();
					double y2 = po.get(j).getLongitude();
					
					double distobj = Math.sqrt(Math.abs((x1-x2)*(x1-x2)-(y1-y2)*(y1-y2)));
					
					if(temp2 > maxDisObj)
					{
						maxDisObj = temp2;
					}
//					Result_queue.add(distobj);//求出对象间的距离并存储起来
				}
			}
		}

		double Cost = a*maxDistance + (1-a)*maxDisObj;//最终的Cost是查询与对象之家的最大距离与对象之间的最大距离的和
		
		return Cost;
	}
	
	public List<IRPoint> CollectiveQuery(final List<IShape> query, List<String> keyword, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();//结果集
		
		for(int i=0;i<query.size();i++)//确定查询集合中的每个查询点都是二维的
		if (query.get(i).getDimension() != m_dimension) throw new IllegalArgumentException("nearestNeighborQuery: IShape dimensionality is incorrect.");
		
		NNComparator nnc = new NNComparator();
		po = CollectiveQuery(query, keyword, v, nnc);//调用下面的方法查找结果集
		return po;
	}

	public List<IRPoint> CollectiveQuery(final List<IShape> query, List<String> keyword, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();//作为最后的查询结果集

		try
		{
			ArrayList queue = new ArrayList();//相当于Type2Appro1中的U

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));//在这里的NNEntry中包含了节点和距离，与U.Enqueue(irTree.root,0)是一样的效果

			ArrayList Result_queue = new ArrayList();//用于存储查询与对象点之间的距离
			
			//keyword在这里是一个未被覆盖的关键字集合
			while (queue.size() != 0)//while U is not empty do
			{
				NNEntry first = (NNEntry) queue.remove(0);//e = U.Dequeue();

				if (first.m_pEntry instanceof Node)//如果first.m_pEntry是一个节点
				{
					n = (Node) first.m_pEntry;
					v.visitNode((INode) n);//访问节点n

					for (int cChild = 0; cChild < n.m_children; cChild++)//(foreach entry e' in node e do)
					{
						if(n.m_pData != null)
						{
							IEntry e;
							//在加入节点的过程中对节点的位置做排序
							if (n.m_level == 0) {									
								String temp = new String(n.m_pData[cChild]);
								
								if(keyword.contains(temp))
								{										
									e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
								
									double d = 0;
									for(int i=0;i<query.size();i++){
										if(d < nnc.getMinimumDistance(query.get(i), e))
											d = nnc.getMinimumDistance(query.get(i), e);//d是第i个查询点与此叶子节点的所有最小距离中的最大距离
									}
									NNEntry e2 = new NNEntry(e, d);
	
									queue.add(e2);
								}								
							} else {
								e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点

								double d = 0;
								double d2 = 0;
								for(int i=0;i<query.size();i++){
									if(d < nnc.getMinimumDistance(query.get(i), e))
										d = nnc.getMinimumDistance(query.get(i), e);
									if(d2<d)
										d2=d;
								}//找出距离查询点的最大距离
								NNEntry e2 = new NNEntry(e, d);
								System.out.println("距离为："+d2);

								queue.add(e2);
							}
						}
					}
				}
				else
				{
					IEntry e = first.m_pEntry;
					if(keyword.size() == 0)//当所有的关键字都已经被查询到时，结束查询
						break;

//					System.out.println("移除前的关键字集合为："+ keyword);
					
					IRPoint temp = new IRPoint();
					temp.m_pointID = first.m_pEntry.getIdentifier();
					temp.Keyword = v.visitData((IData) first.m_pEntry);
					
//					System.out.println("temp.Keyword的值为："+ temp.Keyword);


					if(keyword.contains(temp.Keyword))//在此步找到结果集
					{
						double d = 0;
						double d2 = 0;
						for(int i=0;i<query.size();i++){
							if(d < nnc.getMinimumDistance(query.get(i), e))//一个结果对象与全部查询点之间的
								d = nnc.getMinimumDistance(query.get(i), e);
							if(d2<d)
								d2=d;
						}

						temp.distance = d2;
//						System.out.println("true 在此步去除一个关键字");
						po.add(temp);
						keyword.remove(temp.Keyword);
					}
					else
					{
						continue;
					}
					m_stats.m_queryResults++;
					
					System.out.println("关键字集合为："+ keyword);
					System.out.println();
				}
			}
			System.out.println("查询结束");
		}
		finally
		{
			m_rwLock.read_unlock();
		}

		System.out.println("结果集的长度为："+ po.size());
		double mincost = getCost(po);
		System.out.println("最终结果集对应的最小距离代价为："+ mincost);
		return po;
	}

	public List<IRPoint> CollectiveQuery_Type2(final List<IShape> query, List<String> keyword, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		
		for(int i=0;i<query.size();i++)
		if (query.get(i).getDimension() != m_dimension) throw new IllegalArgumentException("nearestNeighborQuery: IShape dimensionality is incorrect.");
		
		NNComparator nnc = new NNComparator();
		po = CollectiveQuery_Type2(query, keyword, v, nnc);
		return po;
	}
	
	public List<IRPoint> CollectiveQuery_Type2(final List<IShape> query, List<String> keyword, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();

		try
		{
			ArrayList queue = new ArrayList();

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));
			
			po = CollectiveQuery(query,keyword,v);//类型1得到的结果

			double Cost = getCost(po);

			//找出结果集当中最远的节点中包含的关键字

			while (queue.size() != 0)
			{
				NNEntry first = (NNEntry) queue.remove(0);

				if (first.m_pEntry instanceof Node)//if e is not an object then
				{
					n = (Node) first.m_pEntry;
					v.visitNode((INode) n);//访问节点n
					
					//确保结果集中的点包含的distance是存储对象到多个查询点之间的最小距离
					//在这里找出最小距离中的最大距离
					double maxDistance = po.get(0).distance;
					for(int i=1;i<po.size();i++){
						if(maxDistance < po.get(i).distance)
							maxDistance = po.get(i).distance;
					}
					
					if(maxDistance >= Cost)
						break;

					for (int cChild = 0; cChild < n.m_children; cChild++)//(foreach entry e' in node e do)
					{
						if(n.m_pData != null)
						{
							IEntry e;
							//在加入节点的过程中对节点的位置做排序
							if (n.m_level == 0) {									
								String temp = new String(n.m_pData[cChild]);
								
								if(keyword.contains(temp))
								{										
									e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
								
									double d = 0;
									for(int i=0;i<query.size();i++){
										if(d < nnc.getMinimumDistance(query.get(i), e))
											d = nnc.getMinimumDistance(query.get(i), e);
									}
									NNEntry e2 = new NNEntry(e, d);
									//如果e在这里代表着对象，则表示查询点与此对象之间的最小距离；如果e在这里代表着节点，则表示查询点与此节点的MBR之间的距离
	
									queue.add(e2);
								}
								
							} else {
								e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点

								double d = 0;
								for(int i=0;i<query.size();i++){
									if(d < nnc.getMinimumDistance(query.get(i), e))
										d = nnc.getMinimumDistance(query.get(i), e);
								}
								NNEntry e2 = new NNEntry(e, d);

								queue.add(e2);
							}
						}
					}
				}
				else
				{
					//求出查询集合与对象之间的距离的最大值
					IShape s = first.m_pEntry.getShape();
					double Dist=0;
					for(int i=0;i<query.size();i++)
					{
						if(Dist < query.get(i).getMinimumDistance(s))
						Dist = query.get(i).getMinimumDistance(s);
					}
					
					if(Dist >= Cost)
						break;
					
					IRPoint qe = new IRPoint();
					List<IShape> querye = query;
					
					List<IRPoint> poe = CollectiveQuery(querye,keyword,v);;
					
					double Coste = getCost(poe);
					
					if(Coste<Cost){//如果找到了所需Cost更少的一组对象，那么更新Cost值和最优结果集。
						Cost = Coste;
						po = poe;
					}
				}
			}
		}
		finally
		{
			m_rwLock.read_unlock();
		}

		double mincost = getCost(po);
		System.out.println("最终结果集对应的最小距离代价为："+ mincost);
		return po;
	}
	
	
	/**
	 * 返回所有的objects,因为多个bucket，所以是多组object
	 */
	public List<List<IRPoint>> CollectiveQuery_temp(final IShape query, List<String> buckids, final IVisitor v)
	{
		List<List<IRPoint>> po = null;//结果集
//		System.out.println("buckids:");
//		for(int i=0;i<buckids.size();++i){
//			System.out.print(buckids.get(i));
//		}
//		System.out.println();
		if (query.getDimension() != m_dimension) throw new IllegalArgumentException("nearestNeighborQuery: IShape dimensionality is incorrect.");
		
		NNComparator nnc = new NNComparator();
		po = CollectiveQuery_temp(query, buckids, v, nnc);//调用下面的方法查找结果集
		return po;
	}
	
	class Entry{
		public int id;
		public boolean isLeaf;
		
		public Entry(int id,boolean isLeaf) {
			this.id = id;
			this.isLeaf = isLeaf;
		}
	}
	
	public List<List<IRPoint>> CollectiveQuery_temp(final IShape query, List<String> buckids, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<List<IRPoint>> po = new ArrayList<>();//作为最后的查询结果集
		HashMap<String, HashSet> hash = new HashMap<String, HashSet>();//key:buckid,value:objectids 	
		try
		{
			ArrayList queue = new ArrayList();//相当于Type2Appro1中的U
			queue.add(new Entry(n.m_identifier, n.isLeaf()));
			while(!queue.isEmpty()){
				Entry cn = (Entry) queue.remove(0);
//				System.out.println("current node id:" + cn.id);
				
				if(!cn.isLeaf){//非叶子节点
//					System.out.println("current node is not leaf node");
					HashMap<String,ArrayList> invert = newTest.inverted.get(cn.id);
					Iterator<String> iterator = invert.keySet().iterator();
//					System.out.print("inverted keys:" );
//					while(iterator.hasNext()){
//						String next = iterator.next();
//						System.out.print("" + next+",");
//					}
//					System.out.println();
					for(int i=0;i<buckids.size();++i){
//						System.out.println("to find buckid:"+buckids.get(i));
						if(invert.containsKey(buckids.get(i))){//有这个buckid
//							System.out.println("found");
							List l = invert.get(buckids.get(i));
							for(int k=0;k<l.size();++k){
								int id = (int) l.get(k);
								Node n2 = readNode(id);
								queue.add(new Entry(n2.m_identifier, n2.isLeaf()));
//								System.out.println("add node");
							}
						}
					}
				}else{//叶子结点
					HashMap<String,ArrayList> invert = newTest.inverted.get(cn.id);//存的是buckidid和objectids
					
					for(int i=0;i<buckids.size();++i){
						if(invert.containsKey(buckids.get(i))){//有这个buckid
							String next = buckids.get(i);
							if(!hash.containsKey(next)){
								hash.put(next, new HashSet<Integer>());
							}
							hash.get(next).addAll(invert.get(next));
					
						}
					}
					
					
				}
				
				
			}
			
//			queue.add(new NNEntry(n, 0.0));//在这里的NNEntry中包含了节点和距离，与U.Enqueue(irTree.root,0)是一样的效果

//			ArrayList Result_queue = new ArrayList();//用于存储查询与对象点之间的距离

			double knearest = 0.0;
			
			Iterator<String> iterator = hash.keySet().iterator();
			while(iterator.hasNext()){
				String next = iterator.next();//bucketid
				HashSet hashSet = hash.get(next);//objectids
				ArrayList<IRPoint> objects = new ArrayList<IRPoint>();
				Iterator iterator2 = hashSet.iterator();
				while(iterator2.hasNext()){
					int next2 = (int) iterator2.next();
					objects.add(newTest.map.get(next2));
				}
				po.add(objects);
			}
			
		}
		finally
		{
			m_rwLock.read_unlock();
		}

		//System.out.println("结果集的长度为："+ po.size());
		return po;
	}
	
	public List<IRPoint> CollectiveQuery_Type3(final List<IShape> query, List<String> keyword, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		
		for(int i=0;i<query.size();i++)
		if (query.get(i).getDimension() != m_dimension) throw new IllegalArgumentException("nearestNeighborQuery: IShape dimensionality is incorrect.");
		
		NNComparator nnc = new NNComparator();
		po = CollectiveQuery_Type3(query, keyword, v, nnc);
		return po;
	}
	
	public List<IRPoint> CollectiveQuery_Type3(final List<IShape> query, List<String> keyword, final IVisitor v, final INearestNeighborComparator nnc)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		List<IRPoint> poe = new ArrayList<IRPoint>();

		try
		{
			ArrayList queue = new ArrayList();

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));
			
			List<String> keyword2 = new ArrayList<String>(keyword);
//			po = CollectiveQuery_temp(query.get(0),keyword2,v);

			double Cost = getCost(po);
			double Coste = Cost;

			if (query.size() > 0)
			{
				for(int i=1;i<query.size();i++)
				{
					List<String> keyword3 = new ArrayList<String>(keyword);
//					poe = CollectiveQuery_temp(query.get(i),keyword3,v);
					Coste = getCost(poe);
					
					if(Coste < Cost)
					{
						po = poe;
						Cost = Coste;
					}
				}
			}
		}
		finally
		{
			m_rwLock.read_unlock();
		}

		double mincost = getCost(po);
		System.out.println("最终结果集对应的最小距离代价为："+ mincost);
		return po;
	}
	
	//-------------------------------------MoSKQ-----------------------------------------------------------
	
	public List<IRPoint> MoSKQ(final ArrayList<String> query, final IVisitor v)
	{
		List<IRPoint> po = new ArrayList<IRPoint>();
		
		NNComparator nnc = new NNComparator();
//		System.out.println("进入子查询MoSKQ!");
		po = MoSKQ(query, v, nnc);
//		System.out.println("结束子查询MoSKQ!");
		return po;
	}
	
	public List<IRPoint> MoSKQ(final ArrayList<String> query, final IVisitor v, final INearestNeighborComparator nnc)
	{
		ArrayList<IRPoint> po = new ArrayList<IRPoint>();
		List<IRPoint> poe = new ArrayList<IRPoint>();

		try
		{
//			long time1 = System.currentTimeMillis();
			ArrayList queue = new ArrayList();

			ArrayList<ArrayList<Integer>> queryList = new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> resulttemp = new ArrayList<ArrayList<Integer>>();
			ArrayList<Integer> temp = new ArrayList<Integer>();

			for(int count=0 ; count<query.size(); count++){
				temp = n.m_invertedList.get(query.get(count));//temp代表一个桶内的一组id
				queryList.add(temp);//将querylist看做一个二维数组，每一行存储一组id，一共有K列，K表示关键字
			}
//			System.gc();
//			System.out.println("进入MoSKQ函数中的递归！");
			//在这里的递归是通过传入query对应的二维的一个桶号的数组，找出全部的结果集的桶号的组合可能
			resulttemp = permutation(queryList);
//			System.out.println("结束MoSKQ函数中的递归！");
			
			double cost = Double.MAX_VALUE;
			double costtemp = 0.0;
//			System.out.println("结果集的大小为:" + resulttemp.size());
			for(int itemp=0 ; itemp<resulttemp.size(); itemp++){
				for(int jtemp=0 ; jtemp<resulttemp.get(itemp).size() ; jtemp++){
					//已知一个空间文本对象的ID，如何获取其经纬度，并计算与查询点之间的距离
					po.add(newTest.map.get(resulttemp.get(itemp).get(jtemp)));
				}
//				System.gc();
//				System.out.println("计算结果集的距离代价！"+itemp);
				costtemp = getCost(po);
				if(cost > costtemp)
				{
					cost = costtemp;
					Object clone = po.clone();
					poe = (List<IRPoint>) clone;
					
				}
				po.clear();
			}
//			System.gc();
//			System.out.println("结束MoSKQ算法！");
//			long time2 = System.currentTimeMillis();
//			System.out.println("MOSKQ执行一次的时间" + (time2-time1));
		}
		finally
		{
			m_rwLock.read_unlock();
		}

		return poe;
	}
	
	
	public static ArrayList<ArrayList<Integer>> permutation(ArrayList<ArrayList<Integer>> inputList){
		ArrayList<ArrayList<Integer>> resList = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> arr = new ArrayList<Integer>();
		for(int i=0; i<inputList.size(); i++){
			arr.add(0);
		}
	    permutationInt(inputList, resList, 0, arr);//new char[inputList.size()] -> new ArrayList<String>() 
	    return resList;
	}

	public static void permutationInt(ArrayList<ArrayList<Integer>> inputList, ArrayList<ArrayList<Integer>> resList, int ind, ArrayList<Integer> arr) {
	    if(ind == inputList.size()){
	        Object clone = arr.clone();
			resList.add((ArrayList<Integer>) clone);
	        return;
	    }

	    for(int c: inputList.get(ind)){
	        arr.set(ind, c);
	        permutationInt(inputList, resList, ind + 1, arr);
	    }
	}
	
	//-----------------------------------------------------------------------------------------------------

	public void queryStrategy(final IQueryStrategy qs)
	{
		m_rwLock.read_lock();

		Integer next = new Integer(m_rootID);

		try
		{
			while (true)
			{
				Node n = readNode(next.intValue());
				Boolean hasNext = new Boolean(false);
				qs.getNextEntry(n, next, hasNext);
				if (hasNext.booleanValue() == false) break;
			}
		}
		finally
		{
			m_rwLock.read_unlock();
		}
	}

	public PropertySet getIndexProperties()
	{
		PropertySet pRet = new PropertySet();

		// dimension
		pRet.setProperty("Dimension", new Integer(m_dimension));

		// index capacity
		pRet.setProperty("IndexCapacity", new Integer(m_indexCapacity));

		// leaf capacity
		pRet.setProperty("LeafCapacity", new Integer(m_leafCapacity));

		// R-tree variant
		pRet.setProperty("TreeVariant", new Integer(m_treeVariant));

		// fill factor
		pRet.setProperty("FillFactor", new Double(m_fillFactor));

		// near minimum overlap factor
		pRet.setProperty("NearMinimumOverlapFactor", new Integer(m_nearMinimumOverlapFactor));

		// split distribution factor
		pRet.setProperty("SplitDistributionFactor", new Double(m_splitDistributionFactor));

		// reinsert factor
		pRet.setProperty("ReinsertFactor", new Double(m_reinsertFactor));

		return pRet;
	}

	public void addWriteNodeCommand(INodeCommand nc)
	{
		m_writeNodeCommands.add(nc);
	}

	public void addReadNodeCommand(INodeCommand nc)
	{
		m_readNodeCommands.add(nc);
	}

	public void addDeleteNodeCommand(INodeCommand nc)
	{
		m_deleteNodeCommands.add(nc);
	}

	public boolean isIndexValid()
	{
		boolean ret = true;
		Stack st = new Stack();
		Node root = readNode(m_rootID);

		if (root.m_level != m_stats.m_treeHeight - 1)
		{
			System.err.println("Invalid tree height");
			return false;
		}

		HashMap nodesInLevel = new HashMap();
		nodesInLevel.put(new Integer(root.m_level), new Integer(1));

		ValidateEntry e = new ValidateEntry(root.m_nodeMBR, root);
		st.push(e);

		while (! st.empty())
		{
			e = (ValidateEntry) st.pop();

			Region tmpRegion = (Region) m_infiniteRegion.clone();

			for (int cDim = 0; cDim < m_dimension; cDim++)
			{
				tmpRegion.m_pLow[cDim] = Double.POSITIVE_INFINITY;
				tmpRegion.m_pHigh[cDim] = Double.NEGATIVE_INFINITY;

				for (int cChild = 0; cChild < e.m_pNode.m_children; cChild++)
				{
					tmpRegion.m_pLow[cDim] = Math.min(tmpRegion.m_pLow[cDim], e.m_pNode.m_pMBR[cChild].m_pLow[cDim]);
					tmpRegion.m_pHigh[cDim] = Math.max(tmpRegion.m_pHigh[cDim], e.m_pNode.m_pMBR[cChild].m_pHigh[cDim]);
				}
			}

			if (! (tmpRegion.equals(e.m_pNode.m_nodeMBR)))
			{
				System.err.println("Invalid parent information");
				ret = false;
			}
			else if (! (tmpRegion.equals(e.m_parentMBR)))
			{
				System.err.println("Error in parent");
				ret = false;
			}

			if (e.m_pNode.m_level != 0)
			{
				for (int cChild = 0; cChild < e.m_pNode.m_children; cChild++)
				{
					ValidateEntry tmpEntry = new ValidateEntry(e.m_pNode.m_pMBR[cChild], readNode(e.m_pNode.m_pIdentifier[cChild]));

					if (! nodesInLevel.containsKey(new Integer(tmpEntry.m_pNode.m_level)))
					{
						nodesInLevel.put(new Integer(tmpEntry.m_pNode.m_level), new Integer(1));
					}
					else
					{
						int i = ((Integer) nodesInLevel.get(new Integer(tmpEntry.m_pNode.m_level))).intValue();
						nodesInLevel.put(new Integer(tmpEntry.m_pNode.m_level), new Integer(i + 1));
					}

					st.push(tmpEntry);
				}
			}
		}

		int nodes = 0;
		for (int cLevel = 0; cLevel < m_stats.m_treeHeight; cLevel++)
		{
			int i1 = ((Integer) nodesInLevel.get(new Integer(cLevel))).intValue();
			int i2 = ((Integer) m_stats.m_nodesInLevel.get(cLevel)).intValue();
			if (i1 != i2)
			{
				System.err.println("Invalid nodesInLevel information");
				ret = false;
			}

			nodes += i2;
		}

		if (nodes != m_stats.m_nodes)
		{
			System.err.println("Invalid number of nodes information");
			ret = false;
		}

		return ret;
	}

	public IStatistics getStatistics()
	{
		return (IStatistics) m_stats.clone();
	}

	public void flush() throws IllegalStateException
	{
		try
		{
			storeHeader();
			m_pStorageManager.flush();
		}
		catch (IOException e)
		{
			System.err.println(e);
			throw new IllegalStateException("flush failed with IOException");
		}
	}

	//
	// Internals
	//

	private void initNew(PropertySet ps) throws IOException
	{
		Object var;

		// tree variant.
		var = ps.getProperty("TreeVariant");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i != SpatialIndex.RtreeVariantLinear &&  i != SpatialIndex.RtreeVariantQuadratic && i != SpatialIndex.RtreeVariantRstar)
					throw new IllegalArgumentException("Property TreeVariant not a valid variant");
				m_treeVariant = i;
			}
			else
			{
				throw new IllegalArgumentException("Property TreeVariant must be an Integer");
			}
		}

		// fill factor.
		var = ps.getProperty("FillFactor");
		if (var != null)
		{
			if (var instanceof Double)
			{
				double f = ((Double) var).doubleValue();
				if (f <= 0.0f || f >= 1.0f)
					throw new IllegalArgumentException("Property FillFactor must be in (0.0, 1.0)");
				m_fillFactor = f;
			}
			else
			{
				throw new IllegalArgumentException("Property FillFactor must be a Double");
			}
		}

		// index capacity.
		var = ps.getProperty("IndexCapacity");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i < 3) throw new IllegalArgumentException("Property IndexCapacity must be >= 3");
				m_indexCapacity = i;
			}
			else
			{
				throw new IllegalArgumentException("Property IndexCapacity must be an Integer");
			}
		}

		// leaf capacity.
		var = ps.getProperty("LeafCapacity");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i < 3) throw new IllegalArgumentException("Property LeafCapacity must be >= 3");
				m_leafCapacity = i;
			}
			else
			{
				throw new IllegalArgumentException("Property LeafCapacity must be an Integer");
			}
		}

		// near minimum overlap factor.
		var = ps.getProperty("NearMinimumOverlapFactor");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i < 1 || i > m_indexCapacity || i > m_leafCapacity)
					throw new IllegalArgumentException("Property NearMinimumOverlapFactor must be less than both index and leaf capacities");
			m_nearMinimumOverlapFactor = i;
			}
			else
			{
				throw new IllegalArgumentException("Property NearMinimumOverlapFactor must be an Integer");
			}
		}

		// split distribution factor.
		var = ps.getProperty("SplitDistributionFactor");
		if (var != null)
		{
			if (var instanceof Double)
			{
				double f = ((Double) var).doubleValue();
				if (f <= 0.0f || f >= 1.0f)
					throw new IllegalArgumentException("Property SplitDistributionFactor must be in (0.0, 1.0)");
				m_splitDistributionFactor = f;
			}
			else
			{
				throw new IllegalArgumentException("Property SplitDistriburionFactor must be a Double");
			}
		}

		// reinsert factor.
		var = ps.getProperty("ReinsertFactor");
		if (var != null)
		{
			if (var instanceof Double)
			{
				double f = ((Double) var).doubleValue();
				if (f <= 0.0f || f >= 1.0f)
					throw new IllegalArgumentException("Property ReinsertFactor must be in (0.0, 1.0)");
				m_reinsertFactor = f;
			}
			else
			{
				throw new IllegalArgumentException("Property ReinsertFactor must be a Double");
			}
		}

		// dimension
		var = ps.getProperty("Dimension");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i <= 1) throw new IllegalArgumentException("Property Dimension must be >= 1");
				m_dimension = i;
			}
			else
			{
				throw new IllegalArgumentException("Property Dimension must be an Integer");
			}
		}

		m_infiniteRegion.m_pLow = new double[m_dimension];
		m_infiniteRegion.m_pHigh = new double[m_dimension];

		for (int cDim = 0; cDim < m_dimension; cDim++)
		{
			m_infiniteRegion.m_pLow[cDim] = Double.POSITIVE_INFINITY;
			m_infiniteRegion.m_pHigh[cDim] = Double.NEGATIVE_INFINITY;
		}

		m_stats.m_treeHeight = 1;
		m_stats.m_nodesInLevel.add(new Integer(0));

		Leaf root = new Leaf(this, -1);
		m_rootID = writeNode(root);

		storeHeader();
	}

	private void initOld(PropertySet ps) throws IOException
	{
		loadHeader();

		// only some of the properties may be changed.
		// the rest are just ignored.

		Object var;

		// tree variant.
		var = ps.getProperty("TreeVariant");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i != SpatialIndex.RtreeVariantLinear &&  i != SpatialIndex.RtreeVariantQuadratic && i != SpatialIndex.RtreeVariantRstar)
					throw new IllegalArgumentException("Property TreeVariant not a valid variant");
				m_treeVariant = i;
			}
			else
			{
				throw new IllegalArgumentException("Property TreeVariant must be an Integer");
			}
		}

		// near minimum overlap factor.
		var = ps.getProperty("NearMinimumOverlapFactor");
		if (var != null)
		{
			if (var instanceof Integer)
			{
				int i = ((Integer) var).intValue();
				if (i < 1 || i > m_indexCapacity || i > m_leafCapacity)
					throw new IllegalArgumentException("Property NearMinimumOverlapFactor must be less than both index and leaf capacities");
				m_nearMinimumOverlapFactor = i;
			}
			else
			{
				throw new IllegalArgumentException("Property NearMinimumOverlapFactor must be an Integer");
			}
		}

		// split distribution factor.
		var = ps.getProperty("SplitDistributionFactor");
		if (var != null)
		{
			if (var instanceof Double)
			{
				double f = ((Double) var).doubleValue();
				if (f <= 0.0f || f >= 1.0f)
					throw new IllegalArgumentException("Property SplitDistributionFactor must be in (0.0, 1.0)");
				m_splitDistributionFactor = f;
			}
			else
			{
				throw new IllegalArgumentException("Property SplitDistriburionFactor must be a Double");
			}
		}

		// reinsert factor.
		var = ps.getProperty("ReinsertFactor");
		if (var != null)
		{
			if (var instanceof Double)
			{
				double f = ((Double) var).doubleValue();
				if (f <= 0.0f || f >= 1.0f)
					throw new IllegalArgumentException("Property ReinsertFactor must be in (0.0, 1.0)");
				m_reinsertFactor = f;
			}
			else
			{
				throw new IllegalArgumentException("Property ReinsertFactor must be a Double");
			}
		}

		m_infiniteRegion.m_pLow = new double[m_dimension];
		m_infiniteRegion.m_pHigh = new double[m_dimension];

		for (int cDim = 0; cDim < m_dimension; cDim++)
		{
			m_infiniteRegion.m_pLow[cDim] = Double.POSITIVE_INFINITY;
			m_infiniteRegion.m_pHigh[cDim] = Double.NEGATIVE_INFINITY;
		}
	}

	private void storeHeader() throws IOException
	{
		ByteArrayOutputStream bs = new ByteArrayOutputStream();
		DataOutputStream ds = new DataOutputStream(bs);

		ds.writeInt(m_rootID);
		ds.writeInt(m_treeVariant);
		ds.writeDouble(m_fillFactor);
		ds.writeInt(m_indexCapacity);
		ds.writeInt(m_leafCapacity);
		ds.writeInt(m_nearMinimumOverlapFactor);
		ds.writeDouble(m_splitDistributionFactor);
		ds.writeDouble(m_reinsertFactor);
		ds.writeInt(m_dimension);
		ds.writeLong(m_stats.m_nodes);
		ds.writeLong(m_stats.m_data);
		ds.writeInt(m_stats.m_treeHeight);

		for (int cLevel = 0; cLevel < m_stats.m_treeHeight; cLevel++)
		{
			ds.writeInt(((Integer) m_stats.m_nodesInLevel.get(cLevel)).intValue());
		}

		ds.flush();
		m_headerID = m_pStorageManager.storeByteArray(m_headerID, bs.toByteArray());
	}

	private void loadHeader() throws IOException
	{
		byte[] data = m_pStorageManager.loadByteArray(m_headerID);
		DataInputStream ds = new DataInputStream(new ByteArrayInputStream(data));

		m_rootID = ds.readInt();
		m_treeVariant = ds.readInt();
		m_fillFactor = ds.readDouble();
		m_indexCapacity = ds.readInt();
		m_leafCapacity = ds.readInt();
		m_nearMinimumOverlapFactor = ds.readInt();
		m_splitDistributionFactor = ds.readDouble();
		m_reinsertFactor = ds.readDouble();
		m_dimension = ds.readInt();
		m_stats.m_nodes = ds.readLong();
		m_stats.m_data = ds.readLong();
		m_stats.m_treeHeight = ds.readInt();

		for (int cLevel = 0; cLevel < m_stats.m_treeHeight; cLevel++)
		{
			m_stats.m_nodesInLevel.add(new Integer(ds.readInt()));
		}
	}

	protected void insertData_impl(byte[] pData, Region mbr, int id)
	{
		//assert mbr.getDimension() == m_dimension;

		boolean[] overflowTable;

		Stack pathBuffer = new Stack();

		Node root = readNode(m_rootID);

		overflowTable = new boolean[root.m_level];
		for (int cLevel = 0; cLevel < root.m_level; cLevel++) overflowTable[cLevel] = false;

		Node l = root.chooseSubtree(mbr, 0, pathBuffer);
		l.insertData(pData, mbr, id, pathBuffer, overflowTable);

		m_stats.m_data++;
	}

	protected void insertData_impl(byte[] pData, Region mbr, int id, int level, boolean[] overflowTable)
	{
		//assert mbr.getDimension() == m_dimension;

		Stack pathBuffer = new Stack();

		Node root = readNode(m_rootID);
		Node n = root.chooseSubtree(mbr, level, pathBuffer);
		n.insertData(pData, mbr, id, pathBuffer, overflowTable);
	}

	protected boolean deleteData_impl(final Region mbr, int id)
	{
		//assert mbr.getDimension() == m_dimension;

		boolean bRet = false;

		Stack pathBuffer = new Stack();

		Node root = readNode(m_rootID);
		Leaf l = root.findLeaf(mbr, id, pathBuffer);

		if (l != null)
		{
			l.deleteData(id, pathBuffer);
			m_stats.m_data--;
			bRet = true;
		}

		return bRet;
	}

	protected int writeNode(Node n) throws IllegalStateException
	{
		byte[] buffer = null;

		try
		{
			buffer = n.store();
		}
		catch (IOException e)
		{
			System.err.println(e);
			throw new IllegalStateException("writeNode failed with IOException");
		}

		int page;
		if (n.m_identifier < 0) page = IStorageManager.NewPage;
		else page = n.m_identifier;

		try
		{
			page = m_pStorageManager.storeByteArray(page, buffer);
		}
		catch (InvalidPageException e)
		{
			System.err.println(e);
			throw new IllegalStateException("writeNode failed with InvalidPageException");
		}

		if (n.m_identifier < 0)
		{
			n.m_identifier = page;
			m_stats.m_nodes++;
			int i = ((Integer) m_stats.m_nodesInLevel.get(n.m_level)).intValue();
			m_stats.m_nodesInLevel.set(n.m_level, new Integer(i + 1));
		}

		m_stats.m_writes++;

		for (int cIndex = 0; cIndex < m_writeNodeCommands.size(); cIndex++)
		{
			((INodeCommand) m_writeNodeCommands.get(cIndex)).execute(n);
		}

		return page;
	}	
	
	protected Node readNode(int id)
	{
		byte[] buffer;
		DataInputStream ds = null;
		int nodeType = -1;
		Node n = null;

		try
		{
			buffer = m_pStorageManager.loadByteArray(id);
			ds = new DataInputStream(new ByteArrayInputStream(buffer));
			nodeType = ds.readInt();

			if (nodeType == SpatialIndex.PersistentIndex) n = new Index(this, -1, 0);
			else if (nodeType == SpatialIndex.PersistentLeaf) n = new Leaf(this, -1);
			else throw new IllegalStateException("readNode failed reading the correct node type information");

			n.m_pTree = this;
			n.m_identifier = id;
			n.load(buffer);

			m_stats.m_reads++;
		}
		catch (InvalidPageException e)
		{
			System.err.println(e);
			throw new IllegalStateException("readNode failed with InvalidPageException");
		}
		catch (IOException e)
		{
			System.err.println(e);
			throw new IllegalStateException("readNode failed with IOException");
		}

		for (int cIndex = 0; cIndex < m_readNodeCommands.size(); cIndex++)
		{
			((INodeCommand) m_readNodeCommands.get(cIndex)).execute(n);
		}

		return n;
	}

	protected void deleteNode(Node n)
	{
		try
		{
			m_pStorageManager.deleteByteArray(n.m_identifier);
		}
		catch (InvalidPageException e)
		{
			System.err.println(e);
			throw new IllegalStateException("deleteNode failed with InvalidPageException");
		}

		m_stats.m_nodes--;
		int i = ((Integer) m_stats.m_nodesInLevel.get(n.m_level)).intValue();
		m_stats.m_nodesInLevel.set(n.m_level, new Integer(i - 1));

		for (int cIndex = 0; cIndex < m_deleteNodeCommands.size(); cIndex++)
		{
			((INodeCommand) m_deleteNodeCommands.get(cIndex)).execute(n);
		}
	}

	private void rangeQuery(int type, final IShape query, final IVisitor v)
	{
		m_rwLock.read_lock();

		try
		{
			Stack st = new Stack();
			Node root = readNode(m_rootID);

			if (root.m_children > 0 && query.intersects(root.m_nodeMBR)) st.push(root);

			while (! st.empty())
			{
				Node n = (Node) st.pop();

				if (n.m_level == 0)
				{
					v.visitNode((INode) n);

					for (int cChild = 0; cChild < n.m_children; cChild++)
					{
						boolean b;
						if (type == SpatialIndex.ContainmentQuery) b = query.contains(n.m_pMBR[cChild]);
						else b = query.intersects(n.m_pMBR[cChild]);

						if (b)
						{
							Data data = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
							v.visitData(data);
							m_stats.m_queryResults++;
						}
					}
				}
				else
				{
					v.visitNode((INode) n);

					for (int cChild = 0; cChild < n.m_children; cChild++)
					{
						if (query.intersects(n.m_pMBR[cChild]))
						{
							st.push(readNode(n.m_pIdentifier[cChild]));
						}
					}
				}
			}
		}
		finally
		{
			m_rwLock.read_unlock();
		}
	}

	public String toString()
	{
		String s = "Dimension: " + m_dimension + "\n"
						 + "Fill factor: " + m_fillFactor + "\n"
						 + "Index capacity: " + m_indexCapacity + "\n"
						 + "Leaf capacity: " + m_leafCapacity + "\n";

		if (m_treeVariant == SpatialIndex.RtreeVariantRstar)
		{
			s += "Near minimum overlap factor: " + m_nearMinimumOverlapFactor + "\n"
				 + "Reinsert factor: " + m_reinsertFactor + "\n"
				 + "Split distribution factor: " + m_splitDistributionFactor + "\n";
		}

		s += "Utilization: " + 100 * m_stats.getNumberOfData() / (m_stats.getNumberOfNodesInLevel(0) * m_leafCapacity) + "%" + "\n"
			 + m_stats;

		return s;
	}

	class NNEntry
	{
		IEntry m_pEntry;
		double m_minDist;

		NNEntry(IEntry e, double f) { m_pEntry = e; m_minDist = f; }
	}

	class NNEntryComparator implements Comparator
	{
		public int compare(Object o1, Object o2)
		{
			NNEntry n1 = (NNEntry) o1;
			NNEntry n2 = (NNEntry) o2;

			if (n1.m_minDist < n2.m_minDist) return -1;
			if (n1.m_minDist > n2.m_minDist) return 1;
			return 0;
		}
	}

	class NNComparator implements INearestNeighborComparator
	{
		public double getMinimumDistance(IShape query, IEntry e)
		{
			IShape s = e.getShape();
			return query.getMinimumDistance(s);
		}
	}

	class ValidateEntry
	{
		Region m_parentMBR;
		Node m_pNode;

		ValidateEntry(Region r, Node pNode) { m_parentMBR = r; m_pNode = pNode; }
	}

	class Data implements IData
	{
		int m_id;
		Region m_shape;
		byte[] m_pData;

		Data(byte[] pData, Region mbr, int id) { m_id = id; m_shape = mbr; m_pData = pData; }

		public int getIdentifier() { return m_id; }
		public IShape getShape() { return new Region(m_shape); }
		public byte[] getData()
		{
			byte[] data = new byte[m_pData.length];
			System.arraycopy(m_pData, 0, data, 0, m_pData.length);
			return data;
		}
	}
	
	
public static void main(String[] args){
		ArrayList<ArrayList<Integer>> as = new ArrayList<>();
		for(int i = 0;i<3;++i){
			ArrayList<Integer> a = new ArrayList<Integer>();
			a.add(i);
			a.add(i+5);
		}
		permutation(as);
	}

public List<IRPoint> getSubsetInRadius(int radius, Point queryPoint,
		int queryCount)
{
	List<IRPoint> po = getSubsetInRadius(radius,queryPoint,queryCount,new NNComparator());
//	System.out.println(po.size());
	
	
	
	List<IRPoint> result = new ArrayList<>();

	

	//	ArrayList<ArrayList<IRPrayListoint>> result = new ArrayList<ArrayList<IRPoint>>();
////	if (query.getDimension() != m_dimension) throw new IllegalArgumentException("TOPK: IShape dimensionality is incorrect.");
//	NNComparator nnc = new NNComparator();
//	result = getSubsetInRadius(radius, queryPoint, queryCount, nnc);
	return result;
}

	public List<IRPoint> getSubsetInRadius(int radius, Point queryPoint,
		int querycount, final INearestNeighborComparator nnc)
		{
		List<IRPoint> po = new ArrayList<IRPoint>();

		try
		{
			ArrayList queue = new ArrayList();
			ArrayList Result_queue = new ArrayList();

			Node n = readNode(m_rootID);
			queue.add(new NNEntry(n, 0.0));//加入根节点，距离为0

			int count = 0;
			double knearest = 0.0;

			while (queue.size() != 0)
			{
				NNEntry first = (NNEntry) queue.remove(0);

				if (first.m_pEntry instanceof Node)//如果first.m_pEntry是Node的一个实例
				{
					n = (Node) first.m_pEntry;
//					v.visitNode((INode) n);//访问节点n

					for (int cChild = 0; cChild < n.m_children; cChild++)
					{
						IEntry e;

						if (n.m_level == 0)
						{
							e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);
							//如果此节点是叶子节点，则依次读出叶子节点当中的对象的值
							String word = new String(n.m_pData[cChild]);//如果此节点是叶节点，则找出此实例当中的关键字						
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(queryPoint, e));
							//如果e在这里代表着对象，则表示查询点与此对象之间的最小距离；如果e在这里代表着节点，则表示查询点与此节点的MBR之间的距离
							
							// Why don't I use a TreeSet here? See comment above...
							Result_queue.add(nnc.getMinimumDistance(queryPoint, e));//将关键字相同的点之间的距离存到数组当中
							
//							int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
//							if (loc >= 0) queue.add(loc, e2);//增加条件，使得查询到的符合条件的点的关键字与查询的关键字相同
//							else
							queue.add(e2);									
						}
						else
						{
							e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点
							
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(queryPoint, e));
							
//							int loc = Collections.binarySearch(queue, e2, new NNEntryComparator());
//							if (loc >= 0) queue.add(loc, e2);
//							else 
							queue.add(e2);
						}
					}
				}
				else
				{
					// report all nearest neighbors with equal furthest distances.
					// (neighbors can be more than k, if many happen to have the same
					//  furthest distance).
//					if (count >= k && first.m_minDist > knearest) break;
					if(first.m_minDist > radius)
						continue;
//					v.visitData((IData) first.m_pEntry);

					IRPoint temp = new IRPoint();
					temp.m_pointID = first.m_pEntry.getIdentifier();
					//temp.Keyword = keyword;
					temp.distance = first.m_minDist;

					po.add(temp);
					m_stats.m_queryResults++;
					count++;
					knearest = first.m_minDist;
//					System.out.println(first.m_minDist);
				}
			}
			
//			Collections.sort(Result_queue);
//			System.out.println("Key Distances:" + po);
		}
		finally
		{
			m_rwLock.read_unlock();
		}
		return po;
		}

		//递归求所有的子集
		void subsets(int queryCount, ArrayList<IRPoint> allObject, ArrayList<IRPoint> temp,int level, ArrayList<ArrayList<IRPoint>> result)
		{
//		  //如果是叶子节点则加入到result中
//		  if(level == queryCount)
//		  {
//		    result.push_back(temp);
//		    return;
//		  }
//
//		  subsets(queryCount,allObject,temp,level + 1,result);
//
//		  temp.push_back(S[level]);
//		  subsets(queryCount,allObject,temp,level + 1,result);
		}
		
//		public ArrayList<IRPoint> Replace(Point queryPoint, List<IRPoint> po)
//		{
//			ArrayList<IRPoint> result = new ArrayList<IRPoint>();
////			if (query.getDimension() != m_dimension) throw new IllegalArgumentException("TOPK: IShape dimensionality is incorrect.");
//			NNComparator nnc = new NNComparator();
//			result = Replace(queryPoint, po, nnc);
//			return result;
//		}
	
		public void Replace(Point queryPoint, List<IRPoint> po){
			List<IRPoint> result = new ArrayList<IRPoint>();
			IRPoint temp;
			double marginalGain=0;
			
			result = readallIRPoint(queryPoint);
			for(IRPoint irp: po) {
				
			}

		}
		
		public List<IRPoint> readallIRPoint(Point queryPoint)
		{
			List<IRPoint> result = new ArrayList<IRPoint>();
//			if (query.getDimension() != m_dimension) throw new IllegalArgumentException("TOPK: IShape dimensionality is incorrect.");
			NNComparator nnc = new NNComparator();
			result = readallIRPoint(queryPoint, nnc);
			return result;
		}

		public List<IRPoint> readallIRPoint(Point queryPoint, final INearestNeighborComparator nnc) {
			List<IRPoint> po = new ArrayList<IRPoint>();

			try
			{
				ArrayList queue = new ArrayList();
//				ArrayList Result_queue = new ArrayList();

				Node n = readNode(m_rootID);
				queue.add(new NNEntry(n, 0.0));//加入根节点，距离为0

//				int count = 0;

				while (queue.size() != 0)
				{
					NNEntry first = (NNEntry) queue.remove(0);

					if (first.m_pEntry instanceof Node)//如果first.m_pEntry是Node的一个实例
					{
						n = (Node) first.m_pEntry;

						for (int cChild = 0; cChild < n.m_children; cChild++)
						{
							IEntry e;

							if (n.m_level == 0)
							{
								e = new Data(n.m_pData[cChild], n.m_pMBR[cChild], n.m_pIdentifier[cChild]);								
							}
							else
							{
								e = (IEntry) readNode(n.m_pIdentifier[cChild]);//如果此节点不是叶子节点，则依次读叶子节点
							}
							NNEntry e2 = new NNEntry(e, nnc.getMinimumDistance(queryPoint, e));
							queue.add(e);
						}
					}
					else
					{
						IRPoint temp = new IRPoint();
						temp.m_pointID = first.m_pEntry.getIdentifier();
						temp.distance = first.m_minDist;

						po.add(temp);
//						m_stats.m_queryResults++;
//						count++;
					}
				}
				
//				Collections.sort(Result_queue);
				System.out.println("point number:" + po.size());
			}
			finally
			{
				m_rwLock.read_unlock();
			}
			return po;
		}
}


