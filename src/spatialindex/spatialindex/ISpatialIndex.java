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
// License along with this library; if not, write to the Free Software
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

package spatialindex.spatialindex;


import Test.IRPoint;
import spatialindex.storagemanager.PropertySet;

import java.util.ArrayList;
import java.util.List;

public interface ISpatialIndex
{
	public void flush() throws IllegalStateException;
	public void insertData(final byte[] data, final IShape shape, int id);
	public boolean deleteData(final IShape shape, int id);
	public void containmentQuery(final IShape query, final IVisitor v);
	public void intersectionQuery(final IShape query, final IVisitor v);
	public void pointLocationQuery(final IShape query, final IVisitor v);
	public List<IRPoint> nearestNeighborQuery(int k, final IShape query, String keyword, final IVisitor v, INearestNeighborComparator nnc);
	public List<IRPoint> nearestNeighborQuery(int k, final IShape query, String keyword, final IVisitor v);
	public void queryStrategy(final IQueryStrategy qs);
	public PropertySet getIndexProperties();
	public void addWriteNodeCommand(INodeCommand nc);
	public void addReadNodeCommand(INodeCommand nc);
	public void addDeleteNodeCommand(INodeCommand nc);
	public boolean isIndexValid();
	public IStatistics getStatistics();
	
	
	//2016.2.24
	public void BuildInvertedList();
	public double getCost(List<IRPoint> po);

	//2016.6.23
	public List<IRPoint> CollectiveQuery(final List<IShape> query, List<String> keyword, final IVisitor v);
	public List<IRPoint> CollectiveQuery(final List<IShape> query, List<String> keyword, final IVisitor v, INearestNeighborComparator nnc);
	
	//2016.6.28
	public List<IRPoint> CollectiveQuery_Type2(final List<IShape> query, List<String> keyword, final IVisitor v);
	public List<IRPoint> CollectiveQuery_Type2(final List<IShape> query, List<String> keyword, final IVisitor v, INearestNeighborComparator nnc);
	
	//2016.7.13
	public List<IRPoint> CollectiveQuery_Type3(final List<IShape> query, List<String> keyword, final IVisitor v);
	public List<IRPoint> CollectiveQuery_Type3(final List<IShape> query, List<String> keyword, final IVisitor v, INearestNeighborComparator nnc);
	
	public List<List<IRPoint>> CollectiveQuery_temp(final IShape query, List<String> keyword, final IVisitor v);
	public List<List<IRPoint>> CollectiveQuery_temp(final IShape query, List<String> keyword, final IVisitor v, INearestNeighborComparator nnc);
	
	//2016.8.15
	public List<IRPoint> MoSKQ(final ArrayList<String> query, final IVisitor v);
	public List<IRPoint> MoSKQ(final ArrayList<String> query, final IVisitor v, final INearestNeighborComparator nnc);

	//2016.10.24
	public List<IRPoint> TOPK(int k, final IShape query, final IVisitor v, INearestNeighborComparator nnc);
	public List<IRPoint> TOPK(int k, final IShape query, final IVisitor v);
	
	public List<IRPoint> getSubsetInRadius(int radius, Point queryPoint, int queryCount);
	public List<IRPoint> getSubsetInRadius(int radius, Point queryPoint, int querycount, final INearestNeighborComparator nnc);
	public void Replace(Point queryPoint, List<IRPoint> po);
	public List<IRPoint> readallIRPoint(Point queryPoint);
	public List<IRPoint> readallIRPoint(Point queryPoint, final INearestNeighborComparator nnc);
} // ISpatialIndex

