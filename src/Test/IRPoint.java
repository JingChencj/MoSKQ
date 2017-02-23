package Test;

import java.util.ArrayList;

public class IRPoint implements Cloneable{

	int dimention = 20;
	
	public int m_pointID;
	public double[] m_pCoordinate = new double[2];  //每个点的坐标
	public String Keyword;
	public double distance;
	
	public float[] topics = new float[dimention];
	
	public ArrayList<String> bucketnumber = new ArrayList<String>();//用来存储这个空间文本对象对应的桶号

	public double getLatitude(){
		return m_pCoordinate[0];
	}

	public double getLongitude(){
		return m_pCoordinate[1];
	}
}
