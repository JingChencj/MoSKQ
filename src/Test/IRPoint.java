package Test;

import java.util.ArrayList;

public class IRPoint implements Cloneable{

	int dimention = 20;
	
	public int m_pointID;
	public double[] m_pCoordinate = new double[2];  //ÿ���������
	public String Keyword;
	public double distance;
	
	public float[] topics = new float[dimention];
	
	public ArrayList<String> bucketnumber = new ArrayList<String>();//�����洢����ռ��ı������Ӧ��Ͱ��

	public double getLatitude(){
		return m_pCoordinate[0];
	}

	public double getLongitude(){
		return m_pCoordinate[1];
	}
}
