package spatialindex.Ngram;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * K-gram 
 */

public class Ngram {
	public static final int ED_THRESHOLD = 2;
	public static final int DEFAULT_N = 3;
	/**
	 * 鍋滅敤璇�
	 */
	public static Set<String> stopWords = new HashSet<String>(10);
	static {
		stopWords.add("&");
		stopWords.add(" ");
		stopWords.add("\\.");
		stopWords.add("'");
//		// 绉婚櫎鏍囩偣绗﹀彿
		stopWords.add("銆恷銆憒\\(|\\)|锛坾锛墊\\+|-|\\*|/");
	}
	
	/**
	 * 鍒犻櫎鍋滅敤璇�
	 * 
	 * @param word
	 * @return
	 */
	public static String removeStopWords(String word) {
		if (word == null) {
			return null;
		}
		for (String stopWord : stopWords) {
			word = word.replaceAll(stopWord, "");
		}
		return word;
	}

	/**
	 * 鑾峰緱瀛楃涓茬殑K-gram瀛楃涓�
	 * 
	 * @param str
	 * @param k
	 * @return
	 */
	public static List<String> getKGramsList(String str, int k)
	{
		if (str == null || str.trim().isEmpty() || k < 1)//trim鏂规硶杩斿洖鍘婚櫎棣栧熬绌虹櫧瀛楃鐨勫瓧绗︿覆
		{
			return null;
		}
		str = removeStopWords(str);//鍘婚櫎鍋滅敤璇�
		int length = str.length();//鍙栧緱瀛楃涓茬殑闀垮害
		if(k > length)
			return null;
		List<String> grams = new ArrayList<String>(length);
		for (int i = 0; i < length - k + 1; i++) 
		{
            //鎴彇鍘熻瘝()
            grams.add(str.substring(i, i + k));
    }
		return grams;
	}
	
	/**
	 * 鑾峰緱涓嶉噸澶嶇殑K-gram瀛楃涓�
	 * 
	 * @param str
	 * @param k
	 * @return
	 */
	public static Set<String> getKGramsSet(String str, int k) 
	{
		if (str == null || str.isEmpty()) {
			return null;
		}
		return new HashSet<String>(getKGramsList(str, k));
	}
	
	/**
	 * 鑾峰緱鑻ュ共瀛楃涓茬殑鍏叡K-Gram瀛楃涓�
	 * 
	 * @param k
	 * @param str
	 * @return
	 */
	public static Set<String> getNKGramsSet(int k, String...str) 
	{
		if (str == null) {
			return null;
		}
		Set<String> grams = getKGramsSet(str[0], k);
		for (int index = 1 ; index < str.length ; index ++) {
			grams.addAll(getNKGramsSet(k, str[index]));
		}
		return grams;
	}
	
	/**
	 * 鑾峰緱鐩镐氦鐨凨-grams鐨勪釜鏁�
	 * 
	 * @param set1
	 * @param set2
	 * @return
	 */
	public static int getKGramsSize(Set<String> set1, Set<String> set2) {
		if (set1 == null || set2 == null) {
			return 0;
		}
		// 涓嶅彲鏀瑰彉鍘熼泦鍚�
		Set<String> set = new HashSet<String>(set1);
		// 鏌ユ壘鐩镐氦鐨凨-grams涓暟
		set.retainAll(set2);
		return set.size();
	}
	
	public static void main(String[] args) {
		System.out.println(Ngram.getKGramsList("knowledge", 2));
	}
}
