import java.util.*;
import java.util.Map.Entry;
import java.io.*;

public class TopNFeatureSelector {
	public TopNFeatureSelector(){
		
	}
	
	/**
	 * This sets up all variables for the topNFeatureSelection, this must run first
	 * It basically goes through every document in the directory structure, counts the number
	 * of documents, counts the occurances of each word in every document, counts the number
	 * of documents that have at least one occurance of each word, and counts the total
	 * number of words for each document. Words are also considered to be non-header 
	 * alphabetic characters 
	 * @param folderPath
	 * @throws Exception
	 */
	public void loadDocumentsFromFolder(String folderPath, List<Set<Integer>> documentLabels) throws Exception{
		File documentDirectory = new File(folderPath);
		if(!documentDirectory.isDirectory()){
			throw new Exception("Passed Folder Path is not a directory");
		}
		
		int folderNumber = 0;
		int documentNumber = 0;
		//  navigate past group folder to the folders of documets grouped by labels
		for(File labeledFolder : documentDirectory.listFiles()[0].listFiles()){
			if(!labeledFolder.isDirectory()){
				throw new Exception("Incorrect Directory Structure");
			}
			documentLabels.add(new HashSet<Integer>());
			
			System.out.println("Parsing new Labeled Folder");
			

			for(File file : labeledFolder.listFiles()){
				documentParser(file);
				documentLabels.get(folderNumber).add(documentNumber);
				documentNumber++;
			}
			folderNumber++;
		}
		
//		//  after this has execute the documentsWithWord map should have keys for every possible word
//		//  all documentsCounts must be checked then to see how many documents contain each possible word
//		//  all these variables will be used to computer the top N features when that function is called
//		Map<String, Integer> documentsWithWordTemp = new HashMap<String, Integer>();
//		for(Iterator<Entry<String, Integer>> itr = documentsWithWord.entrySet().iterator(); itr.hasNext();){
//			
//			// save the word and add it to the temp list of counts
//			//  a temp list of counts is needed because the original map can't be edited while we're iterating through
//			String word = itr.next().getKey();
//			documentsWithWordTemp.put(word, 0);
//			//  check if the word is in each document
//			for(int i = 0; i < documentCounts.size(); i++){
//				//  if the word is in the document, increment the count of documents with the word
//				if(documentCounts.get(i).containsKey(word)){
//					documentsWithWordTemp.put(word, documentsWithWordTemp.get(word) + 1);
//				}
//			}
//		}
//		documentsWithWord = documentsWithWordTemp;
		
		
		for(int i = 0; i < documentCounts.size(); i++){
			for(Iterator<Entry<String, Integer>> itr = documentCounts.get(i).entrySet().iterator(); itr.hasNext();){
				String word = itr.next().getKey();
				
				documentsWithWord.put(word, documentsWithWord.get(word) + 1);
			}
		}
		
		System.out.println("Total number of documents: " + documentCounts.size());
		System.out.println("Total number of words: " + documentsWithWord.size());
	}
	
	
	
	/**
	 * Takes a list of documents of words and returns a list of documents of word counts
	 * for each document
	 * @param documents
	 * @return
	 * @throws FileNotFoundException 
	 */
	private void documentParser(File file) throws FileNotFoundException{
		m++;   //  total number of documents
		int wordCountForDocument = 0;
		Map<String, Integer> countsForThisDocument = new HashMap<String, Integer>();
		
		Scanner scan = null;
		try{
			
		//  makes a set of all possible words in this document
			scan = new Scanner(file);
			String line = "";
			boolean passedHeaders = false;
			while(scan.hasNextLine()){
				line = scan.nextLine();
				if(passedHeaders){
					//  after headers are passed add every alphabetic word to the dictionary
					for(String token : line.split("\\s+")){
						//  if the token is an alphabetic string then add it
						if(isAlpha(token)){
							wordCountForDocument++;
							
							//  check to see if the word has already been added to global vocabulary
							//  add the word to the map if it has not already been added
							if(!documentsWithWord.containsKey(token)){
								documentsWithWord.put(token, 0);
							}
							//  if this word has not already been see in the document then add it to the map with count 1
							//  otherwise increment it's value to reflect another instance of the word found in the document
							if(!countsForThisDocument.containsKey(token)){
								countsForThisDocument.put(token, 1);
							}else{
								//  increment count
								countsForThisDocument.put(token, countsForThisDocument.get(token) + 1);
							}
						}
					}

				}else{
					//  need to pass an empty line to get past the headers
					if(line.equals("")){
						passedHeaders = true;
					}
				}
			}
			totalWordCountsForDocument.add(wordCountForDocument);  //  add the total number of words in this document to its
																   //  corresponding index in totalWordCountsForDocument
			documentCounts.add(countsForThisDocument);  //  add the map of word counts to this document's index in document counts
			
			
		}finally{
			if(scan != null)
				scan.close();
		}
		
	}
	
	/**
	 * Returns true if the passed string is strictly alphabetic
	 * @param s string to be tested
	 * @return
	 */
	private boolean isAlpha(String s){
		return s.matches("[a-zA-z]+");
	}
	
	
	
	
	public void topNFeatureSelection(int passedN){
		this.n = passedN;
		
		System.out.println("Calculating word weights");
		caluculateDocumentWordWeights();
		System.out.println("Word weights finished");
		
		for(int i = 0; i < documentWordWeights.size(); i++){
			Map<String, Double> unsortedMap = documentWordWeights.get(i);
			
			// Convert Map to List
			//  TODO: this might cause problems
			List<Map.Entry<String, Double>> sortedList = 
				new LinkedList<Map.Entry<String, Double>>(unsortedMap.entrySet());
	 
			// Sort list with comparator, to compare the Map values
			Collections.sort(sortedList, new Comparator<Map.Entry<String, Double>>() {
				public int compare(Map.Entry<String, Double> o1,
	                                           Map.Entry<String, Double> o2) {
//						return (o1.getValue()).compareTo(o2.getValue());
					if(o1.getValue() < o2.getValue()){
						return 1;
					}else if(o1.getValue() > o2.getValue()){
						return -1;
					}else{
						return 0;
					}
				}
			});
			
			//  add top n words to p
			for(int k = 0; k < sortedList.size() && k < n; k++){
				p.add(sortedList.get(k).getKey());
			}
		}
		
		
		System.out.println("Size of P: " + p.size());
//		System.out.println(p);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		totalWordCountsForDocument = null;
		documentsWithWord = null;
		documentWordWeights = null;
		System.gc();
		
		
		
		
		
		//  build a list of top document counts the size of the number of documents
		for(int i = 0; i < documentCounts.size(); i++){
			topNDocumentCounts.add(new ArrayList<Integer>());
		}
		
		
		//  then build topNDocumentCounts and topNWordOrder
		//  for each word in the global p dictionary, go through each document and set their count
		//  for that specific word
		for(Iterator<String> itr = p.iterator(); itr.hasNext();){
			String word = itr.next();  //  next word to get counts for
			//  add to the word Order List
			topNWordOrder.add(word);
			
			//  if the document has occuraces of the word, set its counts accordingly, 0 otherwise
			for(int i = 0; i < documentCounts.size(); i++){

				if(documentCounts.get(i).containsKey(word)){
					topNDocumentCounts.get(i).add(documentCounts.get(i).get(word));
				}else{
					topNDocumentCounts.get(i).add(0);
				}
			}
		}
		
		
		
		
		
		
		
		System.out.println("Finished topN");
	}
	
	
	
	/**
	 * Calculates weights for all word types in each document
	 */
	private void caluculateDocumentWordWeights(){
		//  first calculate all weights/ranks
		//  go through every document and calculate the ranks for every word
		for(int i = 0; i < documentCounts.size(); i++){
			
			Map<String, Double> documentWeights = new HashMap<String, Double>();
			Entry<String, Integer> entry = null;
			double weight = 0.0;
			for(Iterator<Entry<String, Integer>> itr = documentCounts.get(i).entrySet().iterator(); itr.hasNext();){
				entry = itr.next();
				
				//  if the document doesn't have any of this type of word, then set it's weight to 0
				if(entry.getValue() < 1){
					documentWeights.put(entry.getKey(), 0.0);
				}else{
					//  = 1 + frequency of word (count) over total count of words in the documnet 
					//  * m (the number of documents) / the number of documents where word (entry.getKey()) appeaers 
//					weight = ((1.0 + Math.log((double)entry.getValue() / (double)totalWordCountsForDocument.get(i))) * Math.log(m / documentsWithWord.get(entry.getKey())));
					weight = ((1.0 + Math.log((double)entry.getValue())) * Math.log(m / documentsWithWord.get(entry.getKey())));

					documentWeights.put(entry.getKey(), weight);
				}
			}
			documentWordWeights.add(documentWeights);
		}
	}
	
	
	
	
	
	
	
	
	
	
	/**
	 * This returns the top N documents counts, in order for this to be accurate
	 * loadDocumentsFromFolder and topNFeatureSelection must have been called previous
	 * @return a list of counts for each topN word type for each document
	 */
	public List<List<Integer>> getTopNDocumentCounts(){
		return this.topNDocumentCounts;
	}
	
	/**
	 * Converts the list's ints to doubles and returns the same list
	 * but in double format
	 * @return
	 */
	public List<List<Double>> getTopNDocumentCountsDouble(){
		List<List<Double>> doubleListToReturn = new ArrayList<List<Double>>();
		
		for(int i = 0; i < topNDocumentCounts.size(); i++){
			List<Double> documentDoubleList = new ArrayList<Double>();
			for(int j = 0; j < topNDocumentCounts.get(i).size(); j++){
				documentDoubleList.add((double)topNDocumentCounts.get(i).get(j));
			}
			doubleListToReturn.add(documentDoubleList);
		}
		return doubleListToReturn;
	}
	
	/**
	 * This returns the top N word order the top N document counts follows for each document
	 * loadDocumentsFromFolder and topNFeatureSelection must have been called previous to this function
	 * @return the word strings associated with each index of topNDocumentCounts, this can be used
	 * for output or for coordinating preset betas etc since the order returned will be arbitrarily
	 * determined by the internal set used
	 */
	public List<String> getTopNWordOrder(){
		return this.topNWordOrder;
	}
	
	
	
	
	
	
	
	
	
	private List<Map<String, Integer>> documentCounts = new ArrayList<Map<String,Integer>>();
	private List<Integer> totalWordCountsForDocument = new ArrayList<Integer>();
	private int m = 0;  //  document count
	private Map<String, Integer> documentsWithWord = new HashMap<String, Integer>();
	
	
	private int n = 0;  //  top N features from each document
	private List<Map<String, Double>> documentWordWeights = new ArrayList<Map<String, Double>>();
	private Set<String> p = new HashSet<String>();  //  feature pool
	private List<List<Integer>> topNDocumentCounts = new ArrayList<List<Integer>>();
	private List<String> topNWordOrder = new ArrayList<String>();
}




























//
//import java.util.*;
//import java.util.Map.Entry;
//import java.io.*;
//
//public class TopNFeatureSelector {
//	public TopNFeatureSelector(){
//		
//	}
//	
//	/**
//	 * This sets up all variables for the topNFeatureSelection, this must run first
//	 * It basically goes through every document in the directory structure, counts the number
//	 * of documents, counts the occurances of each word in every document, counts the number
//	 * of documents that have at least one occurance of each word, and counts the total
//	 * number of words for each document. Words are also considered to be non-header 
//	 * alphabetic characters 
//	 * @param folderPath
//	 * @throws Exception
//	 */
//	public void loadDocumentsFromFolder(String folderPath, List<Integer> documentLabels) throws Exception{
//		File documentDirectory = new File(folderPath);
//		if(!documentDirectory.isDirectory()){
//			throw new Exception("Passed Folder Path is not a directory");
//		}
//		
//		int folderNumber = 0;
//		//  navigate past group folder to the folders of documets grouped by labels
//		for(File labeledFolder : documentDirectory.listFiles()[0].listFiles()){
//			if(!labeledFolder.isDirectory()){
//				throw new Exception("Incorrect Directory Structure");
//			}
//			
//			System.out.println("Parsing new Labeled Folder");
//			
//
//			for(File file : labeledFolder.listFiles()){
//				documentParser(file);
//				documentLabels.add(folderNumber);
//			}
//			folderNumber++;
//		}
//		
////		//  after this has execute the documentsWithWord map should have keys for every possible word
////		//  all documentsCounts must be checked then to see how many documents contain each possible word
////		//  all these variables will be used to computer the top N features when that function is called
////		Map<String, Integer> documentsWithWordTemp = new HashMap<String, Integer>();
////		for(Iterator<Entry<String, Integer>> itr = documentsWithWord.entrySet().iterator(); itr.hasNext();){
////			
////			// save the word and add it to the temp list of counts
////			//  a temp list of counts is needed because the original map can't be edited while we're iterating through
////			String word = itr.next().getKey();
////			documentsWithWordTemp.put(word, 0);
////			//  check if the word is in each document
////			for(int i = 0; i < documentCounts.size(); i++){
////				//  if the word is in the document, increment the count of documents with the word
////				if(documentCounts.get(i).containsKey(word)){
////					documentsWithWordTemp.put(word, documentsWithWordTemp.get(word) + 1);
////				}
////			}
////		}
////		documentsWithWord = documentsWithWordTemp;
//		
//		
//		for(int i = 0; i < documentCounts.size(); i++){
//			for(Iterator<Entry<String, Integer>> itr = documentCounts.get(i).entrySet().iterator(); itr.hasNext();){
//				String word = itr.next().getKey();
//				
//				documentsWithWord.put(word, documentsWithWord.get(word) + 1);
//			}
//		}
//		
//		System.out.println("Total number of documents: " + documentCounts.size());
//		System.out.println("Total number of words: " + documentsWithWord.size());
//	}
//	
//	
//	
//	/**
//	 * Takes a list of documents of words and returns a list of documents of word counts
//	 * for each document
//	 * @param documents
//	 * @return
//	 * @throws FileNotFoundException 
//	 */
//	private void documentParser(File file) throws FileNotFoundException{
//		m++;   //  total number of documents
//		int wordCountForDocument = 0;
//		Map<String, Integer> countsForThisDocument = new HashMap<String, Integer>();
//		
//		Scanner scan = null;
//		try{
//			
//		//  makes a set of all possible words in this document
//			scan = new Scanner(file);
//			String line = "";
//			boolean passedHeaders = false;
//			while(scan.hasNextLine()){
//				line = scan.nextLine();
//				if(passedHeaders){
//					//  after headers are passed add every alphabetic word to the dictionary
//					for(String token : line.split("\\s+")){
//						//  if the token is an alphabetic string then add it
//						if(isAlpha(token)){
//							wordCountForDocument++;
//							
//							//  check to see if the word has already been added to global vocabulary
//							//  add the word to the map if it has not already been added
//							if(!documentsWithWord.containsKey(token)){
//								documentsWithWord.put(token, 0);
//							}
//							//  if this word has not already been see in the document then add it to the map with count 1
//							//  otherwise increment it's value to reflect another instance of the word found in the document
//							if(!countsForThisDocument.containsKey(token)){
//								countsForThisDocument.put(token, 1);
//							}else{
//								//  increment count
//								countsForThisDocument.put(token, countsForThisDocument.get(token) + 1);
//							}
//						}
//					}
//
//				}else{
//					//  need to pass an empty line to get past the headers
//					if(line.equals("")){
//						passedHeaders = true;
//					}
//				}
//			}
//			totalWordCountsForDocument.add(wordCountForDocument);  //  add the total number of words in this document to its
//																   //  corresponding index in totalWordCountsForDocument
//			documentCounts.add(countsForThisDocument);  //  add the map of word counts to this document's index in document counts
//			
//			
//		}finally{
//			if(scan != null)
//				scan.close();
//		}
//		
//	}
//	
//	/**
//	 * Returns true if the passed string is strictly alphabetic
//	 * @param s string to be tested
//	 * @return
//	 */
//	private boolean isAlpha(String s){
//		return s.matches("[a-zA-z]+");
//	}
//	
//	
//	
//	
//	public void topNFeatureSelection(int passedN){
//		this.n = passedN;
//		
//		System.out.println("Calculating word weights");
//		caluculateDocumentWordWeights();
//		System.out.println("Word weights finished");
//		
//		for(int i = 0; i < documentWordWeights.size(); i++){
//			Map<String, Double> unsortedMap = documentWordWeights.get(i);
//			
//			// Convert Map to List
//			//  TODO: this might cause problems
//			List<Map.Entry<String, Double>> sortedList = 
//				new LinkedList<Map.Entry<String, Double>>(unsortedMap.entrySet());
//	 
//			// Sort list with comparator, to compare the Map values
//			Collections.sort(sortedList, new Comparator<Map.Entry<String, Double>>() {
//				public int compare(Map.Entry<String, Double> o1,
//	                                           Map.Entry<String, Double> o2) {
////						return (o1.getValue()).compareTo(o2.getValue());
//					if(o1.getValue() < o2.getValue()){
//						return 1;
//					}else if(o1.getValue() > o2.getValue()){
//						return -1;
//					}else{
//						return 0;
//					}
//				}
//			});
//			
//			//  add top n words to p
//			for(int k = 0; k < sortedList.size() && k < n; k++){
//				p.add(sortedList.get(k).getKey());
//			}
//		}
//		
//		
//		System.out.println("Size of P: " + p.size());
////		System.out.println(p);
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		
//		totalWordCountsForDocument = null;
//		documentsWithWord = null;
//		documentWordWeights = null;
//		System.gc();
//		
//		
//		topNDocumentCounts = new int[documentCounts.size()][p.size()];
//		topNWordOrder= new String[p.size()];
//		int pIndex = 0;
//		for(Iterator<String> itr = p.iterator(); itr.hasNext();){
//			String word = itr.next();  //  next word to get counts for
//			//  add to the word Order List
////			topNWordOrder.add(word);
//			topNWordOrder[pIndex] = word;
//			
//			//  if the document has occuraces of the word, set its counts accordingly, 0 otherwise
//			for(int i = 0; i < documentCounts.size(); i++){
//	
//				if(documentCounts.get(i).containsKey(word)){
////					topNDocumentCounts.get(i).add(documentCounts.get(i).get(word));
//					topNDocumentCounts[i][pIndex] = documentCounts.get(i).get(word);
//				}else{
////					topNDocumentCounts.get(i).add(0);
//					topNDocumentCounts[i][pIndex] = 0;
//				}
//			}
//			pIndex++;
//		}
//		
//		
//		
//		
//		
//		
////		//  build a list of top document counts the size of the number of documents
////		for(int i = 0; i < documentCounts.size(); i++){
////			topNDocumentCounts.add(new ArrayList<Integer>());
////		}
////		
////		
////		//  then build topNDocumentCounts and topNWordOrder
////		//  for each word in the global p dictionary, go through each document and set their count
////		//  for that specific word
////		for(Iterator<String> itr = p.iterator(); itr.hasNext();){
////			String word = itr.next();  //  next word to get counts for
////			//  add to the word Order List
////			topNWordOrder.add(word);
////			
////			//  if the document has occuraces of the word, set its counts accordingly, 0 otherwise
////			for(int i = 0; i < documentCounts.size(); i++){
////
////				if(documentCounts.get(i).containsKey(word)){
////					topNDocumentCounts.get(i).add(documentCounts.get(i).get(word));
////				}else{
////					topNDocumentCounts.get(i).add(0);
////				}
////			}
////		}
//		
//		
//		
//		
//		
//		
//		
//		System.out.println("Finished topN");
//		
//	}
//	
//	
//	
//	/**
//	 * Calculates weights for all word types in each document
//	 */
//	private void caluculateDocumentWordWeights(){
//		//  first calculate all weights/ranks
//		//  go through every document and calculate the ranks for every word
//		for(int i = 0; i < documentCounts.size(); i++){
//			
//			Map<String, Double> documentWeights = new HashMap<String, Double>();
//			Entry<String, Integer> entry = null;
//			double weight = 0.0;
//			for(Iterator<Entry<String, Integer>> itr = documentCounts.get(i).entrySet().iterator(); itr.hasNext();){
//				entry = itr.next();
//				
//				//  if the document doesn't have any of this type of word, then set it's weight to 0
//				if(entry.getValue() < 1){
//					documentWeights.put(entry.getKey(), 0.0);
//				}else{
//					//  = 1 + frequency of word (count) over total count of words in the documnet 
//					//  * m (the number of documents) / the number of documents where word (entry.getKey()) appeaers 
////					weight = ((1.0 + Math.log((double)entry.getValue() / (double)totalWordCountsForDocument.get(i))) * Math.log(m / documentsWithWord.get(entry.getKey())));
//					weight = ((1.0 + Math.log((double)entry.getValue())) * Math.log(m / documentsWithWord.get(entry.getKey())));
//
//					documentWeights.put(entry.getKey(), weight);
//				}
//			}
//			documentWordWeights.add(documentWeights);
//		}
//	}
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	/**
//	 * This returns the top N documents counts, in order for this to be accurate
//	 * loadDocumentsFromFolder and topNFeatureSelection must have been called previous
//	 * @return a list of counts for each topN word type for each document
//	 */
//	public int[][] getTopNDocumentCounts(){
//		return this.topNDocumentCounts;
//	}
//	
////	/**
////	 * Converts the list's ints to doubles and returns the same list
////	 * but in double format
////	 * @return
////	 */
////	public List<List<Double>> getTopNDocumentCountsDouble(){
////		List<List<Double>> doubleListToReturn = new ArrayList<List<Double>>();
////		
////		for(int i = 0; i < topNDocumentCounts.size(); i++){
////			List<Double> documentDoubleList = new ArrayList<Double>();
////			for(int j = 0; j < topNDocumentCounts.get(i).size(); j++){
////				documentDoubleList.add((double)topNDocumentCounts.get(i).get(j));
////			}
////			doubleListToReturn.add(documentDoubleList);
////		}
////		return doubleListToReturn;
////	}
//	
//	/**
//	 * This returns the top N word order the top N document counts follows for each document
//	 * loadDocumentsFromFolder and topNFeatureSelection must have been called previous to this function
//	 * @return the word strings associated with each index of topNDocumentCounts, this can be used
//	 * for output or for coordinating preset betas etc since the order returned will be arbitrarily
//	 * determined by the internal set used
//	 */
//	public String[]  getTopNWordOrder(){
//		return this.topNWordOrder;
//	}
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	private List<Map<String, Integer>> documentCounts = new ArrayList<Map<String,Integer>>();
//	private List<Integer> totalWordCountsForDocument = new ArrayList<Integer>();
//	private int m = 0;  //  document count
//	private Map<String, Integer> documentsWithWord = new HashMap<String, Integer>();
//	
//	
//	private int n = 0;  //  top N features from each document
//	private List<Map<String, Double>> documentWordWeights = new ArrayList<Map<String, Double>>();
//	private Set<String> p = new HashSet<String>();  //  feature pool
//	private int[][]  topNDocumentCounts = null;
//	private String[]  topNWordOrder = null;
//}













