import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

import jsat.distributions.multivariate.*;
import jsat.linear.DenseVector;
import jsat.linear.Vec;

import org.apache.commons.math3.util.CombinatoricsUtils;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;;

//import java.util.*;
//import org.knowceans.util.*;


public class EmClusterModel {
	
	public static void main(String[] args) throws Exception{
		EmClusterModel model = new EmClusterModel();
		
		List<Set<Integer>> documentLabels = new ArrayList<Set<Integer>>();
		TopNFeatureSelector featureSelector = new TopNFeatureSelector();
//		featureSelector.loadDocumentsFromFolder("20_newsgroups", documentLabels);
//		featureSelector.loadDocumentsFromFolder("20_newsgroups_reduced", documentLabels);
//		featureSelector.loadDocumentsFromFolder("20_newsgroups_Test1", documentLabels);
		featureSelector.loadDocumentsFromFolder("20_newsgroups_Test2", documentLabels);

		
		System.out.println("Documents Parsed!");
		
		featureSelector.topNFeatureSelection(2);
		System.out.println("Top features selected");
		
		List<Double> lambda = new ArrayList<Double>();
		List<List<Double>> beta = new ArrayList<List<Double>>();
		final int K = 6;  //  20 clusters because there are 20 subjects

		
		
		
		
//		double[] emARIs = new double[5];
//		double[] cemARIs = new double[5];
//		
//		for(int i = 0; i < 5; i++){
//			model.generateParameters(K, lambda, beta, featureSelector.getTopNDocumentCounts(), featureSelector.getTopNWordOrder());
//			
//			System.out.println("Lambda: " + lambda);
//			model.logLikelihood(lambda, beta, featureSelector.getTopNDocumentCounts());
//		
//			cemARIs[i] = model.cem(lambda, beta, featureSelector.getTopNDocumentCounts(), documentLabels, K);
//			emARIs[i] = model.em(lambda, beta, featureSelector.getTopNDocumentCounts(), documentLabels, K);
//
//		
//					
//			lambda = new ArrayList<Double>();
//			beta = new ArrayList<List<Double>>();
//			
//		}
//		
//		System.out.println("\n\nEM ARIs: " + Arrays.toString(emARIs));
//		System.out.println("Min: " + model.min(emARIs));
//		System.out.println("Max: " + model.max(emARIs));
//		System.out.println("Mean: " + model.mean(emARIs));
//		System.out.println("Median: " + model.median(emARIs));
//		System.out.println("StandardDeviation: " + model.stddev(emARIs));
//		
//		System.out.println("\n\nCEM ARIs: " + Arrays.toString(cemARIs));
//		System.out.println("Min: " + model.min(cemARIs));
//		System.out.println("Max: " + model.max(cemARIs));
//		System.out.println("Mean: " + model.mean(cemARIs));
//		System.out.println("Median: " + model.median(cemARIs));
//		System.out.println("StandardDeviation: " + model.stddev(cemARIs));
		
		model.generateParameters(K, lambda, beta, featureSelector.getTopNDocumentCounts(), featureSelector.getTopNWordOrder());
		
		System.out.println("Lambda: " + lambda);
		model.logLikelihood(lambda, beta, featureSelector.getTopNDocumentCounts());
	
//		model.cem(lambda, beta, featureSelector.getTopNDocumentCounts(), documentLabels, K);
		model.em(lambda, beta, featureSelector.getTopNDocumentCounts(), documentLabels, K);

	
		model.topNWordsForEachClass(lambda, beta, featureSelector.getTopNWordOrder());
		
		
		
		
		
		
		
		
		
//		//  the model seems to work for these smaller data sets with no danger of underflow
//		//  we'll see if the model does well now
//		
//		EmClusterModel model = new EmClusterModel();
//		
//		List<Double> lambda = new ArrayList<Double>();
//		List<List<Double>> beta = new ArrayList<List<Double>>();
//		List<List<Integer>> documentCounts = new ArrayList<List<Integer>>();
//		List<List<String>> documents = new ArrayList<List<String>>();
//		
//		lambda.add(0.3);
//		lambda.add(0.7);
//		
//		List<Double> beta1 = new ArrayList<Double>();
//		beta1.add(0.3);
//		beta1.add(0.7);
//		List<Double> beta2 = new ArrayList<Double>();
//		beta2.add(0.6);
//		beta2.add(0.4);
//		beta.add(beta1);
//		beta.add(beta2);
//		
//		List<String> sample1 = new ArrayList<String>();
//		List<String> sample2 = new ArrayList<String>();
//		List<String> sample3 = new ArrayList<String>();
//		List<String> sample4 = new ArrayList<String>();
//		List<String> sample5 = new ArrayList<String>();
//		sample1.add("H");
//		sample1.add("H");
//		sample1.add("H");
//		sample2.add("T");
//		sample2.add("T");
//		sample2.add("T");
//		sample3.add("H");
//		sample3.add("H");
//		sample3.add("H");
//		sample4.add("T");
//		sample4.add("T");
//		sample4.add("T");
//		sample5.add("H");
//		sample5.add("H");
//		sample5.add("H");
//		documents.add(sample1);
//		documents.add(sample2);
//		documents.add(sample3);
//		documents.add(sample4);
//		documents.add(sample5);
//		
//		
//		List<Integer> sample1c = new ArrayList<Integer>();
//		List<Integer> sample2c = new ArrayList<Integer>();
//		List<Integer> sample3c = new ArrayList<Integer>();
//		List<Integer> sample4c = new ArrayList<Integer>();
//		List<Integer> sample5c = new ArrayList<Integer>();
//		sample1c.add(3);
//		sample1c.add(0);
//		sample2c.add(0);
//		sample2c.add(3);
//		sample3c.add(3);
//		sample3c.add(0);
//		sample4c.add(0);
//		sample4c.add(3);
//		sample5c.add(3);
//		sample5c.add(0);
//		
//		documentCounts.add(sample1c);
//		documentCounts.add(sample2c);
//		documentCounts.add(sample3c);
//		documentCounts.add(sample4c);
//		documentCounts.add(sample5c);
//		
//		
//		
//		//  doesn't correspond to preset betas so flip it
//		//  this will be better when we have to randomly determine betas etc
////		documentCounts = model.documentListToCountMatrix(documents);
//		
//		System.out.println("Lambda: " + lambda);
//		System.out.println("Beta: " + beta);
//		System.out.println("Documents: " + documents);
//		System.out.println("Documents Counts: " + documentCounts);
//		model.logLikelihood(lambda, beta, documentCounts);
//		
//		
//		model.em(lambda, beta, documentCounts);
//		
		
		
		
		
		
		
		
//		EmClusterModel model = new EmClusterModel();
//		
//		
//		List<Double> lambda = new ArrayList<Double>();
//		List<List<Double>> beta = new ArrayList<List<Double>>();
//		List<List<Integer>> documentCounts = new ArrayList<List<Integer>>();
//		List<List<String>> documents = new ArrayList<List<String>>();
//		
//		lambda.add(0.25); //  cluster1
//		lambda.add(0.25);  //  cluster2
//		lambda.add(0.25);
//		lambda.add(0.25);
//		
//		
//		
////		List<Double> beta1 = new ArrayList<Double>();
////		beta1.add(0.3);  //  how
////		beta1.add(0.1);  //  it
////		beta1.add(0.3);  //  ball
////		beta1.add(0.3);  //  try
////		List<Double> beta2 = new ArrayList<Double>();
////		beta2.add(0.1);
////		beta2.add(0.1);
////		beta2.add(0.1);
////		beta2.add(0.7);
////		List<Double> beta3 = new ArrayList<Double>();
////		beta3.add(0.3);  //  how
////		beta3.add(0.2);  //  it
////		beta3.add(0.3);  //  ball
////		beta3.add(0.2);  //  try
////		List<Double> beta4 = new ArrayList<Double>();
////		beta4.add(0.2);  //  how
////		beta4.add(0.2);  //  it
////		beta4.add(0.1);  //  ball
////		beta4.add(0.5);  //  try
////		beta.add(beta1);
////		beta.add(beta2);
////		beta.add(beta3);
////		beta.add(beta4);
//		beta.add(model.generateUninformedDirichletDistributionList(4));
//		beta.add(model.generateUninformedDirichletDistributionList(4));
//		beta.add(model.generateUninformedDirichletDistributionList(4));
//		beta.add(model.generateUninformedDirichletDistributionList(4));
//		
//		
//		
//		List<String> sample1 = new ArrayList<String>();
//		List<String> sample2 = new ArrayList<String>();
//		List<String> sample3 = new ArrayList<String>();
//		List<String> sample4 = new ArrayList<String>();
//		List<String> sample5 = new ArrayList<String>();
//		List<String> sample6 = new ArrayList<String>();
//		sample1.add("Try");
//		sample1.add("Try");
//		sample1.add("Try");
//		sample1.add("Try");
//		sample1.add("Try");
//		sample2.add("Ball");
//		sample2.add("Ball");
//		sample2.add("It");
//		sample2.add("It");
//		sample2.add("Try");
//		sample3.add("How");
//		sample3.add("It");
//		sample3.add("How");
//		sample3.add("It");
//		sample3.add("Try");
//		sample4.add("How");
//		sample4.add("Try");
//		sample4.add("How");
//		sample4.add("Try");
//		sample4.add("Ball");
//		sample5.add("It");
//		sample5.add("Try");
//		sample5.add("It");
//		sample5.add("Try");
//		sample5.add("Try");
//		sample6.add("It");
//		sample6.add("Try");
//		sample6.add("It");
//		sample6.add("Try");
//		sample6.add("Try");
//		documents.add(sample1);
//		documents.add(sample2);
//		documents.add(sample3);
//		documents.add(sample4);
//		documents.add(sample5);
//		documents.add(sample6);
//		
//		
//		List<Integer> sample1c = new ArrayList<Integer>();
//		List<Integer> sample2c = new ArrayList<Integer>();
//		List<Integer> sample3c = new ArrayList<Integer>();
//		List<Integer> sample4c = new ArrayList<Integer>();
//		List<Integer> sample5c = new ArrayList<Integer>();
//		List<Integer> sample6c = new ArrayList<Integer>();
//		sample1c.add(0);
//		sample1c.add(0);
//		sample1c.add(0);
//		sample1c.add(5);
//		sample2c.add(0);
//		sample2c.add(2);
//		sample2c.add(2);
//		sample2c.add(1);
//		sample3c.add(2);
//		sample3c.add(2);
//		sample3c.add(0);
//		sample3c.add(1);
//		sample4c.add(2);
//		sample4c.add(0);
//		sample4c.add(1);
//		sample4c.add(2);
//		sample5c.add(0);
//		sample5c.add(2);
//		sample5c.add(0);
//		sample5c.add(3);
//		sample6c.add(0);
//		sample6c.add(2);
//		sample6c.add(0);
//		sample6c.add(3);
//		
//		documentCounts.add(sample1c);
//		documentCounts.add(sample2c);
//		documentCounts.add(sample3c);
//		documentCounts.add(sample4c);
//		documentCounts.add(sample5c);
//		documentCounts.add(sample6c);
//		
//		
//		
//		//  doesn't correspond to preset betas so flip it
//		//  this will be better when we have to randomly determine betas etc
////		documentCounts = model.documentListToCountMatrix(documents);
//		
//		System.out.println("Lambda: " + lambda);
//		System.out.println("Beta: " + beta);
//		System.out.println("Documents: " + documents);
//		System.out.println("Documents Counts: " + documentCounts);
//		model.logLikelihood(lambda, beta, documentCounts);
//		
//		
//		model.em(lambda, beta, documentCounts);
		
	
		
		
	}
	
	public EmClusterModel(){	}
	
	
	
	
	private void topNWordsForEachClass(List<Double> lambda, List<List<Double>>beta, List<String> topNWordOrder){
		for(int i = 0; i < lambda.size(); i++){
			System.out.println("\n\nTop 2 words for class " + i + ":");
			
//			double denomiator = 0.0;
//			for(int k = 0; k < beta.get(i).size(); k++){
//				denomiator += ((double)beta.get(i).get(k) * (double)lambda.get(i));
//			}
			
			double highestValue = Double.NEGATIVE_INFINITY;
			int highestIndex = -1;
			double secondHighestValue = Double.NEGATIVE_INFINITY;
			int secondHighestIndex = -1;
			for(int k = 0; k < beta.get(i).size(); k++){
				double newvalue = (double)beta.get(i).get(k) * (double)lambda.get(i);
				if(newvalue > highestValue){
					highestValue = newvalue;
					highestIndex = k;
				}
			}
			
			for(int k = 0; k < beta.get(i).size(); k++){
				double newvalue = (double)beta.get(i).get(k) * (double)lambda.get(i);
				if(newvalue != highestValue && k != highestIndex && newvalue > secondHighestValue){
					secondHighestValue = newvalue;
					secondHighestIndex = k;
				}
			}
			
			System.out.println(topNWordOrder.get(highestIndex) + ", " + topNWordOrder.get(secondHighestIndex));
			
			
		}
	}
	
	
	
	
	
	 /**
     * Returns the maximum value in the array a[], -infinity if no such value.
     */
   public double max(double[] a) {
       double max = Double.NEGATIVE_INFINITY;
       for (int i = 0; i < a.length; i++) {
           if (Double.isNaN(a[i])) return Double.NaN;
           if (a[i] > max) max = a[i];
       }
       return max;
   }


  /**
    * Returns the minimum value in the array a[], +infinity if no such value.
    */
   public double min(double[] a) {
       double min = Double.POSITIVE_INFINITY;
       for (int i = 0; i < a.length; i++) {
           if (Double.isNaN(a[i])) return Double.NaN;
           if (a[i] < min) min = a[i];
       }
       return min;
   }


  /**
    * Returns the average value in the array a[], NaN if no such value.
    */
   public double mean(double[] a) {
       if (a.length == 0) return Double.NaN;
       double sum = sum(a);
       return sum / a.length;
   }
   
   /**
    * Returns the sum of all values in the array a[].
    */
   public double sum(double[] a) {
       double sum = 0.0;
       for (int i = 0; i < a.length; i++) {
           sum += a[i];
       }
       return sum;
   }
   
   /**
    * Returns the median value of the passed array
    * @param a
    * @return
    */
   public double median(double[] a){
	   Arrays.sort(a);
	   
	   return a[a.length / 2];
   }

  /**
    * Returns the sample variance in the array a[], NaN if no such value.
    */
   public double var(double[] a) {
       if (a.length == 0) return Double.NaN;
       double avg = mean(a);
       double sum = 0.0;
       for (int i = 0; i < a.length; i++) {
           sum += (a[i] - avg) * (a[i] - avg);
       }
       return sum / (a.length - 1);
   }


  /**
    * Returns the sample standard deviation in the array a[], NaN if no such value.
    */
   public double stddev(double[] a) {
       return Math.sqrt(var(a));
   }
	
	
	
	
	
	
	
	/**
	 * This will be used later when I need to generate uniform parameters
	 * or use a dirichlet prior or something like that. I'm just going to hard-code all this for now
	 */
	private void generateParameters(final int K, List<Double> lambda, List<List<Double>> beta, 
			List<List<Integer>> documentCounts, List<String> documentWordOrder){

		
		List<Double> tmpLambda = generateUniformlyDistributedList(K);  //  K values of lambda for K clusters
		for(int i = 0; i < tmpLambda.size(); i++){
			lambda.add(tmpLambda.get(i));
		}
		
		
		
		for(int i = 0; i < K; i++){  //  for each class
			//  generate a uniform list of probabilities across all words in the vocabulary
//			beta.add(generateUniformlyDistributedList(documentWordOrder.size()));
			beta.add(generateUninformedDirichletDistributionList(documentWordOrder.size()));
			
		}
	}
	
	/**
	 * Generates a uniform list (all of their probabilities together sum to 1) of passed size
	 * @param size
	 * @return the uniform list
	 */
	private List<Double> generateUniformlyDistributedList(int size){
		List<Double> uniformList = new ArrayList<Double>();
		
		double uniformProbability = 1.0 / (double)size;
		//  make sure list will exactly sum to 1
		double lastProbability = 1.0 - (uniformProbability * (size - 1)); 
		for(int i = 0; i < size - 1; i++){
			uniformList.add(uniformProbability);
		}
		uniformList.add(lastProbability);
		
		return uniformList;
	}
	
	
	/**
	 * Should generate values taken from a unifrom dirichlet distribution
	 * @param size
	 * @return
	 */
	private List<Double> generateUninformedDirichletDistributionList(int size){
		
		List<Double> alphas = new ArrayList<Double>();
		for(int i = 0; i < size; i++){
			alphas.add(0.1);
		}
		Dirichlet d = new Dirichlet(new DenseVector(alphas));
		List<Vec> vectorSamples = d.sample(1, new Random(System.nanoTime()));
		
		//  put vector double sample values into a usable list of doubles
		List<Double> samples = new ArrayList<Double>();
//		double cumulativeProbability = 0.0;
		for(int i = 0; i < size; i++){
			samples.add(vectorSamples.get(0).get(i));
//			cumulativeProbability += vectorSamples.get(0).get(i);
		}
//		samples.add(1 - cumulativeProbability);
//		System.out.println(cumulativeProbability);
		
		
		return samples;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//  TODO:  need to do all of this in logspace
	//  TODO:   add 1 smoothing 
	//  TODO:   don't output partial counts after I've debugged
	//  DONE:   random dirichlet distribution
	//  DONE:  this em can handle multiple cluster classification, it works perfectly with 4 classes
	
	
	
	
	
	
	/**
	 * Parameters these must all be array lists to ensure speed, linked lists not allowed.
	 * (Indexes in Beta correspond to word category indexes found in documentCounts)
	 * @param lambda  
	 * @param beta
	 * @param documentCounts
	 */
	private double em(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts, List<Set<Integer>> documentLabels, int k){
		if(lambda.size() != beta.size() || beta.size() <= 0 || documentCounts.size() <= 0 || beta.get(0).size() != documentCounts.get(0).size()){
			//  do not output just return becaue the parameters are not correct
			System.out.println("EM PARAMETERS ARE INCORRECT!");
			return 0.0;
		}
		
//		List<List<Double>> documentCounts = convertToLogSpace(lambda, beta, documentCountsInt);

		double oldLogLikelihood = Double.NEGATIVE_INFINITY, newLogLikelihood = Double.NEGATIVE_INFINITY;
		List<List<Double>> finalPartialCounts = new ArrayList<List<Double>>();
//		for(int i = 0; i < hardcodedIterationNumber; i++){
		for(int i = 0;; i++){

			System.out.println("\n\nIteration #" + i + ":");
			System.out.println("Lambda: " + lambda);
			System.out.println("Beta: " + beta);
			
			
			newLogLikelihood = logLikelihood(lambda, beta, documentCounts);
			
			
			
			finalPartialCounts = emIterationLog(lambda, beta, documentCounts);
			
			//  if percent changed on log likelihood is less then 3.5%
			double percentChange  = Math.abs(((Math.abs(newLogLikelihood) - Math.abs(oldLogLikelihood)) / Math.abs(oldLogLikelihood)));
			if(oldLogLikelihood != Double.NEGATIVE_INFINITY && percentChange < 0.00035){
				break;
			}
			oldLogLikelihood = newLogLikelihood;
			
		}
		
		System.out.println("\nCluster Assignments:");
		return clusterAssignment(finalPartialCounts, documentLabels, k);
		
		
		
	}
	
	
	
	/**
	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
	 * It goes through each class computing the partial counts for each word along the way and computers the
	 * new lambda's and beta array for each class as it finished it's pass through every document
	 * @param lambda
	 * @param beta
	 * @param documentCounts
	 */
	private List<List<Double>> emIteration(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
		List<Double> newLambda = new ArrayList<Double>(); 
		List<List<Double>> newBeta = new ArrayList<List<Double>>();
		List<List<Double>> partialCountsByClass = new ArrayList<List<Double>>();
		
		List<Double> denominators = new ArrayList<Double>();
		//  for every class Cn
		for(int i = 0; i < lambda.size(); i++){
			
			List<Double> partialCountsForClass = new ArrayList<Double>();
			//  for every document DOCn
			for(int j = 0; j < documentCounts.size(); j++){
				
				//  P(C=ci)
				double partialCountTotal = lambda.get(i);
				//  for every word in the vocabulary
				for(int k = 0; k < beta.get(i).size(); k++){
					//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
					double aa = beta.get(i).get(k);
					double bb = documentCounts.get(j).get(k);
					partialCountTotal *= Math.pow(aa, bb); 
				}
				
				//  this will happen once for each document, denominators are the same only
				//  across the same class
				if(denominators.size() <= j){
					double tempDenominator = 0.0;
					//  do above for every class given the document's counts and sum together for denominator
					for(int c = 0; c < lambda.size(); c++){
						double tempPartialCountForClass = lambda.get(c);
						//  for every word in the vocabulary
						for(int k = 0; k < beta.get(c).size(); k++){
							//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
							tempPartialCountForClass *= Math.pow(beta.get(c).get(k), documentCounts.get(j).get(k)); 
						}
						tempDenominator += tempPartialCountForClass;
					}
					denominators.add(tempDenominator);
				}
				
				partialCountsForClass.add(partialCountTotal / denominators.get(j));
			}
			//  save each partial count to assign clusters after
			partialCountsByClass.add(partialCountsForClass);
			
//			System.out.println("PartialCounts for class " + i + ": " + partialCountsForClass);
			
			
			double totalCountForAllWordsForClass = 0.0;
			//  sum counts of every word given the class Ci
			for(int k = 0; k < beta.get(i).size(); k++){
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals up all words counts for class Ci
					totalCountForAllWordsForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
				}
			}
			
			List<Double> betaValuesForClass = new ArrayList<Double>();
			double totalCountForOneWordForClass;
			for(int k = 0; k < beta.get(i).size(); k++){
				totalCountForOneWordForClass = 0.0;
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals accross all documents just for this word
					totalCountForOneWordForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
				}
				//  nextValueAddedTo Betas for this class
				betaValuesForClass.add(totalCountForOneWordForClass / totalCountForAllWordsForClass);
			}
			newBeta.add(betaValuesForClass);
			
			//  Lambda value for class is simple all the partial counts for the class added together 
			double lambaValueForClass = 0.0;
			for(int l = 0; l < partialCountsForClass.size(); l++){
				lambaValueForClass += partialCountsForClass.get(l);
			}
			newLambda.add(lambaValueForClass / (double)documentCounts.size());
		
			
			
		}
		
		
		//  set the original lambda and beta to the new values
		for(int i = 0; i < lambda.size(); i++){
			lambda.set(i, newLambda.get(i));
		}
		for(int i = 0; i < beta.size(); i++){
			for(int j = 0; j < beta.get(i).size(); j++){
				beta.get(i).set(j, newBeta.get(i).get(j));
			}
		}
		
		
		return partialCountsByClass;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
	 * It goes through each class computing the partial counts for each word along the way and computers the
	 * new lambda's and beta array for each class as it finished it's pass through every document
	 * @param lambda
	 * @param beta
	 * @param documentCounts
	 */
	private List<List<Double>> emIterationLog(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
		List<Double> newLambda = new ArrayList<Double>(); 
		List<List<Double>> newBeta = new ArrayList<List<Double>>();
		List<List<Double>> partialCountsByClass = new ArrayList<List<Double>>();
		
		List<Double> denominators = new ArrayList<Double>();
		//  for every class Cn
		for(int i = 0; i < lambda.size(); i++){
			
			List<Double> partialCountsForClass = new ArrayList<Double>();
			//  for every document DOCn
			for(int j = 0; j < documentCounts.size(); j++){
				
				//  log(P(C=ci))
				double partialCountTotal = Math.log(lambda.get(i));
				//  for every word in the vocabulary
				for(int k = 0; k < beta.get(i).size(); k++){
					//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
					double aa = beta.get(i).get(k);
					double bb = documentCounts.get(j).get(k);
					partialCountTotal += (bb * Math.log(aa)); 
				}
				
				//  this will happen once for each document, denominators are the same only
				//  across the same class
				if(denominators.size() <= j){
					List<Double> tempDenominator = new ArrayList<Double>();
					//  do above for every class given the document's counts and sum together for denominator
					for(int c = 0; c < lambda.size(); c++){
						double tempPartialCountForClass = Math.log(lambda.get(c));
						//  for every word in the vocabulary
						for(int k = 0; k < beta.get(c).size(); k++){
							//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
							tempPartialCountForClass += (documentCounts.get(j).get(k) * Math.log(beta.get(c).get(k))); 
						}
						tempDenominator.add(tempPartialCountForClass);
					}
					denominators.add(logSum(tempDenominator));
				}
				
				partialCountsForClass.add(partialCountTotal - denominators.get(j));
			}
			//  save each partial count to assign clusters after
			partialCountsByClass.add(partialCountsForClass);
			
//			System.out.println("PartialCounts for class " + i + ": " + partialCountsForClass);
			
			
			//  covertPartialCountsBackIntoNormalSpace by exponentiating each value separately
			for(int p = 0; p < partialCountsForClass.size(); p++){
				partialCountsForClass.set(p, Math.exp(partialCountsForClass.get(p)));
			}
			
			double totalCountForAllWordsForClass = 0.0;
			//  sum counts of every word given the class Ci
			for(int k = 0; k < beta.get(i).size(); k++){
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals up all words counts for class Ci
					totalCountForAllWordsForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
				}
			}
			
			List<Double> betaValuesForClass = new ArrayList<Double>();
			double totalCountForOneWordForClass;
			for(int k = 0; k < beta.get(i).size(); k++){
				totalCountForOneWordForClass = 0.0;
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals accross all documents just for this word
					totalCountForOneWordForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
				}
				//  nextValueAddedTo Betas for this class
				//  add one smoothing for betas
				betaValuesForClass.add(((double)totalCountForOneWordForClass + 1.0) / ((double)totalCountForAllWordsForClass + (double)beta.get(0).size()));
			}
			newBeta.add(betaValuesForClass);
			
			//  Lambda value for class is simple all the partial counts for the class added together 
			double lambaValueForClass = 0.0;
			for(int l = 0; l < partialCountsForClass.size(); l++){
				lambaValueForClass += partialCountsForClass.get(l);
			}
			//  should not smooth lambda, it yeilds incorrect values
//			newLambda.add(((double)lambaValueForClass + 1.0) / ((double)documentCounts.size() + (double)beta.get(0).size()));
			newLambda.add(((double)lambaValueForClass) / ((double)documentCounts.size()));

			
			
		}
		
		
		//  set the original lambda and beta to the new values
		for(int i = 0; i < lambda.size(); i++){
			lambda.set(i, newLambda.get(i));
		}
		for(int i = 0; i < beta.size(); i++){
			for(int j = 0; j < beta.get(i).size(); j++){
				beta.get(i).set(j, newBeta.get(i).get(j));
			}
		}
		
		
		return partialCountsByClass;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * Parameters these must all be array lists to ensure speed, linked lists not allowed.
	 * (Indexes in Beta correspond to word category indexes found in documentCounts)
	 * @param lambda  
	 * @param beta
	 * @param documentCounts
	 */
	private double cem(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts, List<Set<Integer>> documentLabels,  int k){
		if(lambda.size() != beta.size() || beta.size() <= 0 || documentCounts.size() <= 0 || beta.get(0).size() != documentCounts.get(0).size()){
			//  do not output just return becaue the parameters are not correct
			System.out.println("EM PARAMETERS ARE INCORRECT!");
			return 0.0;
		}
		

		double oldLogLikelihood = Double.NEGATIVE_INFINITY, newLogLikelihood = Double.NEGATIVE_INFINITY;
//		List<List<Double>> finalPartialCounts = new ArrayList<List<Double>>();
		List<Integer> highestClassIndexForDocument = new ArrayList<Integer>();
//		for(int i = 0; i < hardcodedIterationNumber; i++){
		for(int i = 0;; i++){

			System.out.println("\n\nIteration #" + i + ":");
			System.out.println("Lambda: " + lambda);
			System.out.println("Beta: " + beta);
			
			
			newLogLikelihood = logLikelihood(lambda, beta, documentCounts);
			
			highestClassIndexForDocument = cemIterationLog(lambda, beta, documentCounts);
			
			//  if percent changed on log likelihood is less then 3.5%
			double percentChange  = Math.abs(((Math.abs(newLogLikelihood) - Math.abs(oldLogLikelihood)) / Math.abs(oldLogLikelihood)));
			if(oldLogLikelihood != Double.NEGATIVE_INFINITY && percentChange < 0.00035){
				break;
			}
			oldLogLikelihood = newLogLikelihood;
			
		}
		
		System.out.println("\nCluster Assignments:");
		return clusterAssignmentCem(highestClassIndexForDocument, documentLabels, k);
		
	}
	
	
	
	
	/**
	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
	 * It goes through each class computing the partial counts for each word along the way and computers the
	 * new lambda's and beta array for each class as it finished it's pass through every document
	 * @param lambda
	 * @param beta
	 * @param documentCounts
	 */
	private List<Integer> cemIterationLog(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
		List<Double> newLambda = new ArrayList<Double>(); 
		List<List<Double>> newBeta = new ArrayList<List<Double>>();
		
		//  This saves the index and value of the highest class partial count for each document
		List<Integer> highestClassIndexForDocument = new ArrayList<Integer>();
		List<Double> highestClassPartialCountForDocument = new ArrayList<Double>();
		for(int i = 0; i < documentCounts.size(); i++){
			highestClassIndexForDocument.add(-1);
			highestClassPartialCountForDocument.add(Double.NEGATIVE_INFINITY);
		}
		
		//  for every class Cn
		for(int i = 0; i < lambda.size(); i++){
			
			//  for every document DOCn
			for(int j = 0; j < documentCounts.size(); j++){
				
				//  log(P(C=ci))
				double partialCountTotal = Math.log(lambda.get(i));
				//  for every word in the vocabulary
				for(int k = 0; k < beta.get(i).size(); k++){
					//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
					double aa = beta.get(i).get(k);
					double bb = documentCounts.get(j).get(k);
					partialCountTotal += (bb * Math.log(aa)); 
				}
				
				if(partialCountTotal > highestClassPartialCountForDocument.get(j)){
					highestClassPartialCountForDocument.set(j, partialCountTotal);
					highestClassIndexForDocument.set(j, i);
				}
				
			}
			//  save each partial count to assign clusters after, TODO: might need to fix this
			
		}
		
		
		
		//  set all lambda according to list generated from cem
		for(int i = 0; i < lambda.size(); i++){
			int classCount = 0;
			for(int j = 0; j < highestClassIndexForDocument.size(); j++){
				if(highestClassIndexForDocument.get(j) == i){
					classCount++;
				}
			}
			newLambda.add((double)classCount / (double)highestClassIndexForDocument.size());
		}
		
		
		
	
		
		//  go through every class
		for(int c = 0; c < lambda.size(); c++){
			
			
			
			double totalCountForAllWordsForClass = 0.0;
			//  sum counts of every word given the class Ci
			for(int k = 0; k < beta.get(0).size(); k++){
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals up all words counts for class Ci
					if(highestClassIndexForDocument.get(d) == c){
						totalCountForAllWordsForClass += (documentCounts.get(d).get(k));
					}
				}
			}
			
			
			
			
			
			
			List<Double> betaValuesForClass = new ArrayList<Double>();
			double totalCountForOneWordForClass;
			for(int k = 0; k < beta.get(0).size(); k++){
				totalCountForOneWordForClass = 0.0;
				for(int d = 0; d < documentCounts.size(); d++){
					//  this totals across all documents just for this word
					
					
					if(highestClassIndexForDocument.get(d) == c){
						totalCountForOneWordForClass += (documentCounts.get(d).get(k));
					}
				}
				//  nextValueAddedTo Betas for this class
				//  add one smoothing for betas
				betaValuesForClass.add(((double)totalCountForOneWordForClass + 1.0) / ((double)totalCountForAllWordsForClass + (double)beta.get(0).size()));
			}
			newBeta.add(betaValuesForClass);
		}

		
		
		
		
		
		
		
		
		//  set the original lambda and beta to the new values
		for(int i = 0; i < lambda.size(); i++){
			lambda.set(i, newLambda.get(i));
		}
		for(int i = 0; i < beta.size(); i++){
			for(int j = 0; j < beta.get(i).size(); j++){
				beta.get(i).set(j, newBeta.get(i).get(j));
			}
		}
		
		
		return highestClassIndexForDocument;
	}
	
	
	
	
	
	
//	/**
//	 * Converts given lambda, beta, and documentCounts into logspace
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 * @return document counts using doubles instead of ints since it's it logspace
//	 */
//	private List<List<Double>> convertToLogSpace(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCountsInt){
//		for(int i = 0; i < beta.size(); i++){
//			lambda.set(i, Math.log(lambda.get(i)));
//			
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, Math.log(beta.get(i).get(j)));
//			}
//		}
//		
//		List<List<Double>> documentCountsDouble = new ArrayList<List<Double>>();
//		for(int i = 0; i < documentCountsInt.size(); i++){
//			
//			List<Double> doubleCounts = new ArrayList<Double>();
//			for(int j = 0; j < documentCountsInt.get(i).size(); j++){
//				doubleCounts.add(Math.log(documentCountsInt.get(i).get(j)));
//				
//			}
//			documentCountsDouble.add(doubleCounts);
//		}
//		return documentCountsDouble;
//	}
	
	
	
	
	/**
	 * Adds to logs already in log space and returns an answer still in log space
	 * @param logX
	 * @param logY
	 * @return
	 */
	 private double logAdd(double logX, double logY){
	     // 1. make X the max
	     if (logY > logX){
	         double temp = logX;
	         logX = logY;
	         logY = temp;
	     }
	     // 2. now X is bigger
	     if (logX == Double.NEGATIVE_INFINITY){
	         return logX;
	     }
	     // 3. how far "down" (think decibels) is logY from logX?
	     //    if it's really small (20 orders of magnitude smaller), then ignore
	     double negDiff = logY - logX;
	     if (negDiff < -20){
	         return logX;
	     }
	     // 4. otherwise use some nice algebra to stay in the log domain
	     //    (except for negDiff)
	     return logX + Math.log(1.0 + Math.exp(negDiff));
	 }
	
	
	/**
	 * Takes a list of logarithmic quantities
	 * @param L
	 * @return
	 */
	 private double logSum(List<Double> L){
		double logResult = Double.NEGATIVE_INFINITY;
		for(double logX : L){
			logResult = logAdd(logResult, logX);
		}
		return logResult;
	}
	
	
	
	
	/**
	 * Calculates the logLikelihood of the current iteration
	 * @param lambda
	 * @param beta
	 * @param documentCounts
	 */
	private double logLikelihood(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
		double logLikelihood = 0.0;
		List<Integer> totalWordsByDocument = new ArrayList<Integer>();
		
		
		for(int i = 0; i < documentCounts.size(); i++){
			int tempSum = 0;  // sum up the total number of words for each document
			for(int j = 0; j < documentCounts.get(0).size(); j++){
				tempSum += documentCounts.get(i).get(j);
			}
			totalWordsByDocument.add(tempSum);  //  then add this total at the document's position in totalWords..
		}
		
		//  sum across all documents
		for(int i = 0; i < documentCounts.size(); i++){
			List<Double> tempValuesByClass = new ArrayList<Double>();
			//  find intermediary values for each class
			for(int c = 0; c < lambda.size(); c++){
				double classValue = 0.0;
				classValue += Math.log(lambda.get(c));
//				classValue += Math.log(factorial(totalWordsByDocument.get(i)));
				classValue += CombinatoricsUtils.factorialLog(totalWordsByDocument.get(i));
				
				double sumAcrossVocabulary = 0.0;
				for(int j = 0; j < documentCounts.get(i).size(); j++){
					sumAcrossVocabulary += (documentCounts.get(i).get(j) * Math.log(beta.get(c).get(j)));
//					sumAcrossVocabulary -= Math.log(factorial(documentCounts.get(i).get(j)));
					sumAcrossVocabulary -= CombinatoricsUtils.factorialLog(documentCounts.get(i).get(j));
				}
				classValue += sumAcrossVocabulary;
				tempValuesByClass.add(classValue);
			}
			
			
			//  this only works when there's no risk of serious underflow
//			//  this part is log[e^a+.....+e^n] once a-n class values have all been found
//			double innerExpSum = 0.0;
//			for(int c = 0; c < tempValuesByClass.size(); c++){
//				innerExpSum += Math.pow(Math.E, tempValuesByClass.get(c));
//			}
//			//  add this value to the sum across all documents
//			logLikelihood += Math.log(innerExpSum);
			
			List<Double> logSumTerms = new ArrayList<Double>();
			for(int c = 0; c < tempValuesByClass.size(); c++){
				logSumTerms.add(tempValuesByClass.get(c));
			}
			logLikelihood += logSum(logSumTerms);
		}
		
		
		System.out.println("Log-Likelihood:" + logLikelihood);
		return logLikelihood;
	}
	 
	 
	
	/**
	 * Returns the factorial of passed number
	 * @param n
	 * @return
	 */
	private double factorial(double n){
		double result = 1.0;
		for (double i = 1.0; i <= n; i++) {
		   result = result * i;
		}
		return result;
	}
	
	/**
	 * Takes the final partial counts and uses them to assign each document
	 * To a specific cluster
	 * @param partialCounts
	 */
	private double clusterAssignment(List<List<Double>> partialCounts, List<Set<Integer>> documentLabels, int k){
		int classAssignment = -1;
		List<Set<Integer>> clusterAssignments = new ArrayList<Set<Integer>>();
		for(int i = 0; i < k; i++){
			clusterAssignments.add(new HashSet<Integer>());
		}
		

		double maxValue = Double.NEGATIVE_INFINITY;
		if(partialCounts.size() <= 0)
			return 0.0;
		
		//  go through each document
		for(int i = 0; i < partialCounts.get(0).size(); i++){
			//  go through each class for this document
			for(int j = 0; j < partialCounts.size(); j++){
				if(partialCounts.get(j).get(i) > maxValue){
					classAssignment = j;
					maxValue = partialCounts.get(j).get(i);
				}
			}
			//  assign each value to a cluster
			System.out.println("Document " + i + " assigned to cluster " + classAssignment);
			clusterAssignments.get(classAssignment).add(i);
			maxValue = Double.NEGATIVE_INFINITY;
			classAssignment = -1;
		}
		
		
		double ariValue = calculateARI(documentLabels, clusterAssignments, k);
		return ariValue;
	}
	

	
	/**
	 * Outputs the cluster assignments for each document for cem
	 * @param highestClassIndexForDocument
	 */
	private double clusterAssignmentCem(List<Integer> highestClassIndexForDocument, List<Set<Integer>> documentLabels,  int k){
		List<Set<Integer>> clusterAssignments = new ArrayList<Set<Integer>>();
		for(int i = 0; i < k; i++){
			clusterAssignments.add(new HashSet<Integer>());
		}
		
		for(int i = 0; i < highestClassIndexForDocument.size(); i++){
			System.out.println("Document " + i + " assigned to " + highestClassIndexForDocument.get(i));
			clusterAssignments.get(highestClassIndexForDocument.get(i)).add(i);
		}
		
		
		double ariValue = calculateARI(documentLabels, clusterAssignments, k);
		return ariValue;
	}
	
	
	/**
	 * Calculates ARI value for the two passed clusterings, also outputs the contingency
	 * matrix used to find the ARI value
	 * @param documentLabels
	 * @param clusterAssignments
	 * @param k
	 * @return ARI value
	 */
	private double calculateARI(List<Set<Integer>> documentLabels, List<Set<Integer>> clusterAssignments, int k){
		double ariValue = 0.0;
		List<List<Integer>> contigencyMatrix = new ArrayList<List<Integer>>(); 
		for(int i = 0; i < k + 1; i++){
			List<Integer> labeledRow = new ArrayList<Integer>(); 
			for(int j = 0; j < k + 1; j++){
				labeledRow.add(0);
			}
			contigencyMatrix.add(labeledRow);
		}
		
		//  find intersections among all clusters in the contingency matrix
		for(int i = 0; i < k; i++){
			for(int j = 0; j < k; j++){
				Set<Integer> tmpEM = new HashSet<Integer>(clusterAssignments.get(i));
				Set<Integer> tmpGS = new HashSet<Integer>(documentLabels.get(j));
				tmpEM.retainAll(tmpGS);  //  intersection of the two
				contigencyMatrix.get(i).set(j, tmpEM.size());
			}
		}
		
		//   sum the rows
		for(int i = 0; i < k; i++){
			int rowSum = 0;
			for(int j = 0; j < k; j++){
				rowSum += contigencyMatrix.get(i).get(j);
			}
			contigencyMatrix.get(i).set(k, rowSum);
		}
		
		//   sum the columns
		for(int i = 0; i < k + 1; i++){
			int colSum = 0;
			for(int j = 0; j < k; j++){
				colSum += contigencyMatrix.get(j).get(i);
			}
			contigencyMatrix.get(k).set(i, colSum);
		}
		
		for(int cm = 0; cm < contigencyMatrix.size(); cm++){
			System.out.println(contigencyMatrix.get(cm));
		}
		
		
		double a = 0.0;  //  goes through every count
		for(int i = 0; i < k; i++){
			for(int j = 0; j < k; j++){
				a += mChooseN(contigencyMatrix.get(i).get(j), 2);
			}
		}
		
		//  goes through every rowSum and colSum
		double colSumValue = 0.0;
		for(int j = 0; j < k; j++){
			colSumValue += mChooseN(contigencyMatrix.get(k).get(j), 2);
		}
		double rowSumValue = 0.0;
		for(int j = 0; j < k; j++){
			rowSumValue += mChooseN(contigencyMatrix.get(j).get(k), 2);
		}
		
		double b = colSumValue * rowSumValue;
		
		double c = (mChooseN(contigencyMatrix.get(k).get(k), 2));
		
		double d = colSumValue + rowSumValue;
		
		
		ariValue = (((double)a - ((double)b / (double)c)) / (((double)d / 2.0) - ((double)b / (double)c)));
		
		
		System.out.println("ARI: " + ariValue);
		
		return ariValue;
	}
	
	
	/**
	 * Simple n choose m function
	 * @param m
	 * @param n
	 * @return the result of m choose n as a double
	 */
	private double mChooseN(int m, int n){
		if(m < n || n < 0){
			return 0.0;
		}
		
		//  compute in logspace
		double a = CombinatoricsUtils.factorialLog(m);
		double b = CombinatoricsUtils.factorialLog(n) + CombinatoricsUtils.factorialLog(m - n);
		
		return (Math.exp(((double)a - (double)b)));
	}
	
}





























































































































//
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Random;
//
//import jsat.distributions.multivariate.*;
//import jsat.linear.DenseVector;
//import jsat.linear.Vec;
//import org.apache.commons.math3.util.CombinatoricsUtils;
//
////import java.util.*;
////import org.knowceans.util.*;
//
//
//public class EmClusterModel {
//	
//	public static void main(String[] args) throws Exception{
//		EmClusterModel model = new EmClusterModel();
//		
//		List<Integer> documentLabels = new ArrayList<Integer>();
//		TopNFeatureSelector featureSelector = new TopNFeatureSelector();
//		featureSelector.loadDocumentsFromFolder("20_newsgroups", documentLabels);
////		featureSelector.loadDocumentsFromFolder("20_newsgroups_Test1", documentLabels);
//
//		System.out.println("Documents Parsed!");
//		
//		featureSelector.topNFeatureSelection(2);
//		System.out.println("Top features selected");
//		
//		
//		
//		List<Double> lambda = new ArrayList<Double>();
//		List<List<Double>> beta = new ArrayList<List<Double>>();
//		final int K = 3;  //  20 clusters because there are 20 subjects
//
//		System.out.println(featureSelector.getTopNDocumentCounts());
//		System.out.println(featureSelector.getTopNWordOrder());
//		
//		model.generateParameters(K, lambda, beta, featureSelector.getTopNDocumentCounts(), featureSelector.getTopNWordOrder());
//		
//		System.out.println("Lambda: " + lambda);
//		
//		
//		model.logLikelihood(lambda, beta, featureSelector.getTopNDocumentCounts());
//	
//		
//		
//		model.cem(lambda, beta, featureSelector.getTopNDocumentCounts());
////		model.em(lambda, beta, featureSelector.getTopNDocumentCounts(), documentLabels, K);
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
////		//  the model seems to work for these smaller data sets with no danger of underflow
////		//  we'll see if the model does well now
////		
////		EmClusterModel model = new EmClusterModel();
////		
////		List<Double> lambda = new ArrayList<Double>();
////		List<List<Double>> beta = new ArrayList<List<Double>>();
////		List<List<Integer>> documentCounts = new ArrayList<List<Integer>>();
////		List<List<String>> documents = new ArrayList<List<String>>();
////		
////		lambda.add(0.3);
////		lambda.add(0.7);
////		
////		List<Double> beta1 = new ArrayList<Double>();
////		beta1.add(0.3);
////		beta1.add(0.7);
////		List<Double> beta2 = new ArrayList<Double>();
////		beta2.add(0.6);
////		beta2.add(0.4);
////		beta.add(beta1);
////		beta.add(beta2);
////		
////		List<String> sample1 = new ArrayList<String>();
////		List<String> sample2 = new ArrayList<String>();
////		List<String> sample3 = new ArrayList<String>();
////		List<String> sample4 = new ArrayList<String>();
////		List<String> sample5 = new ArrayList<String>();
////		sample1.add("H");
////		sample1.add("H");
////		sample1.add("H");
////		sample2.add("T");
////		sample2.add("T");
////		sample2.add("T");
////		sample3.add("H");
////		sample3.add("H");
////		sample3.add("H");
////		sample4.add("T");
////		sample4.add("T");
////		sample4.add("T");
////		sample5.add("H");
////		sample5.add("H");
////		sample5.add("H");
////		documents.add(sample1);
////		documents.add(sample2);
////		documents.add(sample3);
////		documents.add(sample4);
////		documents.add(sample5);
////		
////		
////		List<Integer> sample1c = new ArrayList<Integer>();
////		List<Integer> sample2c = new ArrayList<Integer>();
////		List<Integer> sample3c = new ArrayList<Integer>();
////		List<Integer> sample4c = new ArrayList<Integer>();
////		List<Integer> sample5c = new ArrayList<Integer>();
////		sample1c.add(3);
////		sample1c.add(0);
////		sample2c.add(0);
////		sample2c.add(3);
////		sample3c.add(3);
////		sample3c.add(0);
////		sample4c.add(0);
////		sample4c.add(3);
////		sample5c.add(3);
////		sample5c.add(0);
////		
////		documentCounts.add(sample1c);
////		documentCounts.add(sample2c);
////		documentCounts.add(sample3c);
////		documentCounts.add(sample4c);
////		documentCounts.add(sample5c);
////		
////		
////		
////		//  doesn't correspond to preset betas so flip it
////		//  this will be better when we have to randomly determine betas etc
//////		documentCounts = model.documentListToCountMatrix(documents);
////		
////		System.out.println("Lambda: " + lambda);
////		System.out.println("Beta: " + beta);
////		System.out.println("Documents: " + documents);
////		System.out.println("Documents Counts: " + documentCounts);
////		model.logLikelihood(lambda, beta, documentCounts);
////		
////		
////		model.em(lambda, beta, documentCounts);
////		
//		
//		
//		
//		
//		
//		
//		
////		EmClusterModel model = new EmClusterModel();
////		
////		
////		List<Double> lambda = new ArrayList<Double>();
////		List<List<Double>> beta = new ArrayList<List<Double>>();
////		List<List<Integer>> documentCounts = new ArrayList<List<Integer>>();
////		List<List<String>> documents = new ArrayList<List<String>>();
////		
////		lambda.add(0.25); //  cluster1
////		lambda.add(0.25);  //  cluster2
////		lambda.add(0.25);
////		lambda.add(0.25);
////		
////		
////		
//////		List<Double> beta1 = new ArrayList<Double>();
//////		beta1.add(0.3);  //  how
//////		beta1.add(0.1);  //  it
//////		beta1.add(0.3);  //  ball
//////		beta1.add(0.3);  //  try
//////		List<Double> beta2 = new ArrayList<Double>();
//////		beta2.add(0.1);
//////		beta2.add(0.1);
//////		beta2.add(0.1);
//////		beta2.add(0.7);
//////		List<Double> beta3 = new ArrayList<Double>();
//////		beta3.add(0.3);  //  how
//////		beta3.add(0.2);  //  it
//////		beta3.add(0.3);  //  ball
//////		beta3.add(0.2);  //  try
//////		List<Double> beta4 = new ArrayList<Double>();
//////		beta4.add(0.2);  //  how
//////		beta4.add(0.2);  //  it
//////		beta4.add(0.1);  //  ball
//////		beta4.add(0.5);  //  try
//////		beta.add(beta1);
//////		beta.add(beta2);
//////		beta.add(beta3);
//////		beta.add(beta4);
////		beta.add(model.generateUninformedDirichletDistributionList(4));
////		beta.add(model.generateUninformedDirichletDistributionList(4));
////		beta.add(model.generateUninformedDirichletDistributionList(4));
////		beta.add(model.generateUninformedDirichletDistributionList(4));
////		
////		
////		
////		List<String> sample1 = new ArrayList<String>();
////		List<String> sample2 = new ArrayList<String>();
////		List<String> sample3 = new ArrayList<String>();
////		List<String> sample4 = new ArrayList<String>();
////		List<String> sample5 = new ArrayList<String>();
////		List<String> sample6 = new ArrayList<String>();
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample2.add("Ball");
////		sample2.add("Ball");
////		sample2.add("It");
////		sample2.add("It");
////		sample2.add("Try");
////		sample3.add("How");
////		sample3.add("It");
////		sample3.add("How");
////		sample3.add("It");
////		sample3.add("Try");
////		sample4.add("How");
////		sample4.add("Try");
////		sample4.add("How");
////		sample4.add("Try");
////		sample4.add("Ball");
////		sample5.add("It");
////		sample5.add("Try");
////		sample5.add("It");
////		sample5.add("Try");
////		sample5.add("Try");
////		sample6.add("It");
////		sample6.add("Try");
////		sample6.add("It");
////		sample6.add("Try");
////		sample6.add("Try");
////		documents.add(sample1);
////		documents.add(sample2);
////		documents.add(sample3);
////		documents.add(sample4);
////		documents.add(sample5);
////		documents.add(sample6);
////		
////		
////		List<Integer> sample1c = new ArrayList<Integer>();
////		List<Integer> sample2c = new ArrayList<Integer>();
////		List<Integer> sample3c = new ArrayList<Integer>();
////		List<Integer> sample4c = new ArrayList<Integer>();
////		List<Integer> sample5c = new ArrayList<Integer>();
////		List<Integer> sample6c = new ArrayList<Integer>();
////		sample1c.add(0);
////		sample1c.add(0);
////		sample1c.add(0);
////		sample1c.add(5);
////		sample2c.add(0);
////		sample2c.add(2);
////		sample2c.add(2);
////		sample2c.add(1);
////		sample3c.add(2);
////		sample3c.add(2);
////		sample3c.add(0);
////		sample3c.add(1);
////		sample4c.add(2);
////		sample4c.add(0);
////		sample4c.add(1);
////		sample4c.add(2);
////		sample5c.add(0);
////		sample5c.add(2);
////		sample5c.add(0);
////		sample5c.add(3);
////		sample6c.add(0);
////		sample6c.add(2);
////		sample6c.add(0);
////		sample6c.add(3);
////		
////		documentCounts.add(sample1c);
////		documentCounts.add(sample2c);
////		documentCounts.add(sample3c);
////		documentCounts.add(sample4c);
////		documentCounts.add(sample5c);
////		documentCounts.add(sample6c);
////		
////		
////		
////		//  doesn't correspond to preset betas so flip it
////		//  this will be better when we have to randomly determine betas etc
//////		documentCounts = model.documentListToCountMatrix(documents);
////		
////		System.out.println("Lambda: " + lambda);
////		System.out.println("Beta: " + beta);
////		System.out.println("Documents: " + documents);
////		System.out.println("Documents Counts: " + documentCounts);
////		model.logLikelihood(lambda, beta, documentCounts);
////		
////		
////		model.em(lambda, beta, documentCounts);
//		
//	
//		
//		
//	}
//	
//	public EmClusterModel(){
//		
//	}
//	
//	
//	
//	/**
//	 * This will be used later when I need to generate uniform parameters
//	 * or use a dirichlet prior or something like that. I'm just going to hard-code all this for now
//	 */
//	private void generateParameters(final int K, List<Double> lambda, List<List<Double>> beta, 
//			int[][] documentCounts, String[] documentWordOrder){
//		
//		List<Double> tmpLambda = generateUniformlyDistributedList(K);  //  K values of lambda for K clusters
//		for(int i = 0; i < tmpLambda.size(); i++){
//			lambda.add(tmpLambda.get(i));
//		}
//		
//		
//		
//		for(int i = 0; i < K; i++){  //  for each class
//			//  generate a uniform list of probabilities across all words in the vocabulary
////			beta.add(generateUniformlyDistributedList(documentWordOrder.size()));
//			beta.add(generateUninformedDirichletDistributionList(documentWordOrder.length));
//			
//		}
//	}
//	
//	/**
//	 * Generates a uniform list (all of their probabilities together sum to 1) of passed size
//	 * @param size
//	 * @return the uniform list
//	 */
//	private List<Double> generateUniformlyDistributedList(int size){
//		List<Double> uniformList = new ArrayList<Double>();
//		
//		double uniformProbability = 1.0 / (double)size;
//		//  make sure list will exactly sum to 1
//		double lastProbability = 1.0 - (uniformProbability * (size - 1)); 
//		for(int i = 0; i < size - 1; i++){
//			uniformList.add(uniformProbability);
//		}
//		uniformList.add(lastProbability);
//		
//		return uniformList;
//	}
//	
//	
//	/**
//	 * Should generate values taken from a unifrom dirichlet distribution
//	 * @param size
//	 * @return
//	 */
//	private List<Double> generateUninformedDirichletDistributionList(int size){
//		
//		List<Double> alphas = new ArrayList<Double>();
//		for(int i = 0; i < size; i++){
//			alphas.add(0.1);
//		}
//		Dirichlet d = new Dirichlet(new DenseVector(alphas));
//		List<Vec> vectorSamples = d.sample(1, new Random(System.nanoTime()));
//		
//		//  put vector double sample values into a usable list of doubles
//		List<Double> samples = new ArrayList<Double>();
////		double cumulativeProbability = 0.0;
//		for(int i = 0; i < size; i++){
//			samples.add(vectorSamples.get(0).get(i));
////			cumulativeProbability += vectorSamples.get(0).get(i);
//		}
////		samples.add(1 - cumulativeProbability);
////		System.out.println(cumulativeProbability);
//		
//		
//		return samples;
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
//	
//	
//	
//	
//	
//	
//	
//	//  TODO:  need to do all of this in logspace
//	//  TODO:   add 1 smoothing 
//	//  TODO:   don't output partial counts after I've debugged
//	//  DONE:   random dirichlet distribution
//	//  DONE:  this em can handle multiple cluster classification, it works perfectly with 4 classes
//	
//	
//	
//	
//	
//	
//	/**
//	 * Parameters these must all be array lists to ensure speed, linked lists not allowed.
//	 * (Indexes in Beta correspond to word category indexes found in documentCounts)
//	 * @param lambda  
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private void em(List<Double> lambda, List<List<Double>> beta, int[][] documentCounts, List<Integer> documentLabels, int k){
//		if(lambda.size() != beta.size() || beta.size() <= 0 || documentCounts.length <= 0 || beta.get(0).size() != documentCounts[0].length){
//			//  do not output just return becaue the parameters are not correct
//			System.out.println("EM PARAMETERS ARE INCORRECT!");
//			return;
//		}
//		
////		List<List<Double>> documentCounts = convertToLogSpace(lambda, beta, documentCountsInt);
//
//		double oldLogLikelihood = Double.NEGATIVE_INFINITY, newLogLikelihood = Double.NEGATIVE_INFINITY;
//		List<List<Double>> finalPartialCounts = new ArrayList<List<Double>>();
////		for(int i = 0; i < hardcodedIterationNumber; i++){
//		for(int i = 0;; i++){
//
//			System.out.println("\n\nIteration #" + i + ":");
//			System.out.println("Lambda: " + lambda);
//			System.out.println("Beta: " + beta);
//			
//			
//			//************************************************************************************
//			//***********************log conversion***********************************************
//			//************************************************************************************
//			
//			
//			newLogLikelihood = logLikelihood(lambda, beta, documentCounts);
//			
//			
//			
//			finalPartialCounts = emIterationLog(lambda, beta, documentCounts);
//			
//			//  if percent changed on log likelihood is less then 3.5%
//			double percentChange  = Math.abs(((Math.abs(newLogLikelihood) - Math.abs(oldLogLikelihood)) / Math.abs(oldLogLikelihood)));
//			if(oldLogLikelihood != Double.NEGATIVE_INFINITY && percentChange < 0.00035){
//				break;
//			}
//			oldLogLikelihood = newLogLikelihood;
//			
//		}
//		
//		System.out.println("\nCluster Assignments:");
//		clusterAssignment(finalPartialCounts, documentLabels, k);
//		
//		
//		
//	}
//	
//	
//	
//	/**
//	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
//	 * It goes through each class computing the partial counts for each word along the way and computers the
//	 * new lambda's and beta array for each class as it finished it's pass through every document
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private List<List<Double>> emIteration(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
//		List<Double> newLambda = new ArrayList<Double>(); 
//		List<List<Double>> newBeta = new ArrayList<List<Double>>();
//		List<List<Double>> partialCountsByClass = new ArrayList<List<Double>>();
//		
//		List<Double> denominators = new ArrayList<Double>();
//		//  for every class Cn
//		for(int i = 0; i < lambda.size(); i++){
//			
//			List<Double> partialCountsForClass = new ArrayList<Double>();
//			//  for every document DOCn
//			for(int j = 0; j < documentCounts.size(); j++){
//				
//				//  P(C=ci)
//				double partialCountTotal = lambda.get(i);
//				//  for every word in the vocabulary
//				for(int k = 0; k < beta.get(i).size(); k++){
//					//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
//					double aa = beta.get(i).get(k);
//					double bb = documentCounts.get(j).get(k);
//					partialCountTotal *= Math.pow(aa, bb); 
//				}
//				
//				//  this will happen once for each document, denominators are the same only
//				//  across the same class
//				if(denominators.size() <= j){
//					double tempDenominator = 0.0;
//					//  do above for every class given the document's counts and sum together for denominator
//					for(int c = 0; c < lambda.size(); c++){
//						double tempPartialCountForClass = lambda.get(c);
//						//  for every word in the vocabulary
//						for(int k = 0; k < beta.get(c).size(); k++){
//							//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
//							tempPartialCountForClass *= Math.pow(beta.get(c).get(k), documentCounts.get(j).get(k)); 
//						}
//						tempDenominator += tempPartialCountForClass;
//					}
//					denominators.add(tempDenominator);
//				}
//				
//				partialCountsForClass.add(partialCountTotal / denominators.get(j));
//			}
//			//  save each partial count to assign clusters after
//			partialCountsByClass.add(partialCountsForClass);
//			
//			System.out.println("PartialCounts for class " + i + ": " + partialCountsForClass);
//			
//			
//			double totalCountForAllWordsForClass = 0.0;
//			//  sum counts of every word given the class Ci
//			for(int k = 0; k < beta.get(i).size(); k++){
//				for(int d = 0; d < documentCounts.size(); d++){
//					//  this totals up all words counts for class Ci
//					totalCountForAllWordsForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
//				}
//			}
//			
//			List<Double> betaValuesForClass = new ArrayList<Double>();
//			double totalCountForOneWordForClass;
//			for(int k = 0; k < beta.get(i).size(); k++){
//				totalCountForOneWordForClass = 0.0;
//				for(int d = 0; d < documentCounts.size(); d++){
//					//  this totals accross all documents just for this word
//					totalCountForOneWordForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
//				}
//				//  nextValueAddedTo Betas for this class
//				betaValuesForClass.add(totalCountForOneWordForClass / totalCountForAllWordsForClass);
//			}
//			newBeta.add(betaValuesForClass);
//			
//			//  Lambda value for class is simple all the partial counts for the class added together 
//			double lambaValueForClass = 0.0;
//			for(int l = 0; l < partialCountsForClass.size(); l++){
//				lambaValueForClass += partialCountsForClass.get(l);
//			}
//			newLambda.add(lambaValueForClass / (double)documentCounts.size());
//		
//			
//			
//		}
//		
//		
//		//  set the original lambda and beta to the new values
//		for(int i = 0; i < lambda.size(); i++){
//			lambda.set(i, newLambda.get(i));
//		}
//		for(int i = 0; i < beta.size(); i++){
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, newBeta.get(i).get(j));
//			}
//		}
//		
//		
//		return partialCountsByClass;
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
//	
//	
//	
//	
//	/**
//	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
//	 * It goes through each class computing the partial counts for each word along the way and computers the
//	 * new lambda's and beta array for each class as it finished it's pass through every document
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private List<List<Double>> emIterationLog(List<Double> lambda, List<List<Double>> beta, int[][] documentCounts){
//		List<Double> newLambda = new ArrayList<Double>(); 
//		List<List<Double>> newBeta = new ArrayList<List<Double>>();
//		List<List<Double>> partialCountsByClass = new ArrayList<List<Double>>();
//		
//		List<Double> denominators = new ArrayList<Double>();
//		//  for every class Cn
//		for(int i = 0; i < lambda.size(); i++){
//			
//			List<Double> partialCountsForClass = new ArrayList<Double>();
//			//  for every document DOCn
//			for(int j = 0; j < documentCounts.length; j++){
//				
//				//  log(P(C=ci))
//				double partialCountTotal = Math.log(lambda.get(i));
//				//  for every word in the vocabulary
//				for(int k = 0; k < beta.get(i).size(); k++){
//					//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
//					double aa = beta.get(i).get(k);
//					double bb = documentCounts[j][k];
//					partialCountTotal += (bb * Math.log(aa)); 
//				}
//				
//				//  this will happen once for each document, denominators are the same only
//				//  across the same class
//				if(denominators.size() <= j){
//					List<Double> tempDenominator = new ArrayList<Double>();
//					//  do above for every class given the document's counts and sum together for denominator
//					for(int c = 0; c < lambda.size(); c++){
//						double tempPartialCountForClass = Math.log(lambda.get(c));
//						//  for every word in the vocabulary
//						for(int k = 0; k < beta.get(c).size(); k++){
//							//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
//							tempPartialCountForClass += (documentCounts[j][k] * Math.log(beta.get(c).get(k))); 
//						}
//						tempDenominator.add(tempPartialCountForClass);
//					}
//					denominators.add(logSum(tempDenominator));
//				}
//				
//				partialCountsForClass.add(partialCountTotal - denominators.get(j));
//			}
//			//  save each partial count to assign clusters after
//			partialCountsByClass.add(partialCountsForClass);
//			
//			System.out.println("PartialCounts for class " + i + ": " + partialCountsForClass);
//			
//			
//			//  covertPartialCountsBackIntoNormalSpace by exponentiating each value separately
//			for(int p = 0; p < partialCountsForClass.size(); p++){
//				partialCountsForClass.set(p, Math.exp(partialCountsForClass.get(p)));
//			}
//			
//			double totalCountForAllWordsForClass = 0.0;
//			//  sum counts of every word given the class Ci
//			for(int k = 0; k < beta.get(i).size(); k++){
//				for(int d = 0; d < documentCounts.length; d++){
//					//  this totals up all words counts for class Ci
//					totalCountForAllWordsForClass += (partialCountsForClass.get(d) * documentCounts[d][k]);
//				}
//			}
//			
//			List<Double> betaValuesForClass = new ArrayList<Double>();
//			double totalCountForOneWordForClass;
//			for(int k = 0; k < beta.get(i).size(); k++){
//				totalCountForOneWordForClass = 0.0;
//				for(int d = 0; d < documentCounts.length; d++){
//					//  this totals accross all documents just for this word
//					totalCountForOneWordForClass += (partialCountsForClass.get(d) * documentCounts[d][k]);
//				}
//				//  nextValueAddedTo Betas for this class
//				//  add one smoothing for betas
//				betaValuesForClass.add(((double)totalCountForOneWordForClass + 1.0) / ((double)totalCountForAllWordsForClass + (double)beta.get(0).size()));
//			}
//			newBeta.add(betaValuesForClass);
//			
//			//  Lambda value for class is simple all the partial counts for the class added together 
//			double lambaValueForClass = 0.0;
//			for(int l = 0; l < partialCountsForClass.size(); l++){
//				lambaValueForClass += partialCountsForClass.get(l);
//			}
//			//  should not smooth lambda, it yeilds incorrect values
////			newLambda.add(((double)lambaValueForClass + 1.0) / ((double)documentCounts.size() + (double)beta.get(0).size()));
//			newLambda.add(((double)lambaValueForClass) / ((double)documentCounts.length));
//
//			
//			
//		}
//		
//		
//		//  set the original lambda and beta to the new values
//		for(int i = 0; i < lambda.size(); i++){
//			lambda.set(i, newLambda.get(i));
//		}
//		for(int i = 0; i < beta.size(); i++){
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, newBeta.get(i).get(j));
//			}
//		}
//		
//		
//		return partialCountsByClass;
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
//	
//	
//	/**
//	 * Parameters these must all be array lists to ensure speed, linked lists not allowed.
//	 * (Indexes in Beta correspond to word category indexes found in documentCounts)
//	 * @param lambda  
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private void cem(List<Double> lambda, List<List<Double>> beta, int[][] documentCounts){
//		if(lambda.size() != beta.size() || beta.size() <= 0 || documentCounts.length <= 0 || beta.get(0).size() != documentCounts[0].length){
//			//  do not output just return becaue the parameters are not correct
//			System.out.println("EM PARAMETERS ARE INCORRECT!");
//			return;
//		}
//		
//
//		double oldLogLikelihood = Double.NEGATIVE_INFINITY, newLogLikelihood = Double.NEGATIVE_INFINITY;
////		List<List<Double>> finalPartialCounts = new ArrayList<List<Double>>();
//		List<Integer> highestClassIndexForDocument = new ArrayList<Integer>();
////		for(int i = 0; i < hardcodedIterationNumber; i++){
//		for(int i = 0;; i++){
//
//			System.out.println("\n\nIteration #" + i + ":");
//			System.out.println("Lambda: " + lambda);
//			System.out.println("Beta: " + beta);
//			
//			
//			newLogLikelihood = logLikelihood(lambda, beta, documentCounts);
//			
//			highestClassIndexForDocument = cemIterationLog(lambda, beta, documentCounts);
//			
//			//  if percent changed on log likelihood is less then 3.5%
//			double percentChange  = Math.abs(((Math.abs(newLogLikelihood) - Math.abs(oldLogLikelihood)) / Math.abs(oldLogLikelihood)));
//			if(oldLogLikelihood != Double.NEGATIVE_INFINITY && percentChange < 0.00035){
//				break;
//			}
//			oldLogLikelihood = newLogLikelihood;
//			
//		}
//		
//		System.out.println("\nCluster Assignments:");
//		clusterAssignmentCem(highestClassIndexForDocument);
//		
//		
//		
//	}
//	
//	
//	
//	
//	/**
//	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
//	 * It goes through each class computing the partial counts for each word along the way and computers the
//	 * new lambda's and beta array for each class as it finished it's pass through every document
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private List<Integer> cemIterationLog(List<Double> lambda, List<List<Double>> beta, int[][] documentCounts){
//		List<Double> newLambda = new ArrayList<Double>(); 
//		List<List<Double>> newBeta = new ArrayList<List<Double>>();
//		
//		//  This saves the index and value of the highest class partial count for each document
//		List<Integer> highestClassIndexForDocument = new ArrayList<Integer>();
//		List<Double> highestClassPartialCountForDocument = new ArrayList<Double>();
//		for(int i = 0; i < documentCounts.length; i++){
//			highestClassIndexForDocument.add(-1);
//			highestClassPartialCountForDocument.add(Double.NEGATIVE_INFINITY);
//		}
//		
//		//  for every class Cn
//		for(int i = 0; i < lambda.size(); i++){
//			
//			//  for every document DOCn
//			for(int j = 0; j < documentCounts.length; j++){
//				
//				//  log(P(C=ci))
//				double partialCountTotal = Math.log(lambda.get(i));
//				//  for every word in the vocabulary
//				for(int k = 0; k < beta.get(i).size(); k++){
//					//  log(P(W=wk | C=ci) ^ number of times the word appears in the Document n) => count * log(P(W|C))
//					double aa = beta.get(i).get(k);
//					double bb = documentCounts[j][k];
//					partialCountTotal += (bb * Math.log(aa)); 
//				}
//				
//				if(partialCountTotal > highestClassPartialCountForDocument.get(j)){
//					highestClassPartialCountForDocument.set(j, partialCountTotal);
//					highestClassIndexForDocument.set(j, i);
//				}
//				
//			}
//			//  save each partial count to assign clusters after, TODO: might need to fix this
//			
//		}
//		
//		
//		
//		//  set all lambda according to list generated from cem
//		for(int i = 0; i < lambda.size(); i++){
//			int classCount = 0;
//			for(int j = 0; j < highestClassIndexForDocument.size(); j++){
//				if(highestClassIndexForDocument.get(j) == i){
//					classCount++;
//				}
//			}
//			newLambda.add((double)classCount / (double)highestClassIndexForDocument.size());
//		}
//		
//		
//		
//	
//		
//		//  go through every class
//		for(int c = 0; c < lambda.size(); c++){
//			
//			
//			
//			double totalCountForAllWordsForClass = 0.0;
//			//  sum counts of every word given the class Ci
//			for(int k = 0; k < beta.get(0).size(); k++){
//				for(int d = 0; d < documentCounts.length; d++){
//					//  this totals up all words counts for class Ci
//					if(highestClassIndexForDocument.get(d) == c){
//						totalCountForAllWordsForClass += (documentCounts[d][k]);
//					}
//				}
//			}
//			
//			
//			
//			
//			
//			
//			List<Double> betaValuesForClass = new ArrayList<Double>();
//			double totalCountForOneWordForClass;
//			for(int k = 0; k < beta.get(0).size(); k++){
//				totalCountForOneWordForClass = 0.0;
//				for(int d = 0; d < documentCounts.length; d++){
//					//  this totals across all documents just for this word
//					
//					
//					if(highestClassIndexForDocument.get(d) == c){
//						totalCountForOneWordForClass += (documentCounts[d][k]);
//					}
//				}
//				//  nextValueAddedTo Betas for this class
//				//  add one smoothing for betas
//				betaValuesForClass.add(((double)totalCountForOneWordForClass + 1.0) / ((double)totalCountForAllWordsForClass + (double)beta.get(0).size()));
//			}
//			newBeta.add(betaValuesForClass);
//		}
//
//		
//		
//		
//		
//		
//		
//		
//		
//		//  set the original lambda and beta to the new values
//		for(int i = 0; i < lambda.size(); i++){
//			lambda.set(i, newLambda.get(i));
//		}
//		for(int i = 0; i < beta.size(); i++){
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, newBeta.get(i).get(j));
//			}
//		}
//		
//		
//		return highestClassIndexForDocument;
//	}
//	
//	
//	
//	
//	
//	
////	/**
////	 * Converts given lambda, beta, and documentCounts into logspace
////	 * @param lambda
////	 * @param beta
////	 * @param documentCounts
////	 * @return document counts using doubles instead of ints since it's it logspace
////	 */
////	private List<List<Double>> convertToLogSpace(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCountsInt){
////		for(int i = 0; i < beta.size(); i++){
////			lambda.set(i, Math.log(lambda.get(i)));
////			
////			for(int j = 0; j < beta.get(i).size(); j++){
////				beta.get(i).set(j, Math.log(beta.get(i).get(j)));
////			}
////		}
////		
////		List<List<Double>> documentCountsDouble = new ArrayList<List<Double>>();
////		for(int i = 0; i < documentCountsInt.size(); i++){
////			
////			List<Double> doubleCounts = new ArrayList<Double>();
////			for(int j = 0; j < documentCountsInt.get(i).size(); j++){
////				doubleCounts.add(Math.log(documentCountsInt.get(i).get(j)));
////				
////			}
////			documentCountsDouble.add(doubleCounts);
////		}
////		return documentCountsDouble;
////	}
//	
//	
//	
//	
//	/**
//	 * Adds to logs already in log space and returns an answer still in log space
//	 * @param logX
//	 * @param logY
//	 * @return
//	 */
//	 private double logAdd(double logX, double logY){
//	     // 1. make X the max
//	     if (logY > logX){
//	         double temp = logX;
//	         logX = logY;
//	         logY = temp;
//	     }
//	     // 2. now X is bigger
//	     if (logX == Double.NEGATIVE_INFINITY){
//	         return logX;
//	     }
//	     // 3. how far "down" (think decibels) is logY from logX?
//	     //    if it's really small (20 orders of magnitude smaller), then ignore
//	     double negDiff = logY - logX;
//	     if (negDiff < -20){
//	         return logX;
//	     }
//	     // 4. otherwise use some nice algebra to stay in the log domain
//	     //    (except for negDiff)
//	     return logX + Math.log(1.0 + Math.exp(negDiff));
//	 }
//	
//	
//	/**
//	 * Takes a list of logarithmic quantities
//	 * @param L
//	 * @return
//	 */
//	 private double logSum(List<Double> L){
//		double logResult = Double.NEGATIVE_INFINITY;
//		for(double logX : L){
//			logResult = logAdd(logResult, logX);
//		}
//		return logResult;
//	}
//	
//	
//	
//	
//	/**
//	 * Calculates the logLikelihood of the current iteration
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private double logLikelihood(List<Double> lambda, List<List<Double>> beta, int[][] documentCounts){
//		double logLikelihood = 0.0;
//		List<Integer> totalWordsByDocument = new ArrayList<Integer>();
//		
//		
//		for(int i = 0; i < documentCounts.length; i++){
//			int tempSum = 0;  // sum up the total number of words for each document
//			for(int j = 0; j < documentCounts[0].length; j++){
//				tempSum += documentCounts[i][j];
//			}
//			totalWordsByDocument.add(tempSum);  //  then add this total at the document's position in totalWords..
//		}
//		
//		//  sum across all documents
//		for(int i = 0; i < documentCounts.length; i++){
//			List<Double> tempValuesByClass = new ArrayList<Double>();
//			//  find intermediary values for each class
//			for(int c = 0; c < lambda.size(); c++){
//				double classValue = 0.0;
//				classValue += Math.log(lambda.get(c));
////				classValue += Math.log(factorial(totalWordsByDocument.get(i)));
//				classValue += CombinatoricsUtils.factorialLog(totalWordsByDocument.get(i));
//				
//				double sumAcrossVocabulary = 0.0;
//				for(int j = 0; j < documentCounts[i].length; j++){
//					sumAcrossVocabulary += (documentCounts[i][j] * Math.log(beta.get(c).get(j)));
////					sumAcrossVocabulary -= Math.log(factorial(documentCounts.get(i).get(j)));
//					sumAcrossVocabulary -= CombinatoricsUtils.factorialLog(documentCounts[i][j]);
//				}
//				classValue += sumAcrossVocabulary;
//				tempValuesByClass.add(classValue);
//			}
//			
//			
//			//  this only works when there's no risk of serious underflow
////			//  this part is log[e^a+.....+e^n] once a-n class values have all been found
////			double innerExpSum = 0.0;
////			for(int c = 0; c < tempValuesByClass.size(); c++){
////				innerExpSum += Math.pow(Math.E, tempValuesByClass.get(c));
////			}
////			//  add this value to the sum across all documents
////			logLikelihood += Math.log(innerExpSum);
//			
//			List<Double> logSumTerms = new ArrayList<Double>();
//			for(int c = 0; c < tempValuesByClass.size(); c++){
//				logSumTerms.add(tempValuesByClass.get(c));
//			}
//			logLikelihood += logSum(logSumTerms);
//		}
//		
//		
//		System.out.println("Log-Likelihood:" + logLikelihood);
//		return logLikelihood;
//	}
//	 
//	 
//	
//	/**
//	 * Returns the factorial of passed number
//	 * @param n
//	 * @return
//	 */
//	private double factorial(double n){
//		double result = 1.0;
//		for (double i = 1.0; i <= n; i++) {
//		   result = result * i;
//		}
//		return result;
//	}
//	
//	/**
//	 * Takes the final partial counts and uses them to assign each document
//	 * To a specific cluster
//	 * @param partialCounts
//	 */
//	private void clusterAssignment(List<List<Double>> partialCounts, List<Integer> documentLabels, int k){
//		int classAssignment = -1;
//		List<List<Integer>> contigencyMatrix = new ArrayList<List<Integer>>(); 
//		for(int i = 0; i < k; i++){
//			List<Integer> labeledRow = new ArrayList<Integer>(); 
//			for(int j = 0; j < k; j++){
//				labeledRow.add(0);
//			}
//			contigencyMatrix.add(labeledRow);
//		}
//		double maxValue = Double.NEGATIVE_INFINITY;
//		if(partialCounts.size() <= 0)
//			return;
//		
//		//  go through each document
//		for(int i = 0; i < partialCounts.get(0).size(); i++){
//			//  go through each class for this document
//			for(int j = 0; j < partialCounts.size(); j++){
//				if(partialCounts.get(j).get(i) > maxValue){
//					classAssignment = j;
//					maxValue = partialCounts.get(j).get(i);
//				}
//			}
//			System.out.println("Document " + i + " assigned to cluster " + classAssignment);
//			maxValue = Double.NEGATIVE_INFINITY;
//			contigencyMatrix.get(documentLabels.get(i)).set(classAssignment, contigencyMatrix.get(documentLabels.get(i)).get(classAssignment) + 1);
//			classAssignment = -1;
//		}
//		
//		for(int cm = 0; cm < contigencyMatrix.size(); cm++){
//			System.out.println(contigencyMatrix.get(cm));
//		}
////		System.out.println(contigencyMatrix);
//	}
//	
//	
//	
//	/**
//	 * Outputs the cluster assignments for each document for cem
//	 * @param highestClassIndexForDocument
//	 */
//	private void clusterAssignmentCem(List<Integer> highestClassIndexForDocument){
//		for(int i = 0; i < highestClassIndexForDocument.size(); i++){
//			System.out.println("Document " + i + " assigned to " + highestClassIndexForDocument.get(i));
//		}
//	}
//}
//
//











//import java.io.*;
//import java.lang.*;
//import java.math.*;
//import java.util.ArrayList;
//import java.util.List;
//import java.util.Random;
//
//import jsat.distributions.multivariate.*;
//import jsat.linear.DenseVector;
//import jsat.linear.Vec;
//
////import java.util.*;
////import org.knowceans.util.*;
//
//
//public class EmClusterModel {
//	
//	public static void main(String[] args) throws Exception{
//		EmClusterModel model = new EmClusterModel();
//				
//		TopNFeatureSelector featureSelector = new TopNFeatureSelector();
////		featureSelector.loadDocumentsFromFolder("20_newsgroups");
//		featureSelector.loadDocumentsFromFolder("20_newsgroups_Test1");
//
//		System.out.println("Documents Parsed!");
//		
//		featureSelector.topNFeatureSelection(2);
//		System.out.println("Top features selected");
//		
//		List<Double> lambda = new ArrayList<Double>();
//		List<List<Double>> beta = new ArrayList<List<Double>>();
//		final int K = 20;  //  20 clusters because there are 20 subjects
//
//		model.generateParameters(K, lambda, beta, featureSelector.getTopNDocumentCounts(), featureSelector.getTopNWordOrder());
//		
//		System.out.println("Lambda: " + lambda);
//		
//		
//		model.logLikelihood(lambda, beta, featureSelector.getTopNDocumentCounts());
//	
//	
//		model.em(lambda, beta, featureSelector.getTopNDocumentCounts());
//	
//	
//		
//		
//		
////		EmClusterModel model = new EmClusterModel();
////		
////		
////		List<Double> lambda = new ArrayList<Double>();
////		List<List<Double>> beta = new ArrayList<List<Double>>();
////		List<List<Integer>> documentCounts = new ArrayList<List<Integer>>();
////		List<List<String>> documents = new ArrayList<List<String>>();
////		
////		lambda.add(0.25); //  cluster1
////		lambda.add(0.25);  //  cluster2
////		lambda.add(0.25);
////		lambda.add(0.25);
////		
////		
////		
////		List<Double> beta1 = new ArrayList<Double>();
////		beta1.add(0.3);  //  how
////		beta1.add(0.1);  //  it
////		beta1.add(0.3);  //  ball
////		beta1.add(0.3);  //  try
////		List<Double> beta2 = new ArrayList<Double>();
////		beta2.add(0.1);
////		beta2.add(0.1);
////		beta2.add(0.1);
////		beta2.add(0.7);
////		List<Double> beta3 = new ArrayList<Double>();
////		beta3.add(0.3);  //  how
////		beta3.add(0.2);  //  it
////		beta3.add(0.3);  //  ball
////		beta3.add(0.2);  //  try
////		List<Double> beta4 = new ArrayList<Double>();
////		beta4.add(0.2);  //  how
////		beta4.add(0.2);  //  it
////		beta4.add(0.1);  //  ball
////		beta4.add(0.5);  //  try
////		beta.add(beta1);
////		beta.add(beta2);
////		beta.add(beta3);
////		beta.add(beta4);
////
////		
////		
////		
////		List<String> sample1 = new ArrayList<String>();
////		List<String> sample2 = new ArrayList<String>();
////		List<String> sample3 = new ArrayList<String>();
////		List<String> sample4 = new ArrayList<String>();
////		List<String> sample5 = new ArrayList<String>();
////		List<String> sample6 = new ArrayList<String>();
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample1.add("Try");
////		sample2.add("Ball");
////		sample2.add("Ball");
////		sample2.add("It");
////		sample2.add("It");
////		sample2.add("Try");
////		sample3.add("How");
////		sample3.add("It");
////		sample3.add("How");
////		sample3.add("It");
////		sample3.add("Try");
////		sample4.add("How");
////		sample4.add("Try");
////		sample4.add("How");
////		sample4.add("Try");
////		sample4.add("Ball");
////		sample5.add("It");
////		sample5.add("Try");
////		sample5.add("It");
////		sample5.add("Try");
////		sample5.add("Try");
////		sample6.add("It");
////		sample6.add("Try");
////		sample6.add("It");
////		sample6.add("Try");
////		sample6.add("Try");
////		documents.add(sample1);
////		documents.add(sample2);
////		documents.add(sample3);
////		documents.add(sample4);
////		documents.add(sample5);
////		documents.add(sample6);
////		
////		
////		List<Integer> sample1c = new ArrayList<Integer>();
////		List<Integer> sample2c = new ArrayList<Integer>();
////		List<Integer> sample3c = new ArrayList<Integer>();
////		List<Integer> sample4c = new ArrayList<Integer>();
////		List<Integer> sample5c = new ArrayList<Integer>();
////		List<Integer> sample6c = new ArrayList<Integer>();
////		sample1c.add(0);
////		sample1c.add(0);
////		sample1c.add(0);
////		sample1c.add(5);
////		sample2c.add(0);
////		sample2c.add(2);
////		sample2c.add(2);
////		sample2c.add(1);
////		sample3c.add(2);
////		sample3c.add(2);
////		sample3c.add(0);
////		sample3c.add(1);
////		sample4c.add(2);
////		sample4c.add(0);
////		sample4c.add(1);
////		sample4c.add(2);
////		sample5c.add(0);
////		sample5c.add(2);
////		sample5c.add(0);
////		sample5c.add(3);
////		sample6c.add(0);
////		sample6c.add(2);
////		sample6c.add(0);
////		sample6c.add(3);
////		
////		documentCounts.add(sample1c);
////		documentCounts.add(sample2c);
////		documentCounts.add(sample3c);
////		documentCounts.add(sample4c);
////		documentCounts.add(sample5c);
////		documentCounts.add(sample6c);
////		
////		
////		
////		//  doesn't correspond to preset betas so flip it
////		//  this will be better when we have to randomly determine betas etc
//////		documentCounts = model.documentListToCountMatrix(documents);
////		
////		System.out.println("Lambda: " + lambda);
////		System.out.println("Beta: " + beta);
////		System.out.println("Documents: " + documents);
////		System.out.println("Documents Counts: " + documentCounts);
//////		model.logLikelihood(lambda, beta, documentCounts);
////		
////		
////		model.em(lambda, beta, documentCounts);
////		
//	
//		
//		
//	}
//	
//	public void EmClusterModel(){
//		
//	}
//	
//	
//	
//	/**
//	 * This will be used later when I need to generate uniform parameters
//	 * or use a dirichlet prior or something like that. I'm just going to hard-code all this for now
//	 */
//	private void generateParameters(final int K, List<Double> lambda, List<List<Double>> beta, 
//			List<List<Integer>> documentCounts, List<String> documentWordOrder){
//		
//		List<Double> tmpLambda = generateUniformlyDistributedList(K);  //  K values of lambda for K clusters
//		for(int i = 0; i < tmpLambda.size(); i++){
//			lambda.add(tmpLambda.get(i));
//		}
//		
//		
//		
//		for(int i = 0; i < K; i++){  //  for each class
//			//  generate a uniform list of probabilities across all words in the vocabulary
////			beta.add(generateUniformlyDistributedList(documentWordOrder.size()));
//			beta.add(generateUninformedDirichletDistributionList(documentWordOrder.size()));
//			
//		}
//	}
//	
//	/**
//	 * Generates a uniform list (all of their probabilities together sum to 1) of passed size
//	 * @param size
//	 * @return the uniform list
//	 */
//	private List<Double> generateUniformlyDistributedList(int size){
//		List<Double> uniformList = new ArrayList<Double>();
//		
//		double uniformProbability = 1.0 / (double)size;
//		//  make sure list will exactly sum to 1
//		double lastProbability = 1.0 - (uniformProbability * (size - 1)); 
//		for(int i = 0; i < size - 1; i++){
//			uniformList.add(uniformProbability);
//		}
//		uniformList.add(lastProbability);
//		
//		return uniformList;
//	}
//	
//	
//	/**
//	 * Should generate values taken from a unifrom dirichlet distribution
//	 * @param size
//	 * @return
//	 */
//	private List<Double> generateUninformedDirichletDistributionList(int size){
//		
//		List<Double> alphas = new ArrayList<Double>();
//		for(int i = 0; i < size; i++){
//			alphas.add(1.0);
//		}
//		Dirichlet d = new Dirichlet(new DenseVector(alphas));
//		List<Vec> vectorSamples = d.sample(1, new Random(1));
//		
//		//  put vector double sample values into a usable list of doubles
//		List<Double> samples = new ArrayList<Double>();
////		double cumulativeProbability = 0.0;
//		for(int i = 0; i < size; i++){
//			samples.add(vectorSamples.get(0).get(i));
////			cumulativeProbability += vectorSamples.get(0).get(i);
//		}
////		samples.add(1 - cumulativeProbability);
////		System.out.println(cumulativeProbability);
//		
//		
//		return samples;
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
//	
//	
//	
//	
//	
//	
//	
//	//  TODO:  need to do all of this in logspace
//	//  TODO:   add 1 smoothing 
//	//  TODO:   don't output partial counts after I've debugged
//	//  DONE:   random dirichlet distribution
//	//  DONE:  this em can handle multiple cluster classification, it works perfectly with 4 classes
//	
//	
//	
//	
//	
//	
//	/**
//	 * Parameters these must all be array lists to ensure speed, linked lists not allowed.
//	 * (Indexes in Beta correspond to word category indexes found in documentCounts)
//	 * @param lambda  
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private void em(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
//		if(lambda.size() != beta.size() || beta.size() <= 0 || documentCounts.size() <= 0 || beta.get(0).size() != documentCounts.get(0).size()){
//			//  do not output just return becaue the parameters are not correct
//			System.out.println("EM PARAMETERS ARE INCORRECT!");
//			return;
//		}
//		
////		List<List<Double>> documentCounts = convertToLogSpace(lambda, beta, documentCountsInt);
//
//		double oldLogLikelihood = Double.NEGATIVE_INFINITY, newLogLikelihood = Double.NEGATIVE_INFINITY;
//		List<List<Double>> finalPartialCounts = new ArrayList<List<Double>>();
////		for(int i = 0; i < hardcodedIterationNumber; i++){
//		for(int i = 0;; i++){
//
//			System.out.println("\n\nIteration #" + i + ":");
//			System.out.println("Lambda: " + lambda);
//			System.out.println("Beta: " + beta);
//			
//			
//			//************************************************************************************
//			//***********************log conversion***********************************************
//			//************************************************************************************
//			
//			
//			
//			newLogLikelihood = logLikelihood(lambda, beta, documentCounts);
//			
//			return;
//			
////			finalPartialCounts = emIteration(lambda, beta, documentCounts);
////			
////			//  if percent changed on log likelihood is less then 3.5%
////			double percentChange  = Math.abs(((Math.abs(newLogLikelihood) - Math.abs(oldLogLikelihood)) / Math.abs(oldLogLikelihood)));
////			if(oldLogLikelihood != Double.NEGATIVE_INFINITY && percentChange < 0.035){
////				break;
////			}
////			oldLogLikelihood = newLogLikelihood;
//			
//		}
//		
////		System.out.println("\nCluster Assignments:");
////		clusterAssignment(finalPartialCounts);
//		
//		
//		
//	}
//	
//	
//	
//	/**
//	 * Iterates through one entire round of em given, lambda, beta, and the document counts lists
//	 * It goes through each class computing the partial counts for each word along the way and computers the
//	 * new lambda's and beta array for each class as it finished it's pass through every document
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private List<List<Double>> emIteration(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
//		List<Double> newLambda = new ArrayList<Double>(); 
//		List<List<Double>> newBeta = new ArrayList<List<Double>>();
//		List<List<Double>> partialCountsByClass = new ArrayList<List<Double>>();
//		
//		List<Double> denominators = new ArrayList<Double>();
//		//  for every class Cn
//		for(int i = 0; i < lambda.size(); i++){
//			
//			List<Double> partialCountsForClass = new ArrayList<Double>();
//			//  for every document DOCn
//			for(int j = 0; j < documentCounts.size(); j++){
//				
//				//  P(C=ci)
//				double partialCountTotal = lambda.get(i);
//				//  for every word in the vocabulary
//				for(int k = 0; k < beta.get(i).size(); k++){
//					//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
//					double aa = beta.get(i).get(k);
//					double bb = documentCounts.get(j).get(k);
//					partialCountTotal *= Math.pow(aa, bb); 
//				}
//				
//				//  this will happen once for each document, denominators are the same only
//				//  across the same class
//				if(denominators.size() <= j){
//					double tempDenominator = 0.0;
//					//  do above for every class given the document's counts and sum together for denominator
//					for(int c = 0; c < lambda.size(); c++){
//						double tempPartialCountForClass = lambda.get(c);
//						//  for every word in the vocabulary
//						for(int k = 0; k < beta.get(c).size(); k++){
//							//  P(W=wk | C=ci) ^ number of times the word appears in the Document n
//							tempPartialCountForClass *= Math.pow(beta.get(c).get(k), documentCounts.get(j).get(k)); 
//						}
//						tempDenominator += tempPartialCountForClass;
//					}
//					denominators.add(tempDenominator);
//				}
//				
//				partialCountsForClass.add(partialCountTotal / denominators.get(j));
//			}
//			//  save each partial count to assign clusters after
//			partialCountsByClass.add(partialCountsForClass);
//			
//			System.out.println("PartialCounts for class " + i + ": " + partialCountsForClass);
//			
//			
//			double totalCountForAllWordsForClass = 0.0;
//			//  sum counts of every word given the class Ci
//			for(int k = 0; k < beta.get(i).size(); k++){
//				for(int d = 0; d < documentCounts.size(); d++){
//					//  this totals up all words counts for class Ci
//					totalCountForAllWordsForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
//				}
//			}
//			
//			List<Double> betaValuesForClass = new ArrayList<Double>();
//			double totalCountForOneWordForClass;
//			for(int k = 0; k < beta.get(i).size(); k++){
//				totalCountForOneWordForClass = 0.0;
//				for(int d = 0; d < documentCounts.size(); d++){
//					//  this totals accross all documents just for this word
//					totalCountForOneWordForClass += (partialCountsForClass.get(d) * documentCounts.get(d).get(k));
//				}
//				//  nextValueAddedTo Betas for this class
//				betaValuesForClass.add(totalCountForOneWordForClass / totalCountForAllWordsForClass);
//			}
//			newBeta.add(betaValuesForClass);
//			
//			//  Lambda value for class is simple all the partial counts for the class added together 
//			double lambaValueForClass = 0.0;
//			for(int l = 0; l < partialCountsForClass.size(); l++){
//				lambaValueForClass += partialCountsForClass.get(l);
//			}
//			newLambda.add(lambaValueForClass / (double)documentCounts.size());
//		
//			
//			
//		}
//		
//		
//		//  set the original lambda and beta to the new values
//		for(int i = 0; i < lambda.size(); i++){
//			lambda.set(i, newLambda.get(i));
//		}
//		for(int i = 0; i < beta.size(); i++){
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, newBeta.get(i).get(j));
//			}
//		}
//		
//		
//		return partialCountsByClass;
//	}
//	
//	
//	/**
//	 * Converts given lambda, beta, and documentCounts into logspace
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 * @return document counts using doubles instead of ints since it's it logspace
//	 */
//	private List<List<Double>> convertToLogSpace(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCountsInt){
//		for(int i = 0; i < beta.size(); i++){
//			lambda.set(i, Math.log(lambda.get(i)));
//			
//			for(int j = 0; j < beta.get(i).size(); j++){
//				beta.get(i).set(j, Math.log(beta.get(i).get(j)));
//			}
//		}
//		
//		List<List<Double>> documentCountsDouble = new ArrayList<List<Double>>();
//		for(int i = 0; i < documentCountsInt.size(); i++){
//			
//			List<Double> doubleCounts = new ArrayList<Double>();
//			for(int j = 0; j < documentCountsInt.get(i).size(); j++){
//				doubleCounts.add(Math.log(documentCountsInt.get(i).get(j)));
//				
//			}
//			documentCountsDouble.add(doubleCounts);
//		}
//		return documentCountsDouble;
//	}
//	
//	
//	
//	
//	/**
//	 * Adds to logs already in log space and returns an answer still in log space
//	 * @param logX
//	 * @param logY
//	 * @return
//	 */
//	 private double logAdd(double logX, double logY){
//	     // 1. make X the max
//	     if (logY > logX){
//	         double temp = logX;
//	         logX = logY;
//	         logY = temp;
//	     }
//	     // 2. now X is bigger
//	     if (logX == Double.NEGATIVE_INFINITY){
//	         return logX;
//	     }
//	     // 3. how far "down" (think decibels) is logY from logX?
//	     //    if it's really small (20 orders of magnitude smaller), then ignore
//	     double negDiff = logY - logX;
//	     if (negDiff < -20){
//	         return logX;
//	     }
//	     // 4. otherwise use some nice algebra to stay in the log domain
//	     //    (except for negDiff)
//	     return logX + Math.log(1.0 + Math.exp(negDiff));
//	 }
//	
//	
//	/**
//	 * Takes a list of logarithmic quantities
//	 * @param L
//	 * @return
//	 */
//	 private double logSum(List<Double> L){
//		double logResult = Double.NEGATIVE_INFINITY;
//		for(double logX : L){
//			logResult = logAdd(logResult, logX);
//		}
//		return logResult;
//	}
//	
//	
//	
//	
//	/**
//	 * Calculates the logLikelihood of the current iteration
//	 * @param lambda
//	 * @param beta
//	 * @param documentCounts
//	 */
//	private double logLikelihood(List<Double> lambda, List<List<Double>> beta, List<List<Integer>> documentCounts){
//		double logLikelihood = 0.0;
//		List<Integer> totalWordsByDocument = new ArrayList<Integer>();
//		
//		
//		//************************************************************************************
//		//***********************log conversion***********************************************
//		//************************************************************************************
//		
//		
//		
//		for(int i = 0; i < documentCounts.size(); i++){
//			int tempSum = 0;  // sum up the total number of words for each document
//			for(int j = 0; j < documentCounts.get(0).size(); j++){
//				tempSum += documentCounts.get(i).get(j);
//			}
//			totalWordsByDocument.add(tempSum);  //  then add this total at the document's position in totalWords..
//		}
//		
//		//  sum across all documents
//		for(int i = 0; i < documentCounts.size(); i++){
//			List<Double> tempValuesByClass = new ArrayList<Double>();
//			//  find intermediary values for each class
//			for(int c = 0; c < lambda.size(); c++){
//				double classValue = 0.0;
//				classValue += Math.log(lambda.get(c));
//				classValue += Math.log(factorial(totalWordsByDocument.get(i)));
//				
//				double sumAcrossVocabulary = 0.0;
//				for(int j = 0; j < documentCounts.get(i).size(); j++){
//					sumAcrossVocabulary += (documentCounts.get(i).get(j) * Math.log(beta.get(c).get(j)));
//					sumAcrossVocabulary -= Math.log(factorial(documentCounts.get(i).get(j)));
//				}
//				classValue += sumAcrossVocabulary;
//				tempValuesByClass.add(classValue);
//			}
//			
//			
//			
//			
//			//  this part is log[e^a+.....+e^n] once a-n class values have all been found
//			double innerExpSum = 0.0;
//			for(int c = 0; c < tempValuesByClass.size(); c++){
//				innerExpSum += Math.pow(Math.E, tempValuesByClass.get(c));
//			}
//			//  add this value to the sum across all documents
//			logLikelihood += Math.log(innerExpSum);
//		}
//		
//		
//		System.out.println("Log-Likelihood:" + logLikelihood);
//		return logLikelihood;
//	}
//	
//	/**
//	 * Returns the factorial of passed number
//	 * @param n
//	 * @return
//	 */
//	private int factorial(int n){
//		int result = 1;
//		for (int i = 1; i <= n; i++) {
//		   result = result * i;
//		}
//		return result;
//	}
//	
//	/**
//	 * Takes the final partial counts and uses them to assign each document
//	 * To a specific cluster
//	 * @param partialCounts
//	 */
//	private void clusterAssignment(List<List<Double>> partialCounts){
//		int classAssignment = -1;
//		double maxValue = Double.NEGATIVE_INFINITY;
//		if(partialCounts.size() <= 0)
//			return;
//		
//		//  go through each document
//		for(int i = 0; i < partialCounts.get(0).size(); i++){
//			//  go through each class for this document
//			for(int j = 0; j < partialCounts.size(); j++){
//				if(partialCounts.get(j).get(i) > maxValue){
//					classAssignment = j;
//					maxValue = partialCounts.get(j).get(i);
//				}
//			}
//			System.out.println("Document " + i + " assigned to cluster " + classAssignment);
//			classAssignment = -1;
//			maxValue = Double.NEGATIVE_INFINITY;
//		}
//	}
//}

