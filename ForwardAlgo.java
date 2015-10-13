import java.util.ArrayList;


public class ForwardAlgo {

	String test;
	ArrayList <Integer> emission; //the emission converted to generic index
	ArrayList <NodeState> stateArray;
	ArrayList <Double> pi;
	ArrayList <Double> piLog;
	ArrayList <String> emissionSet; //set of observation symbols
	ArrayList <String> sequences;//list the sequences to be analyzed
	ArrayList <String> descriptions; //list the description of the sequences to be analyzed
	HMM currHMM;
	
	public ForwardAlgo(HMM hmm, String fileName){
		
			
		emission = new ArrayList<Integer>();
		stateArray = new ArrayList<NodeState>();
		pi = new ArrayList<Double>();
		piLog = new ArrayList<Double>();
		emissionSet = new ArrayList<String>();
		currHMM = hmm;
		stateArray = currHMM.getStateArray();
		pi = currHMM.getPiArray();
		piLog = currHMM.getPiLogArray();
		emissionSet = currHMM.getEmissionSetArray();//list of possible emissions
	    
		//get sequence here
		FASTAReader fr = new FASTAReader(fileName);
		sequences = fr.getSequences();
		descriptions = fr.getDescriptions();//will use later on
		
		
		for(int i = 0; i <sequences.size(); i++){//iterates through all the sequences to be analyze
			
			test = sequences.get(i);
			ObservationConverter oc = new ObservationConverter();
			emission = oc.processObservation(test, emissionSet);
			//calculateAlpha(i); // computes alpha for each state for each time point
			//calculateLogAlpha(i);
			//calculateLogAlphaTest(i);
			calculateLogAlphaTest2(i);
			
		}
	}//end ForwardAlgo(HMM hmm, String fileName)
	
	private void calculateLogAlphaTest2(int seqIndex){
		NodeState currNode;
		NodeState prevNode;
		int currE;
		double currLogEP;
		double logAlpha;
		double firstPrevLogAlpha = 0;
		double prevLogAlpha;
		double firstLogTransition = 0;
		double logTransition;
		double eSum;
		int firstJ = 0;
		double first = 0;
		double canFirst;
		for(int t = 0; t < emission.size(); t++){	
			currE = emission.get(t);
			if(t == 0){
				for(int i0 = 0; i0 < stateArray.size(); i0++){
					currNode = stateArray.get(i0);
					currLogEP = currNode.getLogEmission(currE);
					logAlpha = currLogEP+piLog.get(i0);
					//System.out.println("LogAlpha at t = "+t+" for i = "+i0+" is: "+logAlpha);
					currNode.setLogAlpha(seqIndex, logAlpha);
				}	
			}else{
				for(int i = 0; i < stateArray.size(); i++){
					currNode = stateArray.get(i);
					currLogEP = currNode.getLogEmission(currE);
					eSum=0;
					//start get max
					firstJ = 0;
					for(int jm = 0; jm < stateArray.size(); jm++){
						prevNode = stateArray.get(jm);
						prevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
						logTransition = prevNode.getLogTransition(i);
						if(jm == 0){
							first = prevLogAlpha + logTransition;
						}else{
							canFirst = prevLogAlpha + logTransition;
							if(canFirst > first){
								first = canFirst;
								firstJ = jm;
							}
						}
					}
					//if(seqIndex == 2 && t < 5){System.out.println("firstJ: "+firstJ);}
					NodeState firstNode = stateArray.get(firstJ);
					firstPrevLogAlpha = firstNode.getLogAlpha(seqIndex, t-1);
					firstLogTransition = firstNode.getLogTransition(i);
					//end get max
					
					for(int j = 0; j <stateArray.size(); j++){
						prevNode = stateArray.get(j);
						prevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
						logTransition = prevNode.getLogTransition(i);
						if(j != firstJ){
							double exponent = prevLogAlpha+logTransition-firstPrevLogAlpha-firstLogTransition;
							//double candFirst = prevLogAlpha+logTransition;//QC code
							//double currFirst = firstPrevLogAlpha+firstLogTransition;//QC code
							/*
							if(seqIndex == 2 && t < 5){
								System.out.println("prevLogAlpha  at t = "+t+" for i = "+i+" and j = "+j+" is: "+prevLogAlpha);
								System.out.println("logTransition at t = "+t+" for i = "+i+" and j = "+j+" is: "+logTransition);
								System.out.println("candFirst     at t = "+t+" for i = "+i+" and j = "+j+" is: "+candFirst);
								System.out.println("currFirst     at t = "+t+" for i = "+i+" and j = "+j+" is: "+currFirst);
							}*/
							double e = Math.exp(exponent);
							if(e == Double.NEGATIVE_INFINITY){
								e = 0;
								System.out.println("NEGATIVE INFINITY");
							}
							if(e == Double.POSITIVE_INFINITY){
								System.out.println("POSITIVE INFINITY FROM ALPHA");
							}
							if(e == Double.NaN){
								System.out.println("NOT A NUMBER");
							}
							eSum = eSum+e;
						}	
					}
					logAlpha = firstPrevLogAlpha+firstLogTransition+Math.log(1+eSum)+currLogEP;
					//QC code
					//if(seqIndex == 2 && t < 10){
						//System.out.println("LogAlpha at t = "+t+" for i = "+i+" is: "+logAlpha);
					//}
					//end QC code	
					currNode.setLogAlpha(seqIndex, logAlpha);
				}
			}//end of else (when t is not 0)
		}	
	}//end calculateLogAlpha(int seqIndex)
	
	private void calculateLogAlphaTest(int seqIndex){
		NodeState currNode;
		NodeState prevNode;
		int currE;
		double currLogEP;
		double logAlpha;
		double firstPrevLogAlpha = 0;
		double prevLogAlpha;
		double firstLogTransition = 0;
		double logTransition;
		double eSum;
		int firstI = 0;
		int firstJ = 0;
		double first = 0;
		double canFirst;
		for(int t = 0; t < emission.size(); t++){	
			currE = emission.get(t);
			if(t == 0){
				for(int i0 = 0; i0 < stateArray.size(); i0++){
					currNode = stateArray.get(i0);
					currLogEP = currNode.getLogEmission(currE);
					logAlpha = currLogEP+piLog.get(i0);
					//System.out.println("LogAlpha at t = "+t+" for i = "+i0+" is: "+logAlpha);
					currNode.setLogAlpha(seqIndex, logAlpha);
				}	
			}else{
				//get max of previous alpha calculation and set it for first
				for(int im = 0; im < stateArray.size(); im++){	
					firstI = 0;
					firstJ = 0;
					for(int jm = 0; jm < stateArray.size(); jm++){
						prevNode = stateArray.get(jm);
						prevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
						logTransition = prevNode.getLogTransition(im);
						if(im == 0 && jm == 0){
							first = prevLogAlpha + logTransition;
							firstI = im;
							firstJ = jm;
						}else{
							canFirst = prevLogAlpha + logTransition;
							if(canFirst > first){
								first = canFirst;
								firstJ = jm;
								firstI = im;
							}
						}
					}
				}
				NodeState firstNode = stateArray.get(firstJ);
				firstPrevLogAlpha = firstNode.getLogAlpha(seqIndex, t-1);
				firstLogTransition = firstNode.getLogTransition(firstI);
				//end get max of previous alpha calculation
				for(int i = 0; i < stateArray.size(); i++){
					currNode = stateArray.get(i);
					currLogEP = currNode.getLogEmission(currE);
					eSum=0;
					for(int j = 0; j <stateArray.size(); j++){
						prevNode = stateArray.get(j);
						prevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
						logTransition = prevNode.getLogTransition(i);
						if(i != firstI || j != firstJ){
							double exponent = prevLogAlpha+logTransition-firstPrevLogAlpha-firstLogTransition;
							//double candFirst = prevLogAlpha+logTransition;//QC code
							//double currFirst = firstPrevLogAlpha+firstLogTransition;//QC code
							
							//if(seqIndex == 2 && t < 5){
								//System.out.println("prevLogAlpha  at t = "+t+" for i = "+i+" and j = "+j+" is: "+prevLogAlpha);
								//System.out.println("logTransition at t = "+t+" for i = "+i+" and j = "+j+" is: "+logTransition);
								//System.out.println("candFirst     at t = "+t+" for i = "+i+" and j = "+j+" is: "+candFirst);
								//System.out.println("currFirst     at t = "+t+" for i = "+i+" and j = "+j+" is: "+currFirst);
							//}
							double e = Math.exp(exponent);
							eSum = eSum+e;
						}	
					}
					logAlpha = firstPrevLogAlpha+firstLogTransition+Math.log(1+eSum)+currLogEP;
					//QC code
					//if(seqIndex == 2 && t < 10){
						//System.out.println("LogAlpha at t = "+t+" for i = "+i+" is: "+logAlpha);
					//}
					//end QC code	
					currNode.setLogAlpha(seqIndex, logAlpha);
				}
			}//end of else (when t is not 0)
		}	
	}//end calculateLogAlpha(int seqIndex)
	
	private void calculateLogAlpha(int seqIndex){
		NodeState currNode;
		int currE;
		double currLogEP;
		double logAlpha;
		double firstPrevLogAlpha = 0;
		double prevLogAlpha;
		double firstLogTransition = 0;
		double logTransition;
		double eSum;
		for(int t = 0; t < emission.size(); t++){	
			currE = emission.get(t);
			eSum = 0;
			for(int j = 0; j < stateArray.size(); j++){
				currNode = stateArray.get(j);
				currLogEP = currNode.getLogEmission(currE);
				if(t==0){
					logAlpha = currLogEP+piLog.get(j);
					//System.out.println("LogAlpha at t = "+0+" for i = "+j+" is: "+logAlpha);
				}else{
					eSum=0;
					for(int k = 0; k <stateArray.size(); k++){
						NodeState prevNode = stateArray.get(k);
						if(j == 0 && k==0){ 
							firstLogTransition = prevNode.getLogTransition(j);
							firstPrevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
						}else{
							prevLogAlpha = prevNode.getLogAlpha(seqIndex, t-1);
							logTransition = prevNode.getLogTransition(j);
							//QC code
							//if(seqIndex == 2 && t < 5){
								//System.out.println("prevLogAlpha  at t = "+t+" for i = "+j+" and j = "+k+" is: "+prevLogAlpha);
								//System.out.println("logTransition at t = "+t+" for i = "+j+" and j = "+k+" is: "+logTransition);
							//}
							//end QC code
							double diff = prevLogAlpha+logTransition-firstPrevLogAlpha-firstLogTransition;
							double e = Math.exp(diff);
							eSum = eSum+e;	
						}
					}
					logAlpha = firstPrevLogAlpha+firstLogTransition+Math.log(1+eSum)+currLogEP;
					//QC code
					//if(seqIndex == 2 && t < 5){
					//System.out.println("LogAlpha at t = "+t+" for i = "+j+" is: "+logAlpha);
					//}
					//end QCcode
				}//end of if-else statement calculating alpha
				currNode.setLogAlpha(seqIndex, logAlpha);
			}
		}	
	}//end calculateLogAlpha(int seqIndex)
	
	
	private void calculateAlpha(int seqIndex){//seqIndex refers to the sequence being analyzed
		//Big O(TN^2)
		NodeState currNode;
		int currE;
		double currEP;
		double alpha;
		double prevAlpha;
		for(int i = 0; i < emission.size(); i++){	
			currE = emission.get(i);
			//System.out.println(currE);
			for(int j = 0; j < stateArray.size(); j++){
				currNode = stateArray.get(j);
				currEP = currNode.getEmission(currE);
				prevAlpha = 0;
				if(i==0){
					alpha = currEP*pi.get(j);
				}else{
					for(int k = 0; k <stateArray.size(); k++){//sums all the previous alpha*a
						NodeState prevNode = stateArray.get(k);
						double transitionProbability = prevNode.getTransition(j);// k->j transition probability
						prevAlpha = prevAlpha + (prevNode.getAlpha(seqIndex, i-1)*transitionProbability);
						//System.out.println("The current alpha sum is: "+prevAlpha);
						//System.out.println("At state j = "+j+" from state i = "+k);
						//System.out.println("The current emission is: "+currE);
						//System.out.println("The current emission probability is: "+currEP);
						//System.out.println("The prevalpha is: "+prevNode.getAlpha(seqIndex,i-1));
						//System.out.println("The transition prob is: "+transitionProbability);
						//System.out.println("The calculated partAlpha is: "+(prevNode.getAlpha(seqIndex, i-1)*transitionProbability));
						//System.out.println("The sum partAlpha is: "+prevAlpha);
					}
					alpha = currEP*prevAlpha;
					//System.out.println("Alpha is: "+alpha);
				}//end of if-else statement calculating alpha
				currNode.setAlpha(seqIndex, alpha);

			}
		}	
	}//end of calculateAlpha()
	
	public void displayScores(){
		//String scores = "";
		String logScores = "";
		for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
			//scores = scores+getScore(seqIndex)+" ";
			logScores = logScores+getLogScore(seqIndex)+" ";
		}
		//System.out.println(scores);
		System.out.println(logScores);
		
	}//end displayScores()
	
	public double getScore(int seqIndex){ // index is for the sequence being analyzed
		double score = 0;
		NodeState currNode;
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			String currSeq = sequences.get(seqIndex);
			int seqLength = currSeq.length();
			score = score + currNode.getAlpha(seqIndex, seqLength-1);//sum of all alpha at seqLength-1 (T)
		}
		//System.out.println("The ForwardAlgo score for sequence "+index+" is "+score);
		return score;
	}//end of getScore()
	/*
	public ArrayList<Double> getLogScore(){ // index is for the sequence being analyzed
		double logScore = 0;
		double firstLogAlpha = 0;
		double currLogAlpha;
		double currE;
		double sumE = 0;
		double exponent;
		NodeState currNode;
		int seqLength;
		for(int stateIndex = 0; stateIndex < stateArray.size(); stateIndex++){
			currNode = stateArray.get(stateIndex);		
			String currSeq = sequences.get(seqIndex);
			seqLength = currSeq.length();//takes into account sequences of various lengths
			currLogAlpha = currNode.getLogAlpha(seqIndex, seqLength-1);//seqLength-1 = T
			if(stateIndex == 0){
				firstLogAlpha = currLogAlpha;
				//System.out.println("the firstLogAlpha is "+firstLogAlpha);
			}else{
				//calculate esum
				//System.out.println("the currLogAlph is "+currLogAlpha);
				exponent = currLogAlpha - firstLogAlpha;
				currE = Math.exp(exponent);
				//System.out.println("the currE is "+currE);
				sumE = sumE+currE;
			}
		}
		//System.out.println("The sumE is "+sumE);
		//System.out.println("The firstLogAlpha is "+firstLogAlpha);
		logScore = firstLogAlpha+Math.log(1+sumE);
		return logScore;
	}//end of getLogScore()
	*/
	public double getLogScore(int seqIndex){ // index is for the sequence being analyzed
		double logScore = 0;
		double firstLogAlpha = 0;
		double currLogAlpha;
		double currE;
		double sumE = 0;
		double exponent;
		NodeState currNode;
		int seqLength;
		for(int stateIndex = 0; stateIndex < stateArray.size(); stateIndex++){
			currNode = stateArray.get(stateIndex);		
			String currSeq = sequences.get(seqIndex);
			seqLength = currSeq.length();//takes into account sequences of various lengths
			currLogAlpha = currNode.getLogAlpha(seqIndex, seqLength-1);//seqLength-1 = T
			if(stateIndex == 0){
				firstLogAlpha = currLogAlpha;
				//System.out.println("the firstLogAlpha is "+firstLogAlpha);
			}else{
				//calculate esum
				//System.out.println("the currLogAlph is "+currLogAlpha);
				exponent = currLogAlpha - firstLogAlpha;
				currE = Math.exp(exponent);
				//System.out.println("the currE is "+currE);
				sumE = sumE+currE;
			}
		}
		//System.out.println("The sumE is "+sumE);
		//System.out.println("The firstLogAlpha is "+firstLogAlpha);
		logScore = firstLogAlpha+Math.log(1+sumE);
		return logScore;
	}//end of getLogScore()
	
	public ArrayList<String> getStateSeq(){//the ViterbiAlgo
		ArrayList<String>stateSeqArray = new ArrayList<String>();
		String stateSeq = "";
		NodeState currNode;
		double maxAlpha;
		int maxState;
		double currAlpha;
		
		//System.out.println("The number of sequences being analyzed is "+sequences.size());
		for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){//iterates through sequence list
			//System.out.println("the length of the sequence being analyzed is "+emission.size());
			stateSeq = "";
			String currSeq = sequences.get(seqIndex);
			int currLength = currSeq.length();
			for(int t = 0; t < currLength; t++){
				maxAlpha = 0;
				maxState = 0;
				for(int i = 0; i < stateArray.size(); i++){
					currNode = stateArray.get(i);
					currAlpha = currNode.getLogAlpha(seqIndex, t);
					if(i == 0){
						maxAlpha = currAlpha;
					}
					if(currAlpha > maxAlpha){
						maxAlpha = currAlpha;
						maxState = i;
					}
				}//end of state iteration
				stateSeq = stateSeq+maxState;
			}//end of emission iteration
			//System.out.println("From Viterbi : "+stateSeq);
			stateSeqArray.add(stateSeq);
		}//end of sequence list iteration
		return stateSeqArray;
	}//end ArrayList<String> getStateSeq()
}//end of class
