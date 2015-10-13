import java.util.ArrayList;


public class BackwardAlgo {

	String test;
	
	ArrayList <Integer> emission; //the emission converted to generic index
	ArrayList <NodeState> stateArray;
	ArrayList <Double> pi;
	ArrayList <Double> piLog;
	ArrayList <String> emissionSet; //set of observation symbols
	HMM currHMM;
	ArrayList <String> sequences;//list the sequences to be analyzed
	ArrayList <String> descriptions; //list the description of the sequences to be analyzed
	
	public BackwardAlgo(HMM hmm, String fileName){
		
		
		emission = new ArrayList<Integer>();
		stateArray = new ArrayList<NodeState>();
		pi = new ArrayList<Double>();
		//emissionArray = new ArrayList<String>();
		emissionSet = new ArrayList<String>();
		currHMM = hmm;
		stateArray = currHMM.getStateArray();
		pi = currHMM.getPiArray();
		piLog = currHMM.getPiLogArray();
		emissionSet = currHMM.getEmissionSetArray();
	    

		//get sequence here
		FASTAReader fr = new FASTAReader(fileName);
		sequences = fr.getSequences();
		descriptions = fr.getDescriptions();
		
		for(int seqIndex = 0; seqIndex <sequences.size(); seqIndex++){//iterates through the sequences to be analyzed
			test = sequences.get(seqIndex);
			ObservationConverter oc = new ObservationConverter();
			emission = oc.processObservation(test, emissionSet);
			//calculateBeta(seqIndex); // computes beta for each state for each time point
			//calculateLogBeta(seqIndex);
			calculateLogBetaTest(seqIndex);
		}
	}//end of BackwardAlgo(HMM hmm, String fileName)
	
	public void calculateLogBetaTest(int seqIndex){//given the observation

		//Big O(TN^2)
		NodeState currNode;//i = state at t
		NodeState nextNode;//j = state at t+1
		int nextE;//emission at t+1
		double nextLogEP;//log emission probability at t+1
		double logBeta;//being calculated for all N states
		double logTransition;
		double nextLogBeta;//previous beta calculation
		double firstLogBeta = 0;//what is being subtracted in the exponent
		double firstLogTransition = 0;
		double firstLogEP = 0;
		double exponent;
		double e = 0;
		double eSum;
		int firstJ;
		double first = 0;
		double canFirst;
		
		int emissionT = emission.size()-1;// T = 0 to observation length-1
		for(int t = emissionT; t >= 0; t--){	
			nextE = emission.get(t);
			for(int i = 0; i < stateArray.size(); i++){//populates values at T
				currNode = stateArray.get(i);
				if(t == emissionT){
					logBeta = Math.log(1);
				}else{
					//Start calculate max
					firstJ = 0;
					for(int jm = 0; jm < stateArray.size(); jm++){
						nextNode = stateArray.get(jm);
						nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
						nextLogEP = nextNode.getLogEmission(nextE);
						logTransition = currNode.getLogTransition(jm);
						if(jm == 0){
							first = nextLogBeta+logTransition+nextLogEP;
						}else{
							canFirst = nextLogBeta+logTransition+nextLogEP;
							if(canFirst > first){
								first = canFirst;
								firstJ = jm;
							}
						}
					}
					
					NodeState firstNode = stateArray.get(firstJ);
					firstLogBeta = firstNode.getLogBeta(seqIndex, t+1);
					firstLogTransition = currNode.getLogTransition(firstJ);
					firstLogEP = firstNode.getLogEmission(nextE);
					
					//end calculate max
					eSum = 0;
					for(int j = 0; j <stateArray.size(); j++){//iterates through all possible future states
						nextNode = stateArray.get(j);
						nextLogEP = nextNode.getLogEmission(nextE);
						logTransition = currNode.getLogTransition(j);// j to k transition probability
						nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
						if(j != firstJ){
							exponent = nextLogBeta+logTransition+nextLogEP-firstLogBeta
								-firstLogTransition-firstLogEP;
							e = Math.exp(exponent);
							if(e == Double.NEGATIVE_INFINITY){
								e = 0;
								System.out.println("NEGATIVE INFINITY");
							}
							if(e == Double.POSITIVE_INFINITY){
								System.out.println("POSITIVE INFINITY FROM BETA");
							}
							if(e == Double.NaN){
								System.out.println("NOT A NUMBER");
							}
							eSum = eSum + e;
						}
					}
					logBeta = firstLogBeta+firstLogTransition+firstLogEP+Math.log(1+eSum);
					//System.out.println("A logBeta"+logBeta);
				}
				currNode.setLogBeta(seqIndex, emission.size(), t, logBeta);	
			}
		}	
	}//end of calculateLogBeta()
	
	public void calculateLogBeta(int seqIndex){//given the observation

		//Big O(TN^2)
		NodeState currNode;//i = state at t
		NodeState nextNode;//j = state at t+1
		int nextE;//emission at t+1
		double nextLogEP;//log emission probability at t+1
		double logBeta;//being calculated for all N states
		double logTransition;
		double nextLogBeta;//previous beta calculation
		double firstLogBeta = 0;//what is being subtracted in the exponent
		double firstLogTransition = 0;
		double firstLogEP = 0;
		double exponent;
		double e = 0;
		double eSum;
		int emissionT = emission.size()-1;// T = 0 to observation length-1
		for(int i = emissionT; i >= 0; i--){	
			nextE = emission.get(i);
			for(int j = 0; j < stateArray.size(); j++){//iterates through all states
				currNode = stateArray.get(j);
				if(i==emissionT){
					logBeta = Math.log(1);
				}else{
					eSum = 0;
					for(int k = 0; k <stateArray.size(); k++){//iterates through all possible future states
						nextNode = stateArray.get(k);
						nextLogEP = nextNode.getLogEmission(nextE);
						logTransition = currNode.getLogTransition(k);// j to k transition probability
						nextLogBeta = nextNode.getLogBeta(seqIndex, i+1);
						if(j == 0 && k == 0){
							firstLogBeta = nextLogBeta;
							firstLogTransition = logTransition;
							firstLogEP = nextLogEP;
						}else{//calculate e to add to eSum
							exponent = nextLogBeta+logTransition+nextLogEP-firstLogBeta-firstLogTransition-firstLogEP;
							e = Math.exp(exponent);
							eSum = eSum + e;
						}
					}
					logBeta = firstLogBeta+firstLogTransition+firstLogEP+Math.log(1+eSum);
					//System.out.println("A logBeta"+logBeta);
				}
				currNode.setLogBeta(seqIndex, emission.size(), i, logBeta);	
			}
		}	
	}//end of calculateLogBeta()
	
	public void calculateBeta(int seqIndex){//given the observation

		//Big O(TN^2)
		NodeState currNode;//i = state at t
		NodeState nextNode;//j = state at t+1
		int nextE;//emission at t+1
		double nextEP;//emission probability at t+1
		double beta;//being calculated for all N states
		double nextbeta;//previous beta calculation
		int emissionT = emission.size()-1;// T = 0 to size -1
		for(int i = emissionT; i >= 0; i--){	
			nextE = emission.get(i);
			for(int j = 0; j < stateArray.size(); j++){//iterates through all states
				currNode = stateArray.get(j);
				nextbeta = 0;
				if(i==emissionT){
					beta = 1;
				}else{
					for(int k = 0; k <stateArray.size(); k++){//iterates through all possible future states
						nextNode = stateArray.get(k);
						nextEP = nextNode.getEmission(nextE);
						double transitionProbability = currNode.getTransition(k);// j to k transition probability; changed j to k
						nextbeta = nextbeta + (nextNode.getBeta(seqIndex, i+1)*transitionProbability*nextEP);
						//System.out.println("At state i: "+j+" to state j: "+k);
						//System.out.println("The next emission is: "+nextE);
						//System.out.println("The next emission probability is: "+nextEP);	
						//System.out.println("The next transition prob is: "+transitionProbability);
						//System.out.println("The next beta is: "+nextNode.getBeta(seqIndex, i+1));
						//System.out.println("The calc beta is: "+(nextNode.getBeta(seqIndex, i+1)*transitionProbability*nextEP));
						//System.out.println("The sum beta is: "+nextbeta);
					}
				beta = nextbeta;
				}//end of if-else statement calculating alpha
				currNode.setBeta(seqIndex, emission.size(), i, beta);	
			}
		}	
	}//end of calculateBeta()
	
	public void displayScores(){
		for(int i = 0; i < sequences.size(); i++){
			//System.out.println("The BackwardAlgo score for sequence "+i+" is "+getScore(i));
			System.out.println("The BackwardAlgo logScore for sequence "+i+" is "+getLogScore(i));
		}
	}//end displayScores()
	
	public double getLogScore(int seqIndex){
		NodeState currNode;
		double currLogBeta;
		double currLogPi;
		double currLogEP;
		double firstLogBeta = 0;
		double firstLogPi = 0;
		double firstLogEP = 0;
		double exponent;
		double e = 0;
		double eSum = 0;
		int firstEmission;
		
		double logScore = 0;
		
		String currSeq = sequences.get(seqIndex);
		ObservationConverter oc = new ObservationConverter();
		ArrayList<Integer>currEmission = oc.processObservation(currSeq, emissionSet);
		firstEmission = currEmission.get(0);
		
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			currLogEP = currNode.getLogEmission(firstEmission);
			currLogPi = piLog.get(i);
			currLogBeta = currNode.getLogBeta(seqIndex, 0); 		

			if(i == 0){
				firstLogEP = currLogEP;
				firstLogPi = currLogPi;
				firstLogBeta = 	currLogBeta;
			}else{//calculate eSum
				exponent = currLogBeta +currLogPi+currLogEP-firstLogBeta-firstLogPi-firstLogEP;
				e =	Math.exp(exponent);	
				eSum = eSum +e;
			}
			
		}
		
		logScore = firstLogBeta+firstLogPi+firstLogEP+Math.log(1+eSum);
		return logScore;
	}//end of getLogScore()
	
	public double getScore(int seqIndex){
		double score = 0;
		NodeState currNode;
		String currSeq = sequences.get(seqIndex);
		ObservationConverter oc = new ObservationConverter();
		ArrayList<Integer>currEmission = oc.processObservation(currSeq, emissionSet);
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			score = score + (currNode.getBeta(seqIndex, 0)*pi.get(i)*currNode
					.getEmission(currEmission.get(0)));
		}
		return score;
	}//end of getScore()	
}// end of class
