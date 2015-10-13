//this class uses alpha and beta calculated by the forward and backward
//algorithms to calculate gamma and xi

import java.util.ArrayList;


public class ForwardBackward {
	
	String test;	
	ArrayList <Integer> emission; //the emission converted to generic int index
	ArrayList <NodeState> stateArray;
	ArrayList <Double> pi;
	ArrayList <String> emissionSet; //set of observation symbols
	HMM currHMM;
	ForwardAlgo fa;
	BackwardAlgo ba;
	ArrayList <String> sequences;//list the sequences to be analyzed
	ArrayList <String> descriptions; //list the description of the sequences to be analyzed
	
	public ForwardBackward(HMM hmm, String fileName){
	
		emission = new ArrayList<Integer>();
		stateArray = new ArrayList<NodeState>();
		pi = new ArrayList<Double>();
		//emissionArray = new ArrayList<String>();
		emissionSet = new ArrayList<String>();
		
		currHMM = hmm;
		stateArray = currHMM.getStateArray();
		pi = currHMM.getPiArray();
		emissionSet = currHMM.getEmissionSetArray();
	    //get sequence here
		
		FASTAReader fr = new FASTAReader(fileName);
		sequences = fr.getSequences();
		descriptions = fr.getDescriptions();
		
		for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
			test = sequences.get(seqIndex);
			ObservationConverter oc = new ObservationConverter();
			emission = oc.processObservation(test, emissionSet);
			//calculateGamma(seqIndex, fileName);
			//calculateLogGamma(seqIndex, fileName);
			calculateLogGammaTest(seqIndex, fileName);
			//calculateXi(seqIndex);
			//calculateLogXi(seqIndex);
			calculateLogXiTest(seqIndex);
		}
	}//end ForwardBackward()
	
	private void calculateLogGammaTest(int seqIndex, String fileName){
		fa = new ForwardAlgo(currHMM, fileName);
		ba = new BackwardAlgo(currHMM, fileName);
		NodeState currNode = null;
		NodeState allNode = null;
		double e;
		double eSum;
		double exponent;
		double logGamma;
		double firstLogAlpha = 0;
		double firstLogBeta = 0;
		double currLogAlpha;
		double currLogBeta;
		double allLogAlpha;
		double allLogBeta;
		int firstIndex = 0;//make sure this is the biggest for numerical stability
		double canFirst;
		double first;
		for(int t = 0; t <emission.size(); t++){
			
			//start to find firstIndex
			first = 0;
			for(int im = 0; im < stateArray.size(); im++){
				currNode = stateArray.get(im);
				if(im == 0){
					first = currNode.getLogAlpha(seqIndex, t)+currNode.getLogBeta(seqIndex, t);
				}else{
					canFirst = currNode.getLogAlpha(seqIndex, t)+currNode.getLogBeta(seqIndex, t);
					if(canFirst > first){
						first = canFirst;
						firstIndex = im;
					}
				}
			}
			NodeState firstState = stateArray.get(firstIndex);
			firstLogAlpha = firstState.getLogAlpha(seqIndex, t);
			firstLogBeta = firstState.getLogBeta(seqIndex, t);
			//System.out.println(maxIndex);
			
			for(int i = 0; i < stateArray.size(); i++){
				currNode = stateArray.get(i);
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				currLogBeta = currNode.getLogBeta(seqIndex, t);
				eSum = 0;
				for(int ai=0; ai < stateArray.size(); ai++){// SUM(j)[alpha(j)*beta(j)]
					allNode = stateArray.get(ai);
					allLogAlpha = allNode.getLogAlpha(seqIndex, t);
					allLogBeta = allNode.getLogBeta(seqIndex, t);
					if(ai != firstIndex){
						exponent = allLogAlpha+allLogBeta-firstLogAlpha-firstLogBeta;
						e = Math.exp(exponent);
						
						
						if(e == Double.NEGATIVE_INFINITY){
							e = 0;
							System.out.println("NEGATIVE INFINITY");
						}
						if(e == Double.POSITIVE_INFINITY){
							System.out.println("POSITIVE INFINITY FROM GAMMA");
						}
						if(e == Double.NaN){
							System.out.println("NOT A NUMBER");
						}
						
						
						eSum = eSum+e;
					}
				}
				logGamma = currLogAlpha+currLogBeta-(firstLogAlpha+firstLogBeta+Math.log(1+eSum));
				//System.out.println("logGamma for i = "+i+" at t = "+t+" is: "+logGamma);//QC code
				currNode.setLogGamma(seqIndex, logGamma);
			}	
		}
	}//end calculateLogGamma()

	
	private void calculateLogGamma(int seqIndex, String fileName){
		fa = new ForwardAlgo(currHMM, fileName);
		ba = new BackwardAlgo(currHMM, fileName);
		NodeState currNode = null;
		NodeState allNode = null;
		double e;
		double eSum;
		double exponent;
		double logGamma;
		double firstLogAlpha = 0;
		double firstLogBeta = 0;
		double currLogAlpha;
		double currLogBeta;
		double allLogAlpha;
		double allLogBeta;
		for(int t =0; t <emission.size(); t++){
			for(int i = 0; i < stateArray.size(); i++){
				currNode = stateArray.get(i);
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				currLogBeta = currNode.getLogBeta(seqIndex, t);
				eSum = 0;
				for(int ai=0; ai < stateArray.size(); ai++){
					allNode = stateArray.get(ai);
					allLogAlpha = allNode.getLogAlpha(seqIndex, t);
					allLogBeta = allNode.getLogBeta(seqIndex, t);
					if(ai ==0){
						firstLogAlpha = allLogAlpha;
						firstLogBeta = allLogBeta;
					}else{//calculate eSum
						
						exponent = allLogAlpha+allLogBeta-firstLogAlpha-firstLogBeta;
						e = Math.exp(exponent);
						eSum = eSum+e;
					}
				}
				logGamma = currLogAlpha+currLogBeta-(firstLogAlpha+firstLogBeta+Math.log(1+eSum));
				currNode.setLogGamma(seqIndex, logGamma);
			}	
		}
	}//end calculateLogGamma()

	private void calculateGamma(int seqIndex, String fileName){
		fa = new ForwardAlgo(currHMM, fileName);
		ba = new BackwardAlgo(currHMM, fileName);
		NodeState currNode;
		double alphaBetaSum;
		double gamma;
		for(int i=0; i<emission.size(); i++){
			alphaBetaSum = 0;
			for(int j = 0; j < stateArray.size(); j++){//calculates the sum of all alpha*beta for a given t
				currNode = stateArray.get(j);
				alphaBetaSum = alphaBetaSum + (currNode.getAlpha(seqIndex, i)*currNode.getBeta(seqIndex, i));
			}
			for(int k=0; k < stateArray.size(); k++){
				currNode = stateArray.get(k);
				gamma = (currNode.getAlpha(seqIndex, i)*currNode.getBeta(seqIndex, i))/alphaBetaSum;
				//System.out.println("gama for sequence "+seqIndex+" for state "+k+" at emission"+i+" is "+gamma);
				currNode.setGamma(seqIndex, gamma);//the probability of being at current state at time i given the sequence
			}
		}
	}//end calculateGamma()
	
	private void calculateLogXiTest(int seqIndex){
		NodeState currNode;
		NodeState nextNode;
		double e;
		double exponent;
		double eSum;
		double eSumSum;
		double firstLogAlpha = 0;
		double firstLogTransition = 0;
		double firstLogEP = 0;
		double firstLogBeta = 0;
		double nextLogEP;
		double currLogAlpha;
		double nextLogBeta;
		double currLogTransition;
		double denominator = 0;
		double logXi;
		
		int maxIndexI;
		int maxIndexJ;
		double currMax = 0;
		double currFirst;
		
		for(int t = 0; t < emission.size()-1; t++){//iterate through the emission	
			int nextE = emission.get(t+1);
			
			maxIndexI = 0;//zero at beginning of i
			maxIndexJ = 0;//zero at beginning of k
			for(int i = 0; i < stateArray.size(); i++){
				currNode = stateArray.get(i);
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				for(int j = 0; j <stateArray.size(); j++){			
					nextNode = stateArray.get(j);
					currLogTransition = currNode.getLogTransition(j);
					nextLogEP = nextNode.getLogEmission(nextE);
					nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);		
					if( i == 0 && j == 0){
						currMax = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta;
					}else{
						currFirst = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta;
						if(currFirst > currMax){
							currMax = currFirst;
							maxIndexI = i;
							maxIndexJ = j;
						}
					}		
				}
			}
			
			NodeState maxI = stateArray.get(maxIndexI); 
			NodeState maxJ =stateArray.get(maxIndexJ);		
			firstLogAlpha = maxI.getLogAlpha(seqIndex, t);
			firstLogTransition = maxI.getLogTransition(maxIndexJ);
			firstLogEP= maxJ.getLogEmission(nextE);
			firstLogBeta = maxJ.getLogBeta(seqIndex, t+1);
			//code for testing
			//double first = firstLogAlpha+firstLogTransition+firstLogEP+firstLogBeta;//for testing
			//if(seqIndex == 2) System.out.println("For t = "+t+" imax = "+maxIndexI+"  jmax = "+maxIndexJ);
			//end of code for testing
			eSumSum = 0;
			//calc denominator for xi
			for(int id = 0; id < stateArray.size(); id++){
				currNode = stateArray.get(id);
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				eSum = 0;
				for(int j = 0; j < stateArray.size(); j++){
					currLogTransition = currNode.getLogTransition(j);
					nextNode = stateArray.get(j);
					nextLogEP = nextNode.getLogEmission(nextE);
					nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
					if(id != maxIndexI || j != maxIndexJ){//wrong rules
						exponent = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta
							-firstLogAlpha-firstLogTransition-firstLogEP-firstLogBeta;
						//code for testing
						//double curr = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta;
						//if(t == 515 && seqIndex == 2) System.out.println("Curr Calcluated for t = 515 i = "+id+"  j = "+j+" is: "+curr);
						//end of code for testing
						e = Math.exp(exponent);
						if(e == Double.NEGATIVE_INFINITY){
							e = 0;
							System.out.println("NEGATIVE INFINITY");
						}
						if(e == Double.POSITIVE_INFINITY){
							System.out.println("POSITIVE INFINITY FROM XI");
						}
						if(e == Double.NaN){
							System.out.println("NOT A NUMBER");
						}
						eSum = eSum+e;
					}
				}
				eSumSum = eSumSum+eSum;
				denominator = firstLogAlpha+firstLogTransition+firstLogEP+firstLogBeta+Math.log(1+eSumSum);
				//System.out.println("xiTestDenominator: "+denominator);
			}
			//calc numerator and xi
			for(int in =0; in < stateArray.size(); in++){
				currNode = stateArray.get(in);
				ArrayList <Double> logXiArray = new ArrayList<Double>();
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				for(int j = 0; j < stateArray.size(); j++){
					currLogTransition = currNode.getLogTransition(j);
					nextNode = stateArray.get(j);
					nextLogEP = nextNode.getLogEmission(nextE);
					nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
					//double numerator = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta;
					logXi = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta-denominator;
					//System.out.println("currLogAlpha: "+currLogAlpha);
					//System.out.println("currLogTransition: "+currLogTransition);
					//System.out.println("nextLogEP: "+nextLogEP);
					//System.out.println("nextLogBeta: "+nextLogBeta);
					//System.out.println("xiTestDenominator: "+denominator);
					//System.out.println("numerator i = "+in+" j = "+j+" at t = "+t+" is: "+numerator);
					//System.out.println("denominator i = "+in+" j = "+j+" at t = "+t+" is: "+denominator);
					//System.out.println("logXi i = "+in+" j = "+j+" at t = "+t+" is: "+logXi);//testing
					logXiArray.add(logXi);
				}
				currNode.setLogXi(seqIndex, logXiArray);
			}
		}//end iteration through emission
	}//end calculateLogXi

	
	private void calculateLogXi(int seqIndex){
		NodeState currNode;
		NodeState nextNode;
		double e;
		double exponent;
		double eSum;
		double eSumSum;
		double firstLogAlpha = 0;
		double firstLogTransition = 0;
		double firstLogEP = 0;
		double firstLogBeta = 0;
		double nextLogEP;
		double currLogAlpha;
		double nextLogBeta;
		double currLogTransition;
		double denominator = 0;
		double logXi;
		for(int t = 0; t < emission.size()-1; t++){//iterate through the emission	
			int currE = emission.get(t);
			int nextE = emission.get(t+1);
			eSumSum = 0;
			//calc denominator for xi
			for(int id = 0; id < stateArray.size(); id++){
				currNode = stateArray.get(id);
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				eSum = 0;
				for(int j = 0; j < stateArray.size(); j++){
					currLogTransition = currNode.getLogTransition(j);
					nextNode = stateArray.get(j);
					nextLogEP = nextNode.getLogEmission(nextE);
					nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
					if(id == 0 && j == 0){
						firstLogAlpha = currLogAlpha;
						firstLogTransition = currLogTransition;
						firstLogEP= nextLogEP;
						firstLogBeta = nextLogBeta;
					}else{//calculate e and eSum
						
						
						exponent = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta-firstLogAlpha-firstLogTransition-firstLogEP-firstLogBeta;
						e = Math.exp(exponent);
						eSum = eSum+e;
					}
				}
				eSumSum = eSumSum+eSum;
				denominator = firstLogAlpha+firstLogTransition+firstLogEP+firstLogBeta+Math.log(1+eSumSum);
			}
			
			//calc numerator and xi
			for(int in =0; in < stateArray.size(); in++){
				currNode = stateArray.get(in);
				ArrayList <Double> logXiArray = new ArrayList<Double>();
				currLogAlpha = currNode.getLogAlpha(seqIndex, t);
				for(int j = 0; j < stateArray.size(); j++){
					currLogTransition = currNode.getLogTransition(j);
					nextNode = stateArray.get(j);
					nextLogEP = nextNode.getLogEmission(nextE);
					nextLogBeta = nextNode.getLogBeta(seqIndex, t+1);
					double numerator = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta;
					logXi = currLogAlpha+currLogTransition+nextLogEP+nextLogBeta-denominator;
					//System.out.println("numerator i = "+in+" j = "+j+" at t = "+t+" is: "+numerator);
					//System.out.println("denominator i = "+in+" j = "+j+" at t = "+t+" is: "+denominator);
					//System.out.println("logXi i = "+in+" j = "+j+" at t = "+t+" is: "+logXi);
					logXiArray.add(logXi);
				}
				currNode.setLogXi(seqIndex, logXiArray);
			}
		}//end iteration through emission
	}//end calculateLogXi
	
	private void calculateXi(int seqIndex){
		NodeState currNode;
		NodeState nextNode;
		double currAlpha;
		double nextBeta;
		int nextE;//emission at t+1
		double nextEP;//emission probability of nextE at state j
		double jk;
		double jkS;
		double jkSS;
		double jkTransition;//the transition probability from j to k
		
		
		for(int i = 0; i < emission.size()-1; i++){//iterate through the emission	
			//calculates the denominator
			jkSS = 0;
			for(int j = 0; j < stateArray.size(); j++){
				currNode = stateArray.get(j);
				currAlpha = currNode.getAlpha(seqIndex, i);
				jkS = 0;
				for(int k = 0; k < stateArray.size(); k++){//iterates through k states
					nextNode = stateArray.get(k);
					nextBeta = nextNode.getBeta(seqIndex, i+1);
					nextE = emission.get(i+1);
					nextEP = nextNode.getEmission(nextE);
					jkTransition = currNode.getTransition(k);
					jk = currAlpha*jkTransition*nextEP*nextBeta;
					//System.out.println("denominator jk"+jk);
					jkS = jkS+jk;//sum over k
				}
				jkSS = jkSS+jkS; //sums over j
				
			}
			//calculates the numerator and xi
			for(int j1 =0; j1 < stateArray.size(); j1++){//iterates through j states
				ArrayList <Double> xiArray = new ArrayList<Double>();
				double xi;
				currNode = stateArray.get(j1);
				currAlpha = currNode.getAlpha(seqIndex, i);
				for(int k = 0; k < stateArray.size(); k++){//iterates through k states
					nextNode = stateArray.get(k);
					nextBeta = nextNode.getBeta(seqIndex, i+1);
					nextE = emission.get(i+1);
					nextEP = nextNode.getEmission(nextE);
					jkTransition = currNode.getTransition(k);
					jk = currAlpha*jkTransition*nextEP*nextBeta;//gives the numerator
					xi = jk/jkSS;
					//System.out.println("jk is "+jk+" jkSS"+jkSS);
					//System.out.println("at t = "+i+" XI = "+xi+" is being added to i = "+j1+" and j = "+k );
					xiArray.add(xi);
				}//end iteration through k states
				currNode.setXi(seqIndex, xiArray); //adds Xi values to currNode
				
			}//end iteration through j states
		}//end iteration through emission
	}//end calculateXi
	
	public String getStateProbabilities(int s, int t){
		String stateProbabilities = "";
		NodeState currNode;
		for(int i = 0; i<stateArray.size(); i++){
			currNode = stateArray.get(i);
			stateProbabilities = stateProbabilities+currNode.getGamma(s, t)+", ";
		}
		return stateProbabilities;
	}
}//end of class
