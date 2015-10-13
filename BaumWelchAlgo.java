//this class adjust the initial transition, transition, and emission probabilities 
//for a model given training sequences
import java.util.ArrayList;

public class BaumWelchAlgo {
	String test;
	ForwardBackward fb;
	ArrayList <Integer> emission; //the emission converted to generic index
	ArrayList <NodeState> stateArray;
	ArrayList <Double> pi;
	ArrayList <String> emissionSet; //set of observation symbols
	HMM currHMM;
	ArrayList <String> sequences;//list the sequences to be analyzed
	ArrayList <String> descriptions; //list the description of the sequences to be analyzed
	
	public BaumWelchAlgo(HMM hmm, String fileName){
		currHMM = hmm;
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
		fb = new ForwardBackward(currHMM, fileName);
		
	}//end BaumWelchAlgo(HMM currHMM)
	
	public String adjustParameters(){
		String newPi = newPiLogTest();//newPi();
		String newA = newALogTest();//newA();
		String newB = newBLogTest();//newB();
		String eSet = "\"BE\":[";
		for(int i = 0; i < emissionSet.size(); i++){
			if (i == emissionSet.size()-1){
				eSet = eSet+emissionSet.get(i)+"]";
			}else if(i < emissionSet.size()){
				eSet = eSet+emissionSet.get(i)+",";}
		}
		String newHMM = "{"+newPi+","+newA+","+newB+","+eSet+"}";
		//System.out.println("The adjusted HMM is : "+newHMM);
		return newHMM;
	}//end adjustParameters()
	//uses log values to calculate new pi parameters over multiple sequences
	private String newPiLogTest(){
		NodeState currNode;
		double currLogGamma;
		double firstLogGamma = 0;
		double e;
		double eSum;
		double logGammaSum = 0;
		double exponent;
		double newPi;
		double first = 0;
		double canFirst;
		int firstIndex;
		String pi= "\"pi\":[";
		
		for(int i=0; i <stateArray.size(); i++){
			currNode = stateArray.get(i);
			//begin find first
			firstIndex = 0;
			for(int seqM = 0; seqM < sequences.size(); seqM++){
				currLogGamma = currNode.getLogGamma(seqM, 0);
				if(seqM == 0){
					first = currLogGamma;
					firstIndex = 0;
				}else{
					canFirst = currLogGamma;
					if(canFirst > first){
						first = canFirst;
						firstIndex = seqM;
					}
				}
			}
			//end find first
			firstLogGamma = currNode.getLogGamma(firstIndex, 0);
			
			eSum = 0;
			for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++ ){
				currLogGamma = currNode.getLogGamma(seqIndex, 0);
				if(seqIndex != firstIndex){ 
					exponent = currLogGamma - firstLogGamma;
					e = Math.exp(exponent);
					if(e == Double.NEGATIVE_INFINITY){
						e = 0;
						System.out.println("NEGATIVE INFINITY");
					}
					if(e == Double.POSITIVE_INFINITY){
						System.out.println("POSITIVE INFINITY NEW PI");
					}
					if(e == Double.NaN){
						System.out.println("NOT A NUMBER");
					}
					eSum = eSum+e;
				}			
			}
			logGammaSum = firstLogGamma+Math.log(1+eSum);
			//newPi = Math.exp(logGammaSum)/sequences.size();
			newPi = logGammaSum - Math.log(sequences.size());
			//the next 4 lines of code formats the Pi data for json
			if (i == stateArray.size()-1){
				pi = pi+newPi+"]";
			}else if(i < stateArray.size()){
					pi = pi+newPi+",";}
			//System.out.println("The new initial transistion probability for state "+i+" is "+newPi);
		}
		//System.out.println("log newPi: "+ pi);
		return pi;
	}//end newPiLog()
	//uses log values to calculate new pi parameters over multiple sequences
	private String newPiLog(){
		NodeState currNode;
		double currLogGamma;
		double firstLogGamma = 0;
		double e;
		double eSum;
		double logGammaSum = 0;
		double exponent;
		double newPi;
		String pi= "\"pi\":[";
		for(int i=0; i <stateArray.size(); i++){
			currNode = stateArray.get(i);
			eSum = 0;
			for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++ ){
				currLogGamma = currNode.getLogGamma(seqIndex, 0);
					
					if(seqIndex == 0){ 
						//first = true;//flag incase of nan first value
						firstLogGamma = currLogGamma;
					}else{
						exponent = currLogGamma - firstLogGamma;
						e = Math.exp(exponent);
						eSum = eSum+e;
					}			
			}
			logGammaSum = firstLogGamma+Math.log(1+eSum);
			//newPi = Math.exp(logGammaSum)/sequences.size();
			newPi = logGammaSum - Math.log(sequences.size());
			//the next 4 lines of code formats the Pi data for json
			if (i == stateArray.size()-1){
				pi = pi+newPi+"]";
			}else if(i < stateArray.size()){
					pi = pi+newPi+",";}
			//System.out.println("The new initial transistion probability for state "+i+" is "+newPi);
		}
		//System.out.println("log newPi: "+ pi);
		return pi;
	}//end newPiLog()

	//works with multiple sequences to adjust the initial transition probabilities
	private String newPi(){
		
		NodeState currNode;
		double sPI;
		double newPi;
		String pi= "\"pi\":[";
		for(int i=0; i <stateArray.size(); i++){
			currNode = stateArray.get(i);
			sPI = 0;
			for(int s = 0; s < sequences.size(); s++ ){
				sPI = sPI+currNode.getGamma(s, 0);
			}
			newPi = sPI/sequences.size();
			//the next 4 lines of code formats the Pi data for json
			if (i == stateArray.size()-1){
				pi = pi+newPi+"]";
			}else if(i < stateArray.size()){
					pi = pi+newPi+",";}
			//System.out.println("The new initial transistion probability for state "+i+" is "+newPi);
		}
		//System.out.println(pi);
		return pi;
	}//end newPi()
	
	//Uses log values to calculate newA parameters over multiple sequences
		private String newALogTest(){
			NodeState currNode;
			String currSeq;
			String allA = "\"A\":[";
			String stateA ="[";
			ArrayList <Integer> currEmission;
			double e;
			double eSum;
			double eSumSum;
			double exponent;
			double currLogGamma;
			double firstLogGamma = 0;
			double logGammaSumSum;
			double currLogXi;
			double firstLogXi = 0;
			double logXiSumSum;
			double logNewTransition;
			double newTransition = 0;
			double firstGamma = 0;
			double canFirstGamma;
			int firstGammaS = 0;
			int firstGammaT = 0;
			double firstXi = 0;
			double canFirstXi;
			int firstXiS = 0;
			int firstXiT = 0;
			ArrayList<Double> logXiArray;
			ArrayList<Double> firstLogXiArray;
			ObservationConverter oc;
			//iterates through every state
			for(int i=0; i < stateArray.size(); i++){
				currNode = stateArray.get(i);	
				
				//begin find firstLogGamma
				for(int seqM = 0; seqM < sequences.size(); seqM++){
					currSeq = sequences.get(seqM); //gets a sequence to be analyzed as a string
					oc = new ObservationConverter();
					currEmission = oc.processObservation(currSeq, emissionSet); 
					//firstGammaS = 0;
					//firstGammaT = 0;
					for(int t = 0; t < currEmission.size()-1; t++){
						currLogGamma = currNode.getLogGamma(seqM, t);
						if(seqM == 0 && t == 0){
							firstGamma = currLogGamma;
							firstGammaS = 0;
							firstGammaT = 0;
						}else{
							canFirstGamma = currLogGamma;
							if(canFirstGamma > firstGamma){
								firstGamma = canFirstGamma;
								firstGammaS = seqM;
								firstGammaT = t;
							}
						}
					}
				}
				firstLogGamma = currNode.getLogGamma(firstGammaS, firstGammaT);
				//end find firstLogGamma
				/*
				The block of code below calc the denominator(logGammaSumSum)
				*/
				eSumSum = 0;
				for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
					currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
					oc = new ObservationConverter();
					currEmission = oc.processObservation(currSeq, emissionSet); 
					eSum = 0;
					for(int t = 0; t < currEmission.size()-1; t++){
						currLogGamma = currNode.getLogGamma(seqIndex, t);
						if(seqIndex != firstGammaS || t != firstGammaT){
							exponent = currLogGamma - firstLogGamma;
							e = Math.exp(exponent);
							if(e == Double.NEGATIVE_INFINITY){
								e = 0;
								System.out.println("NEGATIVE INFINITY");
							}
							if(e == Double.POSITIVE_INFINITY){
								System.out.println("POSITIVE INFINITY NEW A DENOMINATOR");
							}
							if(e == Double.NaN){
								System.out.println("NOT A NUMBER");
							}
							eSum = eSum + e;
						}
					}
					eSumSum = eSumSum + eSum;
				}
				logGammaSumSum = firstLogGamma+Math.log(1+eSumSum);
				/*
				The block of code below calc logXiSumSum (numerator)  
				*/
				stateA ="[";//wipes previous state's data
				for(int j = 0; j < stateArray.size(); j++){// Xi(i to j)
					//begin find firstLogXi
					for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
						currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
						oc = new ObservationConverter();
						// the line below converts currSeq to an array of Integers
						currEmission = oc.processObservation(currSeq, emissionSet); 
						for(int t = 0; t < currEmission.size()-1; t++){
							logXiArray = currNode.getLogXi(seqIndex, t);
							currLogXi = logXiArray.get(j);
							if(seqIndex == 0 && t == 0){
								firstXi = currLogXi;
								firstXiS = 0;
								firstXiT = 0;
							}else{
								canFirstXi = currLogXi;
								if(canFirstXi > firstXi){
									firstXi = canFirstXi;
									firstXiS = seqIndex;
									firstXiT = t;
								}
							}
						}
					}
					firstLogXiArray = currNode.getLogXi(firstXiS, firstXiT);
					firstLogXi = firstLogXiArray.get(j);
					//end find firstLogXi
					eSumSum = 0;
					for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
						currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
						oc = new ObservationConverter();
						// the line below converts currSeq to an array of Integers
						currEmission = oc.processObservation(currSeq, emissionSet); 
						eSum = 0;
						for(int t = 0; t < currEmission.size()-1; t++){
							logXiArray = currNode.getLogXi(seqIndex, t);
							currLogXi = logXiArray.get(j);
							if(seqIndex != firstXiS || t != firstXiT){//takes out firstLogXi
								exponent = currLogXi - firstLogXi;
								//System.out.println("currLog: "+currLogXi);
								//System.out.println("firstLog: "+firstLogXi);
								//System.out.println("exponent: "+exponent);
								
								e = Math.exp(exponent);
								if(e == Double.NEGATIVE_INFINITY){
									e = 0;
									System.out.println("NEGATIVE INFINITY");
								}
								if(e == Double.POSITIVE_INFINITY){
									System.out.println("POSITIVE INFINITY NEW A NUMERATOR");
								}
								if(e == Double.NaN){
									System.out.println("NOT A NUMBER");
								}
								eSum = eSum +e;
							}
						}
						eSumSum = eSumSum + eSum;
					}
					logXiSumSum = firstLogXi + Math.log(1+eSumSum);

					logNewTransition = logXiSumSum - logGammaSumSum;
					//newTransition = Math.exp(logNewTransition);
					newTransition = logNewTransition;
					//System.out.println("logXiSumSum: "+logXiSumSum);
					//System.out.println("logGammaSumSum: "+logGammaSumSum);
					//System.out.println("logNewTransition: "+logNewTransition);
					
					//the next 4 lines of code formats the stateA data for json
					if (j == stateArray.size()-1){
						stateA = stateA+newTransition+"]";
					}else if(j < stateArray.size()){
						stateA = stateA+newTransition+",";}
				}
				//the next 4 lines of code formats the allA data for json
				if (i == stateArray.size()-1){
					allA = allA+stateA+"]";
				}else if(i < stateArray.size()){
					allA = allA+stateA+",";}
			}//end of for all state i
			//System.out.println("New Log transition param: "+allA);
			return allA;
		}

	
	//Uses log values to calculate newA parameters over multiple sequences
	private String newALog(){
		NodeState currNode;
		String currSeq;
		String allA = "\"A\":[";
		String stateA ="[";
		ArrayList <Integer> currEmission;
		double e;
		double eSum;
		double eSumSum;
		double exponent;
		double currLogGamma;
		double firstLogGamma = 0;
		double logGammaSumSum;
		double currLogXi;
		double firstLogXi = 0;
		double logXiSumSum;
		double logNewTransition;
		double newTransition = 0;
		ArrayList<Double> logXiArray;
		ObservationConverter oc;
		//iterates through every state
		for(int i=0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);	
			/*
			The block of code below calc the denominator(logGammaSumSum)
			*/
			eSumSum = 0;
			for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
				currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
				oc = new ObservationConverter();
				currEmission = oc.processObservation(currSeq, emissionSet); 
				eSum = 0;
				for(int t = 0; t < currEmission.size()-1; t++){
					currLogGamma = currNode.getLogGamma(seqIndex, t);
					if(seqIndex == 0 && t == 0){
						firstLogGamma = currLogGamma;
						//if(i == 2)System.out.println("firstLogGamma"+firstLogGamma);
					}else{
						exponent = currLogGamma - firstLogGamma;
						e = Math.exp(exponent);
						eSum = eSum + e;
					}
				}
				eSumSum = eSumSum + eSum;
			}
			logGammaSumSum = firstLogGamma+Math.log(1+eSumSum);
			
			/*
			The block of code below calc logXiSumSum (numerator)  
			*/
			stateA ="[";//wipes previous state's data
			for(int j = 0; j < stateArray.size(); j++){// Xi(i to j)
				eSumSum = 0;
				for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
					currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
					oc = new ObservationConverter();
					// the line below converts currSeq to an array of Integers
					currEmission = oc.processObservation(currSeq, emissionSet); 
					eSum = 0;
					for(int t = 0; t < currEmission.size()-1; t++){
						logXiArray = currNode.getLogXi(seqIndex, t);
						currLogXi = logXiArray.get(j);
						//if(i==2)System.out.println("The currLogXi: "+currLogXi);//for QC
						if(seqIndex == 0 && t == 0){
							firstLogXi = currLogXi;
							//if(i == 2)System.out.println("The firstLogXi: "+firstLogXi);//for QC
						}else{
							exponent = currLogXi - firstLogXi;
							e = Math.exp(exponent);
							eSum = eSum +e;
						}
					}
					eSumSum = eSumSum + eSum;
				}
				logXiSumSum = firstLogXi + Math.log(1+eSumSum);

				logNewTransition = logXiSumSum - logGammaSumSum;
				//newTransition = Math.exp(logNewTransition);
				newTransition = logNewTransition;
				//System.out.println("logXiSumSum: "+logXiSumSum);
				//System.out.println("logGammaSumSum: "+logGammaSumSum);
				//System.out.println("logNewTransition: "+logNewTransition);
				
				//the next 4 lines of code formats the stateA data for json
				if (j == stateArray.size()-1){
					stateA = stateA+newTransition+"]";
				}else if(j < stateArray.size()){
					stateA = stateA+newTransition+",";}
			}//end of for loop that iterates through all states(k) at t+1
			//the next 4 lines of code formats the allA data for json
			if (i == stateArray.size()-1){
				allA = allA+stateA+"]";
			}else if(i < stateArray.size()){
				allA = allA+stateA+",";}
		}//end of for all state i
		//System.out.println("New Log transition param: "+allA);
		return allA;
	}//end newALog()

	/*
	 works with multiple sequences to adjust the transition probabilities
	 will result in NaN when sequence is one emission because of 0/0
	"A":[[0.5,0.5],[0.5,0.5]]
	*/
	private String newA(){
		NodeState currNode;
		String currSeq;
		String allA = "\"A\":[";
		String stateA ="[";
		ArrayList <Integer> currEmission;
		double gammaSum;//sum of gammas for a state over t0 to T
		double gammaSumSeq;//sum of gammas for a state over multiple sequences
		double xiSum = 0;
		double xiSumSeq;//over multiple sequences
		double newTransition;
		ArrayList<Double> stateXi;
		ObservationConverter oc;
		//iterates through every state
		for(int i=0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);	
			/*
			The block of code below calculates gammaSumSeq for state i for all sequences
			to be analyzed, from t0 to T-1. 
			*/
			//iterates through every sequence to be analyzed	
			gammaSumSeq = 0;//sum of all gamma from state i for all sequences to be analyzed 
			for(int s = 0; s < sequences.size(); s++){
				currSeq = sequences.get(s); //gets a sequence to be analyzed as a string
				oc = new ObservationConverter();
				// the line below converts currSeq to an array of Integers
				currEmission = oc.processObservation(currSeq, emissionSet); 
				gammaSum = 0;//sum of all gamma from t0 to T for state i for sequence s 
				for(int j = 0; j < currEmission.size()-1; j++){//t0 to T-1 to get gammaSum
					gammaSum = gammaSum + currNode.getGamma(s, j);
				}//end of for loop to get gammaSum
				gammaSumSeq = gammaSumSeq + gammaSum;
			}//end of for all sequences to be analyzed (to get gammaSumSeq)
			//System.out.println("The gammaSumSeq for state = "+i+" is "+gammaSumSeq);
			//end of block
			/*
			The block of code below calculates xiSumSeq for state i for
			all sequences to be analyzed.   
			*/
			stateA ="[";//wipes previous state's data
			//the for loop below iterates through all states at t+1
			for(int k = 0; k < stateArray.size(); k++){//j
				xiSumSeq = 0;
				for(int s = 0; s < sequences.size(); s++){
					currSeq = sequences.get(s); //gets a sequence to be analyzed as a string
					oc = new ObservationConverter();
					// the line below converts currSeq to an array of Integers
					currEmission = oc.processObservation(currSeq, emissionSet); 
					xiSum = 0;
					for(int l = 0; l < currEmission.size()-1; l++){//t0 to T-1 to get xiSum
						stateXi = currNode.getXi(s, l);//holds xi values for state i for all states at t = l
						//System.out.println("For sequence s = "+s+" the xi value at t = "+l+" "
							//+"from state i ="+i+" to state j = "+k+" is "+stateXi.get(k));
						xiSum = xiSum+stateXi.get(k);
					}//end for loop from t0 to T to get xiSum
					xiSumSeq = xiSumSeq +xiSum;
				}//end for all sequences to be analyzed	
				newTransition = xiSumSeq/gammaSumSeq;
				//System.out.println("The new transition from state i = "+i+" to state j = "+k+" "
					//+ "is "+newTransition);
				//the next 4 lines of code formats the stateA data for json
				if (k == stateArray.size()-1){
					stateA = stateA+newTransition+"]";
				}else if(k < stateArray.size()){
					stateA = stateA+newTransition+",";}
			}//end of for loop that iterates through all states(k) at t+1
			//System.out.println(stateA);
			//the next 4 lines of code formats the allA data for json
			if (i == stateArray.size()-1){
				allA = allA+stateA+"]";
			}else if(i < stateArray.size()){
				allA = allA+stateA+",";}
		}//end of for all state i
		//System.out.println(allA);
		return allA;
	}//end newA()
	
	private String newBLogTest(){
		NodeState currNode;
		String currSeq;
		String allB = "\"B\":[";
		String stateB ="[";
		ObservationConverter oc;
		ArrayList <Integer> currEmission;
		double e;
		double eSum;
		double eSumSum;
		double exponent;
		double currLogGamma;
		double firstLogGamma = 0;
		double firstGamma = 0;
		double canFirstGamma;
		int firstGammaS = 0;
		int firstGammaT = 0;
		double newEmission;
		double numerator = 0;
		double denominator = 0;
		int currE; //current emission
				
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			/*
			The block of code below calculates the denominator over multiple sequences
			*/
			//begin find firstLogGamma
			for(int seqM = 0; seqM < sequences.size(); seqM++){
				currSeq = sequences.get(seqM); //gets a sequence to be analyzed as a string
				oc = new ObservationConverter();
				currEmission = oc.processObservation(currSeq, emissionSet); 
				//firstGammaS = 0;
				//firstGammaT = 0;
				for(int t = 0; t < currEmission.size()-1; t++){
					currLogGamma = currNode.getLogGamma(seqM, t);
					if(seqM == 0 && t == 0){
						firstGamma = currLogGamma;
						firstGammaS = 0;
						firstGammaT = 0;
					}else{
						canFirstGamma = currLogGamma;
						if(canFirstGamma > firstGamma){
							firstGamma = canFirstGamma;
							firstGammaS = seqM;
							firstGammaT = t;
						}
					}
				}
			}
			firstLogGamma = currNode.getLogGamma(firstGammaS, firstGammaT);
			//end find firstLogGamma for denominator
			/*
			The block of code below calc the denominator(logGammaSumSum)
			*/
			eSumSum = 0;
			for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
				currSeq = sequences.get(seqIndex); //gets a sequence to be analyzed as a string
				oc = new ObservationConverter();
				currEmission = oc.processObservation(currSeq, emissionSet); 
				eSum = 0;
				for(int t = 0; t < currEmission.size(); t++){
					currLogGamma = currNode.getLogGamma(seqIndex, t);
					//System.out.println("before minus first, at t = "+t+" es = "+currEmission.get(t)+" for i = "+i+" the currGamma is: "+currLogGamma);
					if(seqIndex != firstGammaS || t != firstGammaT){
						exponent = currLogGamma - firstLogGamma;
						////////////write QC CODE HERE
						
						//System.out.println("After minus first, at t = "+t+" es = "+currEmission.get(t)+" for i = "+i+" the currGamma is: "+currLogGamma);
						//System.out.println("firstT: "+firstGammaT);
						
						
						e = Math.exp(exponent);
						if(e == Double.NEGATIVE_INFINITY){
							e = 0;
							System.out.println("NEGATIVE INFINITY");
						}
						if(e == Double.POSITIVE_INFINITY){
							System.out.println("POSITIVE INFINITY FROM NEW B");
						}
						if(e == Double.NaN){
							System.out.println("NOT A NUMBER");
						}
						eSum = eSum + e;
					}
				}
				eSumSum = eSumSum + eSum;
			}
			denominator = firstLogGamma+Math.log(1+eSumSum);
			//System.out.println("The denominator is: "+denominator);
			
			
			
			/*
			The block of code below calculates the numerator over multiple sequences
			*/
			stateB ="[";//wipes previous state's data
			for (int es =0; es <emissionSet.size(); es++){//iterate through the emission set
				//begin find firstLogSigma
				boolean match = false;// false if no instance of a match for the emission
				double firstSigma = Double.NEGATIVE_INFINITY;
				double canSigma;
				int firstSigmaS = 0;
				int firstSigmaT = 0;
				for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
					currSeq = sequences.get(seqIndex);
					oc = new ObservationConverter();
					currEmission = oc.processObservation(currSeq, emissionSet);	
					for(int t = 0; t < currEmission.size(); t++){//iterates from t0 through T for currSeq
						currE=currEmission.get(t);
						if(currE == es){
							currLogGamma = currNode.getLogGamma(seqIndex, t);
							if(!match){
								match = true;
								firstSigma = currLogGamma;
								firstSigmaS = seqIndex;
								firstSigmaT = t;
							}else{
								canSigma = currLogGamma;
								if(canSigma > firstSigma){
									firstSigma = canSigma;
									firstSigmaS = seqIndex;
									firstSigmaT = t;
								}
							}
						}	
					}
				}
				//end find firstLogSigma
				
				if(match){
					firstLogGamma = currNode.getLogGamma(firstSigmaS, firstSigmaT);
					eSumSum = 0;
					for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
						currSeq = sequences.get(seqIndex);
						oc = new ObservationConverter();
						currEmission = oc.processObservation(currSeq, emissionSet);	
						eSum = 0;
						for(int t = 0; t < currEmission.size(); t++){//iterates from t0 through T for currSeq
							currE=currEmission.get(t);
							if(currE == es){
								currLogGamma = currNode.getLogGamma(seqIndex, t);
								//System.out.println("before minus first, at t = "+t+" es = "+es+" for i = "+i+" the currGamma is: "+currLogGamma);
								if(seqIndex != firstSigmaS || t != firstSigmaT){
									exponent = currLogGamma - firstLogGamma;
									//System.out.println("After minus first, at t = "+t+" es = "+es+" for i = "+i+" the currGamma is: "+currLogGamma);
									//System.out.println("firstT: "+firstSigmaT);
									e = Math.exp(exponent);
									if(e == Double.NEGATIVE_INFINITY){
										e = 0;
										System.out.println("NEGATIVE INFINITY");
									}
									if(e == Double.POSITIVE_INFINITY){
										System.out.println("POSITIVE INFINITY NEW B");
										System.out.println("currSigma: "+currLogGamma);
										System.out.println("firstSigma: "+firstLogGamma);
										System.out.println("exponent: "+exponent);
									}
									if(e == Double.NaN){
										System.out.println("NOT A NUMBER");
									}
									eSum = eSum+e;
								}
							}	
						}
						eSumSum = eSumSum +eSum;
					}//end for each sequence to be analyzed
					numerator = firstLogGamma+Math.log(1+eSumSum);
					newEmission = numerator - denominator;
					//System.out.println("numerator: "+numerator);
					//System.out.println("denominator: "+denominator);
					//System.out.println("newEmission: "+newEmission);
			}else{
				newEmission =-1000;//placeholder for zero probability
			}	
				//the next 4 lines of code formats the stateB data for json
				if (es == emissionSet.size()-1){
					stateB = stateB+newEmission+"]";
				}else if(es < emissionSet.size()){
					stateB = stateB+newEmission+",";}
			}//end for emission set
			//the next 4 lines of code formats the allB data for json
			if (i == stateArray.size()-1){
				allB = allB+stateB+"]";
			}else if(i < stateArray.size()){
				allB = allB+stateB+",";}
		}//iterates through state i
		//System.out.println("NewLogB: "+allB);
		return allB;
	}
	
	private String newBLog(){
		NodeState currNode;
		String currSeq;
		String allB = "\"B\":[";
		String stateB ="[";
		ObservationConverter oc;
		ArrayList <Integer> currEmission;
		double e;
		double eSum;
		double eSumSum;
		double exponent;
		double currLogGamma;
		double firstLogGamma = 0;
		double newEmission;
		double numerator = 0;
		double denominator = 0;
		int currE; //current emission
				
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			/*
			The block of code below calculates the denominator over multiple sequences
			*/
			eSumSum = 0;
			for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
				currSeq = sequences.get(seqIndex); 
				oc = new ObservationConverter();
				currEmission = oc.processObservation(currSeq, emissionSet);  
				eSum = 0;
				for(int t = 0; t < currEmission.size(); t++){
					currLogGamma = currNode.getLogGamma(seqIndex, t);
					if(seqIndex == 0 && t ==0){
						firstLogGamma = currLogGamma;
					}else{
						exponent = currLogGamma - firstLogGamma;
						e = Math.exp(exponent);
						eSum = eSum + e; // over the length of the sequence
					}
				}
				eSumSum = eSumSum +eSum;//over multiple sequences
			}
			denominator = firstLogGamma+Math.log(1+eSumSum);
			//System.out.println("The denominator is: "+denominator);
			/*
			The block of code below calculates the numerator over multiple sequences
			*/
			stateB ="[";//wipes previous state's data
			for (int es =0; es <emissionSet.size(); es++){//iterate through the emission set
				boolean first = true;
				eSumSum = 0;
				for(int seqIndex = 0; seqIndex < sequences.size(); seqIndex++){
					currSeq = sequences.get(seqIndex);
					oc = new ObservationConverter();
					currEmission = oc.processObservation(currSeq, emissionSet);	
					eSum = 0;
					for(int t = 0; t < currEmission.size(); t++){//iterates from t0 through T for currSeq
						currE=currEmission.get(t);
						if(currE == es){
							currLogGamma = currNode.getLogGamma(seqIndex, t);
							if(first){
								first = false;
								firstLogGamma = currLogGamma;
							}else{
								exponent = currLogGamma - firstLogGamma;
								e = Math.exp(exponent);
								eSum = eSum+e;
							}
						}	
					}
					eSumSum = eSumSum +eSum;
				}//end for each sequence to be analyzed
				numerator = firstLogGamma+Math.log(1+eSumSum);
				//System.out.println("The numerator is: "+numerator);
				//newEmission = Math.exp(numerator - denominator);
				newEmission = numerator - denominator;
				//the next 4 lines of code formats the stateB data for json
				if (es == emissionSet.size()-1){
					stateB = stateB+newEmission+"]";
				}else if(es < emissionSet.size()){
					stateB = stateB+newEmission+",";}
			}//end for emission set
			//the next 4 lines of code formats the allB data for json
			if (i == stateArray.size()-1){
				allB = allB+stateB+"]";
			}else if(i < stateArray.size()){
				allB = allB+stateB+",";}
		}//iterates through state i
		//System.out.println("NewLogB: "+allB);
		return allB;
	}//end newBLog()
	
	private String newB(){
		//MAY HAVE TO ACCOUNT FOR ZERO PROBABILITY EMISSION TO HAVE FLEXIBILITY
		NodeState currNode;
		String currSeq;
		String allB = "\"B\":[";
		String stateB ="[";
		ObservationConverter oc;
		ArrayList <Integer> currEmission;
		double gammaSum;
		double gammaSumSeq;
		double currGamma;
		double sigmaSum;
		double sigmaSumSeq = 0;
		double newEmission;
		int currE; //current emission
				
		for(int i = 0; i < stateArray.size(); i++){
			currNode = stateArray.get(i);
			/*
			The block of code below calculates gammaSumSeq for state i for all sequences
			to be analyzed, from t0 to T. 
			*/
			gammaSumSeq = 0;//sum of all gamma from state i for all sequences to be analyzed 
			for(int s = 0; s < sequences.size(); s++){
				currSeq = sequences.get(s); //gets a sequence to be analyzed as a string
				oc = new ObservationConverter();
				// the line below converts currSeq to an array of Integers
				currEmission = oc.processObservation(currSeq, emissionSet); 
				gammaSum = 0;//sum of all gamma from t0 to T for state i for sequence s 
				for(int j = 0; j < currEmission.size(); j++){//t0 to T to get gammaSum
					gammaSum = gammaSum + currNode.getGamma(s, j);
				}//end of for loop to get gammaSum
				gammaSumSeq = gammaSumSeq + gammaSum;
			}//end of for all sequences to be analyzed (to get gammaSumSeq)
			//System.out.println("The gammaSumSeq for state = "+i+" is "+gammaSumSeq);
			//end of block	
			/*
			The block of code below calculates sigmaSumSeq for state i for emission t
			*/
			stateB ="[";//wipes previous state's data
			for (int es =0; es <emissionSet.size(); es++){//iterate through the emission set
				sigmaSumSeq = 0;
				for(int s = 0; s < sequences.size(); s++){
					currSeq = sequences.get(s); //gets a sequence to be analyzed as a string
					oc = new ObservationConverter();
					// the line below converts currSeq to an array of Integers
					currEmission = oc.processObservation(currSeq, emissionSet);	
					sigmaSum = 0;
					for(int l = 0; l < currEmission.size(); l++){//iterates from t0 through T for sequence s
						currE=currEmission.get(l);
						currGamma = currNode.getGamma(s, l);
						if(currE == es)
							sigmaSum = sigmaSum+currGamma;
					}
					sigmaSumSeq = sigmaSumSeq+sigmaSum;

				}//end for each sequence to be analyzed
				//System.out.println("SigmaSumSeq is: "+sigmaSumSeq);
				//System.out.println("GammaSumSeq is: "+gammaSumSeq);
				newEmission = sigmaSumSeq/gammaSumSeq;
				//System.out.println("New emission for state:"+i+"for emission:"+es+" is "+newEmission);
				//the next 4 lines of code formats the stateB data for json
				if (es == emissionSet.size()-1){
					stateB = stateB+newEmission+"]";
				}else if(es < emissionSet.size()){
					stateB = stateB+newEmission+",";}
			}//end for emission set
			//the next 4 lines of code formats the allB data for json
			if (i == stateArray.size()-1){
				allB = allB+stateB+"]";
			}else if(i < stateArray.size()){
				allB = allB+stateB+",";}
		}//iterates through state i
		//System.out.println(allB);
		return allB;
	}//end newB()
}//end class
