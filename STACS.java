import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class STACS {	

	HMM currHMM;//the HMM currently loaded 

	private void menu(){
		
		boolean quit = false;
		Scanner select = new Scanner(System.in);
		
		while(!quit){
			System.out.println("What would you like to do?");
			System.out.println("(A) Upload HMM");
			System.out.println("(B) Test adjusted HMM");
			System.out.println("(C) Get HMM parameters");
			System.out.println("(D) Score sequences");
			System.out.println("(E) Most likely sequence of states");
			System.out.println("(F) Probability of being at a state given the Sequence");
			System.out.println("(G) Adjust Parameters (specify iterations)");
			System.out.println("(H) Adjust Parameters (until convergence)");
			System.out.println("(Q) Quit");
			System.out.print("Enter your selection:");
			String response = select.nextLine();
		
			if (response.equals("A")){
				processHMM();
			}
			if (response.equals("B")){
				buildHMM();
			}
			if (response.equals("C")){
				displayHMMParameters();
			}
			if(response.equals("D")){
				scoreSequences();
			}
			if(response.equals("E")){
				stateSequence();
			}
			if(response.equals("F")){
				stateProbability();
			}
			if(response.equals("G")){
				adjustParameters();
			}
			if(response.equals("H")){
				convergence();
			}
			if(response.equals("Q")){
				quit = true;
			}	
		}//end of while(!quit)		
	}//end of menu()
	
	
	
	private void processHMM(){
		
		Scanner fileName = new Scanner(System.in);
		System.out.println("File Name:");
		String file = fileName.nextLine();
		File hMM = new File(file);
		System.out.println("Building "+file+" HMM from file");
		try {
			Scanner scan = new Scanner(hMM);
			while(scan.hasNextLine()) {
			    String jsonHMM = scan.nextLine();
			    loadHMM(jsonHMM);
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//end processHMM()
	
	//METHOD USED FOR TESTING ADJUSTED HMMS
	private void processAHMM(){
		
		Scanner fileName = new Scanner(System.in);
		System.out.println("File Name:");
		String file = fileName.nextLine();
		File hMM = new File(file);
		System.out.println("Building "+file+" HMM from file");
		try {
			Scanner scan = new Scanner(hMM);
			while(scan.hasNextLine()) {
			    String jsonHMM = scan.nextLine();
			    loadHMM(jsonHMM, false);
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}//end processHMM()
	
	private void buildHMM(){
		processAHMM();
	}//end buildHMM()
	
	private void displayHMMParameters(){
		System.out.println("The emission set is: { "+currHMM.getEmissionSet()+"}");
		System.out.println("The initial transition parameters are: { "+currHMM.getPI()+"}");
		System.out.println("The transition probabilities are: { "+currHMM.getA()+"}");
		System.out.println("The emission probabilities are: { "+currHMM.getB()+"}");
	}//end getHMMParameters()
	
	private void scoreSequences(){
		System.out.print("Enter the name of the file containing the sequences:");
		Scanner file = new Scanner(System.in);
		String fileName = file.nextLine();
		ForwardAlgo fa = new ForwardAlgo(currHMM, fileName);
		fa.displayScores();
		BackwardAlgo ba = new BackwardAlgo(currHMM, fileName);
		ba.displayScores();
	}//end scoreSequences()
	
	private void scoreSequences(String fileName){//this version is used in Baum-Welch iterations
		ForwardAlgo fa = new ForwardAlgo(currHMM, fileName);
		fa.displayScores();
		//BackwardAlgo ba = new BackwardAlgo(currHMM, fileName);
		//ba.displayScores();
		
	}//end scoreSequences()
	
	private void stateSequence(){
		//display score and sequence of states
		System.out.print("Enter the name of the file containing the sequences:");
		Scanner file = new Scanner(System.in);
		String fileName = file.nextLine();
		ForwardAlgo fa = new ForwardAlgo(currHMM, fileName);
		ArrayList<String>stateSeqArray = fa.getStateSeq();
		for(int i = 0; i < stateSeqArray.size(); i++){
			System.out.println(stateSeqArray.get(i));
		}
		
	}//end stateSequence()
	
	private void stateProbability(){
		ForwardBackward fb = new ForwardBackward(currHMM, "file name");
		Scanner sequence = new Scanner(System.in);
		Scanner time = new Scanner(System.in);
		System.out.println("Which sequence would you like to look at?");
		String seq = sequence.nextLine(); 
		System.out.println("What time point would you like to see?");
		String tm = time.nextLine();
		int s = Integer.parseInt(seq);
		int t = Integer.parseInt(tm);
		String sp = fb.getStateProbabilities(s,t);
		System.out.println(sp);
	}
	
	private void adjustParameters(){
		
		System.out.print("Enter the name of the file containing the training sequences:");
		Scanner file = new Scanner(System.in);
		String fileName = file.nextLine();
		System.out.println("How many sessions with the training sequences?");
		Scanner sessions = new Scanner(System.in);
		String times = sessions.nextLine();
		int iterations = Integer.parseInt(times);
		scoreSequences(fileName);
		for(int i = 0; i < iterations; i++){
			BaumWelchAlgo bw = new BaumWelchAlgo(currHMM, fileName);
			String newHMM = bw.adjustParameters();
			loadHMM(newHMM, false);//false indicates that it is an adjusted hmm
			scoreSequences(fileName);
		}
	}
	
	private void convergence(){
		
		System.out.print("Enter the name of the file containing the training sequences:");
		Scanner file = new Scanner(System.in);
		String fileName = file.nextLine();
		scoreSequences(fileName);
		boolean convergence = false;
		while(!convergence){
			BaumWelchAlgo bw = new BaumWelchAlgo(currHMM, fileName);
			String newHMM = bw.adjustParameters();
			loadHMM(newHMM, false);//false indicates that it is an adjusted hmm
			scoreSequences(fileName);
			convergence = determineConvergence();
		}
	}
	private boolean determineConvergence(){
		boolean convergence = false;
		
		return convergence;
	}
	
	private void loadHMM(String json, boolean ori){
		
		//System.out.println(json);
		currHMM = new HMM(json, ori);
	}//end loadHMM(String json)
	
	private void loadHMM(String json){
		
		//System.out.println(json);
		currHMM = new HMM(json);
	}//end loadHMM(String json)
	
	public static void main(String[] args ){
		
		STACS ns = new STACS();
		double a = -999999999+-9999999;
		double b = -9+-9;
		double e = Math.exp(a-b);
		System.out.println(e);
		if(e == Double.NEGATIVE_INFINITY){
			System.out.println("NEGATIVE INFINITY");
		}
		if(e == Double.POSITIVE_INFINITY){
			System.out.println("POSITIVE INFINITY");
		}
		ns.menu();
	}//end main()
}
