import java.io.*;
import java.lang.Math.*;
import java.util.Arrays;
public class PWMModell{
	static double pseudocount = 1;
	static int PWMLength = 30;
	static char[] muster = {'A','T','G','C'};
	static double[] p = {0.27,0.25,0.25,0.23};
	static int sequenceLength = 200;
	static Score[] scoreArray;
	static double[] threshold;
	static int validCount;
	static int invalidCount;
	static boolean codonIncluded = false;
	static int countA = 0;
	static int countT = 0;
	static int countC = 0;
	static int countG = 0;
	
	static double[][] matrix;
	
	public static void main (String[] args){
		threshold = new double[10];
		if(args.length != 0){
			PWMLength = Integer.parseInt(args[0]);
			codonIncluded = true;
		}
		matrix = new double[muster.length][PWMLength];
		
		String name = "TIS-Ecoli.txt";
		File file = new File(name);
		initPWM(file);
		//Zählt die Anzahl der Token, dann erstellet ein Feld der entsprechenden Länge
		//Dann initialisert das Feld mit den jew. scores der Codons
		int temp = initScoreArray(file);
		Arrays.sort(scoreArray);
		//for(int i = 1; i<10; i++)
		System.out.println("Threshold bei 50% ist ca " + String.format("%.2f", getThreshold(50)));
		//System.out.println("ValidTokens: " + cvt);
		System.out.println("TotalTokens: " + countStartCodons(file));
		
		
		
		name = "Trainingssequenzen.txt";
		file = new File(name);
		//Zählt die Anzahl der Token, dann erstellet ein Feld der entsprechenden Länge
		//Dann initialisert das Feld mit den jew scores der Codons
		temp = initScoreArray(file);
		for(int i = 1; i<10; i++){
			System.out.println("Threshold bei " + i*10 + "% ist ca " + String.format("%.2f", getThreshold(i*10)));
			threshold[i-1] = getThreshold(i*10);
			}
		
		
		name = "Testsequenzen.txt";
		file = new File(name);

		//Zählt die Anzahl der Token, dann erstellet ein Feld der entsprechenden Länge
		//Dann initialisert das Feld mit den jew scores der Codons
		temp = initScoreArray(file);
		for(int i = 1; i<10; i++){
			countValid(threshold[i-1]);
			System.out.println("validCount : " + validCount + " invalidCount : " + invalidCount);
		}
		countChars(new File("TIS-Ecoli.txt"));
		
		int sum = countA + countT + countC + countG;
		System.out.println("countA: " + ((double)countA / sum));
		System.out.println("countT: " + ((double)countT / sum));
		System.out.println("countC: " + ((double)countC / sum));
		System.out.println("countG: " + ((double)countG / sum));
		return;
		
	}
	
	//Zählt die Startcodons in der TIS-Ecoli.txt Datei
	private static int countStartCodons(File file){
		int sum = 0;
		String line;
		try (BufferedReader br = new BufferedReader(new FileReader(file))){
				//System.out.println("Anzahl der Startcodons: " + countStartCodons(file, br));
				while((line = br.readLine()) != null){
					//Zähle die möglichen Startcodons pro Zeile zusammen
					sum += checkLine(line);
				}
		} catch (IOException e){}
		return sum;
	}
	
	//Zählt die möglichen Startcodons in einer Sequenz
	private static int checkLine(String line){
		int sum = 0;
		//Startcodon = [A,T,G]TG
		for(int i=0;i<line.length()-2;i++){
			if(line.charAt(i) != 'C' && line.charAt(i+1) == 'T' && line.charAt(i+2) == 'G'){
				sum++;
				//System.out.println(score(i, line.charAt(i))+ " " + i);
			}
		}
		return sum;
	}
	
	private static void initPWM(File file){
		int lineCount = 0;
		String line="";
		try(BufferedReader br = new BufferedReader(new FileReader(file))){
			initPseudoCount();

				//Wertet jedes Zeichen der Sequenz aus und erstellt die PFM
				//System.out.println(line);
				while((line = br.readLine()) !=null){
					lineCount++;
					for(int j = 71; j<101; j++){
						int i=j-71;
						switch(line.charAt(j)){
							case 'A': matrix[0][i] +=1.0;
								break;
							case 'T': matrix[1][i] +=1.0;
								break;
							case 'G': matrix[2][i] +=1.0;
								break;
							case 'C': matrix[3][i] +=1.0;
								break;
						}
					}
				}
		}catch(IOException e){}
		//PFM -> PPM -> PWM
		PFMtoPPM(lineCount);
		PPMtoPWM();
	}
	
	
	//Verarbeitet die PFM zu einer PPM
	private static void PFMtoPPM(int lineCount){
		for(int i=0;i<4;i++)
			for(int j=0;j<PWMLength;j++)
				matrix[i][j] = matrix[i][j]/(double)lineCount;
			
	}
	
	//Verarbeitet die PPM zu einer PWM
	private static void PPMtoPWM(){
		for(int i=0;i<4;i++)
			for(int j=0;j<PWMLength;j++)
				matrix[i][j] = (Math.log(matrix[i][j]/p[i])/Math.log(2));	//log2(a[ij]/0,25)
	}
	
	//Addiert den Pseudocount zu jedem eintrag der PFM
	private static void initPseudoCount(){
		for(int i=0;i<4;i++)
			for(int j=0;j<PWMLength;j++)
				matrix[i][j] += pseudocount;
	}
	
	//Gibt die invertierte Matrix aus mit abschneiden überflüssiger Kommastellen
	private static void printMatrix(){
		for(int j=0;j<PWMLength;j++){
			for(int i=0;i<4;i++){
				System.out.print(String.format("%.2f", matrix[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	//Zählt alle den Kriterien entsprechenden Startcodons in einer Datei
	private static int countValidToken(File file,  boolean initScoreArray){
		int sum = 0;
			String line;
			try (BufferedReader br = new BufferedReader(new FileReader(file))){
					//System.out.println("Anzahl der Startcodons: " + countStartCodons(file, br));
					while((line = br.readLine()) != null)
						//Zähle die möglichen Startcodons pro Zeile zusammen
						sum += checkLineForValid(line, initScoreArray, sum);
			} catch (IOException e){}
			return sum;
	}
	
	//Zählt alle den Kriterien entsprechenden Startcodons in einer Zeile
	private static int checkLineForValid(String line, boolean initScoreArray, int temp){
		double value = 0;
		int sum = 0;
		//PWMLength damit immer mind 30 Basen vor dem Codon beobachtbar sind
		for(int i=PWMLength;i<line.length()-2;i++){
			if(line.charAt(i) != 'C' && line.charAt(i+1) == 'T' && line.charAt(i+2) == 'G'){
				value = score(i, line);
				if(initScoreArray)
					scoreArray[sum+temp] = new Score(codonIncluded ? value +1 : value, i == 100);
				sum++;
			}
		}
		return sum;
	}
	
	//Berechnet den Score eines codons durch multiplikation der einzelscores
	private static double score(int n, String line){
		double score=0;
		int start = n-PWMLength;
		for(int i = start; i<n; i++){
			switch(line.charAt(i)){
				case 'A': score += matrix[0][i-start];
					break;
				case 'T': score += matrix[1][i-start];
					break;
				case 'G': score += matrix[2][i-start];
					break;
				case 'C': score += matrix[3][i-start];
					break;
			}
		}
		
		return score;
	}
	
	//ermittelt die grenze
	private static double getThreshold(double percentage) {
		//abfangen falls kein prozentwert übergeben
		if(percentage > 100 || percentage < 0)
			System.exit(1);
		
		int validCount = 0;
		for (int i =scoreArray.length-1; i>=0; i--) {
			if(scoreArray[i].getValid())
				validCount++;
		}
		int secondCount = 0;
		int position = 0;
		for (int i =scoreArray.length-1; i>=0; i--) {
			if(scoreArray[i].getValid())
				secondCount++;
			if(((double)secondCount / (double)validCount) >= (percentage / 100)){
				position = i;
				break;
			}
		}
			
		System.out.println("Feler " + (scoreArray.length - position));
		return scoreArray[position].getScore();
	}
	
	private static void countValid(double threshold){
		validCount = 0;
		invalidCount = 0;
		int i = scoreArray.length-1;
		while(threshold <= scoreArray[i].getScore()){
			if(scoreArray[i].getValid())
				validCount++;
			else
				invalidCount++;
			i--;
		}	
	}
	
	private static int initScoreArray(File file){
		int cvt = countValidToken(file, false);
		scoreArray = new Score[cvt];
		int temp = countValidToken(file, true);
		Arrays.sort(scoreArray);
		return cvt;
	}
	
	private static void countChars(File file){
		String line = "";
		try(BufferedReader br = new BufferedReader(new FileReader(file))){

				//Wertet jedes Zeichen der Sequenz aus und erstellt die PFM
				//System.out.println(line);
				while((line = br.readLine()) !=null){
					
					for(int j = 0; j<line.length(); j++){
						switch(line.charAt(j)){
							case 'A': countA++;
								break;
							case 'T': countT++;
								break;
							case 'G': countG++;
								break;
							case 'C': countC++;
								break;
						}
					}
				}
		}catch(IOException e){}
		return;
	}
}






























