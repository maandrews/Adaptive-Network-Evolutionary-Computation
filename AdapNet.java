import java.io.*;
import java.util.*;

/*
 * This program uses evolutionary computation to simulate an adaptive network (graph) where its nodes are modelling a 
 * population vulnerable to disease.  Nodes may use different strategies to rewire their connections within the network
 * in order to avoid becoming infected.  The best strategies are replicated from generation to generation, in order to find
 * the optimal strategy that maximizes influence in the network.
 */

public class AdapNet {
	
	
	public static void main(String[] args) throws IOException {
		
		File data = new File("NetValues.txt");
		PrintWriter pwInput = new PrintWriter(data);
		
		final double INF = Double.POSITIVE_INFINITY;
		final int SIZE = 200;
		final int TIMESTEPS = 2000; // length of simulation
		final int GENLEN = 100; // generation length
		
		Random ranGen = new Random();
		
		double [][] network = new double [SIZE][SIZE]; // adjacency matrix
		double [][] networkFW = new double [SIZE][SIZE]; //shortest path network
		int [][] popStatus = new int [SIZE][2]; // Shows infection status, 1=not inf, -1=inf, and how long someone has been infected
		int [][] numcontacts = new int [SIZE][2]; // degree of each node
		double [][] centrality = new double [SIZE][2]; // centrality of each node
		
		final int geneLen = 1; // Strategy set length
		final int NUMSTRATS = 5; // number of strategies used
		int strategies [][] = new int [SIZE][1]; // array of strategies
		double fitness [][] = new double [SIZE][2]; // array of fitnesses and which nodes have which fitness
		
		
		double [] prev = new double [TIMESTEPS+1]; // prevalence
		double [] s0 = new double [(TIMESTEPS/GENLEN)+1];
		double [] s1 = new double [(TIMESTEPS/GENLEN)+1];
		double [] s2 = new double [(TIMESTEPS/GENLEN)+1];
		double [] s3 = new double [(TIMESTEPS/GENLEN)+1];
		double [] s4 = new double [(TIMESTEPS/GENLEN)+1];
		
		final double TRATE = 0.03; // transmission rate
		final double RRATE = 0.05; // rewire rate
		final double MRATE = 0.02; // mutation rate
		
		double totalInf = 0;
		
		/* Strategies:
		 * 
		 * 0 = random node
		 * 1 = highest degree node
		 * 2 = lowest degree node
		 * 3 = highest centrality node
		 * 4 = lowest centrality node
		 * 
		 */
		
		for(int i=0 ; i<network.length ; i++){
			fitness[i][0] = i;
			fitness[i][1] = 0;
			
			strategies[i][0] = ranGen.nextInt(NUMSTRATS);
			
			for(int j=0 ; j<network.length ; j++){
				network[i][j] = 0;
			}
		}
		
		// Constructing random network
		for(int i=0 ; i<SIZE ; i++){
			popStatus[i][0] = 1;
			popStatus[i][1] = 0;
			for(int j=0 ; j<SIZE ; j++){
				double r = Math.random();

				if(i != j && network[i][j] != 1){
					if(r < 0.025){network[i][j] = 1;network[j][i] = 1;}
				}
			}

		}
		
		// initial shortest paths
		for(int i=0 ; i<SIZE ; i++){
			for(int j=0 ; j<SIZE ; j++){
				if(j==i){networkFW[i][j] = 0;}
				else if (network[i][j]==0){networkFW[i][j]=INF;}
				else{networkFW[i][j] = network[i][j];}
			}
		}
		networkFW = FWAlg(networkFW);
		
		
		// Printing out degree of each node if desired
		/*int [] numcontacts1 = new int [SIZE];
		for(int i=0 ; i<SIZE ; i++){numcontacts1[i] =0;}
		
		for(int i=0 ; i<SIZE ; i++){
			for(int j=0 ; j<SIZE ; j++){	
				if(network[i][j] == 1){numcontacts1[i]++;}	
			}
		}
		
		System.out.print("O=[");
		for(int i = 0 ; i < SIZE ; i++){
			System.out.print(numcontacts1[i]);
			if(i != SIZE-1){System.out.print("; ");}
		}
		System.out.println("];");
	
		//System.exit(0);*/

		
		while(totalInf < 5){ // infecting initial people
			
			int i = ranGen.nextInt(SIZE);
			while(popStatus[i][0] < 0){
				i = ranGen.nextInt(SIZE);
			}
			popStatus[i][0] = -1;
			popStatus[i][1] = 1;
			totalInf++;
		}
		
		
		for(int i=0 ; i<SIZE ; i++){		
			if(strategies[i][0]==0){s0[0] = s0[0] + 1;}
			else if(strategies[i][0]==1){s1[0] = s1[0] + 1;}
			else if(strategies[i][0]==2){s2[0] = s2[0] + 1;}
			else if(strategies[i][0]==3){s3[0] = s3[0] + 1;}
			else if(strategies[i][0]==4){s4[0] = s4[0] + 1;}	
		}
		
		// Setting prevalence and strategy usage among population
		prev[0] = totalInf / (double)SIZE;
		s0[0] = s0[0] / (double)SIZE;
		s1[0] = s1[0] / (double)SIZE;
		s2[0] = s2[0] / (double)SIZE;
		s3[0] = s3[0] / (double)SIZE;
		s4[0] = s4[0] / (double)SIZE;
		
		
		
		
		for(int i=0 ; i<SIZE ; i++){ // finding degree of each node
			numcontacts[i][0] = i;
			numcontacts[i][1] = 0;
			for(int j=0 ; j<SIZE ; j++){	
				if(network[i][j] == 1){numcontacts[i][1]++;}	
			}
		}
		
		for(int i = 0 ; i < SIZE ; i++){ // finding centrality of each node
			centrality[i][0] = i;
			centrality[i][1] = 0;
			for(int j = 0 ; j < SIZE ; j++){
				if(i != j && networkFW[i][j] != INF && networkFW[i][j] != 0){
					centrality[i][1] += (1 / networkFW[i][j]);
				}
			}
		}
		
		
		
		
		for(int step = 1 ; step < TIMESTEPS ; step++){
			
			
			for(int i = 0 ; i < SIZE ; i++){ // Update fitness			
				for(int j = 0 ; j < SIZE ; j++){
					
					if(i != j && networkFW[i][j] != INF && networkFW[i][j] != 0){
						fitness[i][1] += (1 / networkFW[i][j]);
					}
					
				}				
			}
			
			
			// End generation, evaluate fitness and change strategies/apply mutations
			if(step % GENLEN == 0){
				
				Arrays.sort(fitness, new Comparator<double[]>() { // order fitnesses best to worst
					@Override public int compare(double[] o1, double[] o2) {
						return Double.valueOf(o2[1]).compareTo(o1[1]);
					}
				});
				
				for(int i = 0 ; i < 50 ; i++){ // replace worst strategies with best strategies
					
					// giving people a new strategy
					strategies[(int)fitness[SIZE-i-1][0]][0] = strategies[(int)fitness[i][0]][0];
					
				}
				
				for(int i = 0 ; i < SIZE ; i++){ // mutate
					double r = Math.random();
					if(r < MRATE){
						int nS = ranGen.nextInt(NUMSTRATS);
						while(nS == strategies[i][0]){nS = ranGen.nextInt(NUMSTRATS);}
						strategies[i][0] = nS;
					}	

				}
				
				// Reset fitness
				for(int i=0 ; i<SIZE ; i++){
					fitness[i][0] = i;
					fitness[i][1] = 0;
				}
				
				
				// Updating prevalence of strategies amongst population
				for(int i=0 ; i<SIZE ; i++){		
					if(strategies[i][0]==0){s0[step/GENLEN] = s0[step/GENLEN] + 1;}
					else if(strategies[i][0]==1){s1[step/GENLEN] = s1[step/GENLEN] + 1;}
					else if(strategies[i][0]==2){s2[step/GENLEN] = s2[step/GENLEN] + 1;}
					else if(strategies[i][0]==3){s3[step/GENLEN] = s3[step/GENLEN] + 1;}
					else if(strategies[i][0]==4){s4[step/GENLEN] = s4[step/GENLEN] + 1;}	
				}
				s0[step/GENLEN] = s0[step/GENLEN] / (double)SIZE;
				s1[step/GENLEN] = s1[step/GENLEN] / (double)SIZE;
				s2[step/GENLEN] = s2[step/GENLEN] / (double)SIZE;
				s3[step/GENLEN] = s3[step/GENLEN] / (double)SIZE;
				s4[step/GENLEN] = s4[step/GENLEN] / (double)SIZE;
				
				
			} // end generation update
			

			
			// finding/ordering highest/lowest degree node.
			Arrays.sort(numcontacts, new Comparator<int[]>() { // order degrees highest to lowest
				@Override public int compare(int[] o1, int[] o2) {
					return Integer.valueOf(o2[1]).compareTo(o1[1]);
				}
			}); 
			// end finding highest/lowest degree node
			
			// finding/ordering highest/lowest centrality node.
			Arrays.sort(centrality, new Comparator<double[]>() { // order centralities highest to lowest
				@Override public int compare(double[] o1, double[] o2) {
					return Double.valueOf(o2[1]).compareTo(o1[1]);
				}
			}); 
			// end ordering highest/lowest centrality node
						
			////////////////////////////////// Rewiring /////////////////////////////////////////
			for(int i = 0 ; i < SIZE ; i++){  // rewire

				if(popStatus[i][0] > 0){ // if someone is not infected

					for(int j = 0 ; j < SIZE ; j++){ 	// they may rewire

						if(network[i][j] > 0 && popStatus[j][1] > 0){ // if j is an infectious contact

							double r = Math.random();
							if(r < RRATE){	// person rewires
								
								boolean rewired = false;
								while(!rewired){
									
									switch(strategies[i][0]){
									
									case 0: 
										int newC = ranGen.nextInt(SIZE);
										while(network[i][newC] > 0 || popStatus[newC][1] > 0 || newC == i){
											newC = ranGen.nextInt(SIZE);
										}
										network[i][j] = 0;
										network[j][i] = 0;
										network[i][newC] = 1;
										network[newC][i] = 1;
										rewired = true;
										break;
									// end case 0
										
									case 1:
										System.out.println(step);
										rewired = true;
										int k = 0;
										while(k < SIZE){
											if(network[i][numcontacts[k][0]] == 0 && popStatus[numcontacts[k][0]][1] == 0 && numcontacts[k][0] != i){  
												network[i][j] = 0;
												network[j][i] = 0;
												network[i][numcontacts[k][0]] = 1;
												network[numcontacts[k][0]][i] = 1;
												
												k = SIZE;
											}
											k++;
										}
										break;
										
									//end case 1
										
									case 2:
										rewired = true;
										k = SIZE-1;
										while(k > 0){
											if(network[i][numcontacts[k][0]] == 0 && popStatus[numcontacts[k][0]][1] == 0 && numcontacts[k][0] != i){  
												network[i][j] = 0;
												network[j][i] = 0;
												network[i][numcontacts[k][0]] = 1;
												network[numcontacts[k][0]][i] = 1;
												
												k = -1;
											}
											k--;
										}
										break;
										
									//end case 2
										
									case 3:
										rewired = true;
										k = 0;
										while(k < SIZE){
											if(network[i][(int)centrality[k][0]] == 0 && popStatus[(int)centrality[k][0]][1] == 0 && (int)centrality[k][0] != i){  
												network[i][j] = 0;
												network[j][i] = 0;
												network[i][(int)centrality[k][0]] = 1;
												network[(int)centrality[k][0]][i] = 1;
												
												k = SIZE;
											}
											k++;
										}
										break;
										
									//end case 3
										
									case 4:
										rewired = true;
										k = SIZE-1;
										while(k > 0){
											if(network[i][(int)centrality[k][0]] == 0 && popStatus[(int)centrality[k][0]][1] == 0 && (int)centrality[k][0] != i){  
												network[i][j] = 0;
												network[j][i] = 0;
												network[i][(int)centrality[k][0]] = 1;
												network[(int)centrality[k][0]][i] = 1;
												
												k = -1;
											}
											k--;
										}
										break;
										
									//end case 4
										
									} // end switch
									
								}
								
							
							} // end rewire

						} // end if infectious loop

					}

				} // end if person rewiring loop

			} // end rewire people
					
			
			//////////////////////////// Infection //////////////////////////////////////////////
			for(int i = 0 ; i < SIZE ; i++){  // Infect people

				if(popStatus[i][0] < 0 && popStatus[i][1] > 0){ // if someone is infected

					for(int j = 0 ; j < SIZE ; j++){ 			// they may infect contacts

						if(network[i][j] > 0 && popStatus[j][0] > 0){

							double r = Math.random();
							if(r < TRATE){	// person becomes infected
								popStatus[j][0] = -1;
							} 

						} // end if

					}

				} // end if

			} // end infect people ///////////////////////////////////////////////
			
			totalInf = 0;

			for(int i = 0 ; i < SIZE ; i++){  // Allowing people to potentially recover and determining prevalence 

				if(popStatus[i][0] < 0){
					popStatus[i][1]++;
					if(popStatus[i][1] > 5){popStatus[i][0] = 1; popStatus[i][1] = 0;}
					else{totalInf++;}
				}

			}
			
			prev[step] = totalInf / (double)SIZE;
			
			if(totalInf == 0){ // if no one is infected anymore...
				
				while(totalInf < 5){ // infect 5 more people
					
					int i = ranGen.nextInt(SIZE);
					while(popStatus[i][0] < 0){
						i = ranGen.nextInt(SIZE);
					}
					popStatus[i][0] = -1;
					popStatus[i][1] = 1;
					totalInf++;
				}
				
			}
			
			
			// Finding shortest paths
			for(int i=0 ; i<SIZE ; i++){
				for(int j=0 ; j<SIZE ; j++){
					if(j==i){networkFW[i][j] = 0;}
					else if (network[i][j]==0){networkFW[i][j]=INF;}
					else{networkFW[i][j] = network[i][j];}
				}
			}
			networkFW = FWAlg(networkFW);
		

			// rewires that make new random connections
			double r = Math.random();
			if(r < 0.01){
				int n1 = ranGen.nextInt(SIZE);
				while(numcontacts[n1][1] == 0){
					n1 = ranGen.nextInt(SIZE);
				}
				int n2 = ranGen.nextInt(SIZE);
				while(n2 == n1 || network[n1][n2] > 0){
					n2 = ranGen.nextInt(SIZE);
				}
				int n3 = ranGen.nextInt(SIZE);
				while(n1 == n3 || n2 == n3 || network[n1][n3] == 0){
					n3 = ranGen.nextInt(SIZE);
				}
				network[n1][n3] = 0;
				network[n3][n1] = 0;
				network[n1][n2] = 1;
				network[n2][n1] = 1;
			}
			
			for(int i=0 ; i<SIZE ; i++){ // finding degree of each node
				numcontacts[i][0] = i;
				numcontacts[i][1] = 0;
				for(int j=0 ; j<SIZE ; j++){	
					if(network[i][j] == 1){numcontacts[i][1]++;}	
				}
			}
			
			for(int i = 0 ; i < SIZE ; i++){ // finding centrality of each node
				centrality[i][0] = i;
				centrality[i][1] = 0;
				for(int j = 0 ; j < SIZE ; j++){
					if(i != j && networkFW[i][j] != INF && networkFW[i][j] != 0){
						centrality[i][1] += (1 / networkFW[i][j]);
					}
				}
			}
			
			
			
		} // end timesteps loop
		
		
		//Printing out degree of each node if desired
		/*int [] numcontacts1 = new int [SIZE];
		for(int i=0 ; i<SIZE ; i++){numcontacts1[i] =0;}
		
		for(int i=0 ; i<SIZE ; i++){
			for(int j=0 ; j<SIZE ; j++){	
				if(network[i][j] == 1){numcontacts1[i]++;}	
			}
		}
		
		pwInput.print("O=[");
		for(int i = 0 ; i < SIZE ; i++){
			pwInput.print(numcontacts1[i]);
			if(i != SIZE-1){pwInput.print("; ");}
		}
		pwInput.println("];");*/
		
		
		
		
		// printing some info
		pwInput.print("I=["); // print infection pevalence
		for(int i = 0 ; i < TIMESTEPS ; i++){
			
			pwInput.print(prev[i] + "; ");
			if(i == TIMESTEPS-1){pwInput.print(prev[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("s0=["); // strat 0 prevalence
		for(int i = 0 ; i < (TIMESTEPS/GENLEN)+1 ; i++){
			
			pwInput.print(s0[i] + "; ");
			if(i == (TIMESTEPS/GENLEN)){pwInput.print(s0[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("s1=["); // strat 1 prevalence
		for(int i = 0 ; i < (TIMESTEPS/GENLEN)+1 ; i++){
			
			pwInput.print(s1[i] + "; ");
			if(i == (TIMESTEPS/GENLEN)){pwInput.print(s1[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("s2=["); // strat 2 prevalence
		for(int i = 0 ; i < (TIMESTEPS/GENLEN)+1 ; i++){
			
			pwInput.print(s2[i] + "; ");
			if(i == (TIMESTEPS/GENLEN)){pwInput.print(s2[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("s3=["); // strat 3 prevalence
		for(int i = 0 ; i < (TIMESTEPS/GENLEN)+1 ; i++){
			
			pwInput.print(s3[i] + "; ");
			if(i == (TIMESTEPS/GENLEN)){pwInput.print(s3[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("s4=["); // strat 4 prevalence
		for(int i = 0 ; i < (TIMESTEPS/GENLEN)+1 ; i++){
			
			pwInput.print(s4[i] + "; ");
			if(i == (TIMESTEPS/GENLEN)){pwInput.print(s4[i]);}
			
		}
		pwInput.println("];");
		
		pwInput.print("t=[");
		for(int i = 0 ; i < s0.length ; i++){
			
			pwInput.print(i + "; ");
			if(i == s0.length-1){pwInput.print(i);}
			
		}
		pwInput.println("];");
		
		pwInput.print("t=[");
		for(int i = 0 ; i < TIMESTEPS ; i++){
			
			pwInput.print(i + "; ");
			if(i == TIMESTEPS-1){pwInput.print(i);}
			
		}
		pwInput.println("];");
		
		
	}
	
	
	public static double[][] FWAlg(double[][] adjMatrix) {
		
		for(int k=0;k<adjMatrix.length;k++){
			for(int i=0;i<adjMatrix.length;i++){
				for(int j=0;j<adjMatrix.length;j++){
					adjMatrix[i][j]=Math.min(adjMatrix[i][j],adjMatrix[i][k]+adjMatrix[k][j]);
				}
			}
		}
		
		return adjMatrix;
	}
	
	
	static void shuffleArray(int[] ar)
	  {
	    Random rnd = new Random();
	    for (int i = ar.length - 1; i >= 0; i--)
	    {
	      int index = rnd.nextInt(i + 1);
	      
	      int a = ar[index];
	      ar[index] = ar[i];
	      ar[i] = a;
	    }
	  }

}




