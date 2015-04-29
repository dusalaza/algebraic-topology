import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Scanner;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.lang.StringBuilder;
import java.io.BufferedWriter;

/**    
*    Furnace Foam Evaluator
*    Copyright (C) 2015  Dido Salazar
*
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License along
*    with this program; if not, write to the Free Software Foundation, Inc.,
*    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*
*   
*/


public class KhovanovJ{
  
  public static int[] mlrab = {0,0,0,0,0};
  public static final int MLRLEN = mlrab.length - 2;
  public static final int DOTS = mlrab.length;
  public static String[] monomial = {"X^{", "Y^{", "Z^{", "U^{", "V^{" };

  //Checks if an array is a cyclic permutation of another (inplace O(n) algorithm)
  public static boolean isRot(int[] a, int[] b){

    int n = a.length; //length of array
    
    //return false if the sizes don’t match
    if(n != b.length){
      return false;
    }
    
    //indices and counters
    
    int i = -1;
    int j = -1;
    int matchp1 = 1; //number of matches + 1
    
    /* The following loop begins by searching for the index of the first element of b inside a.
     * Once it finds it, the inner loop compares the consecutive entries of a to b and counts the matches. 
     * If they all match, the function returns. If there is a mismatch, the bigger loop breaks and the function returns false.
    */
    while((i < n-1) && (j < n-1)){
    
      //track the number of consecutive matches
      while((matchp1 <= n) && (a[(i + matchp1) % n] == b[(j + matchp1) % n])){
        matchp1++;
      } 
      
      //if all match, then return
      if(matchp1 > n){
        return true; 
      }
      
      //Search for the index of b[0] inside a
    
      if(a[(i + matchp1) % n] > b[(j + matchp1) % n]){
        i += matchp1;     
      } else{
        j += matchp1;
      }
    }
    return false;
  }
  
  //Khovanov’s Theta-Foam Evaluation
  public static int kTheta(int[] abc){

    //check if a + b + c = 3
    
    int sum = 0;
    
    for(int i = 0; i < abc.length; i++){
      sum += abc[i];
    }
    
    //return 0 if the sum != 3
    if(sum != 3){ 
      return 0;
    }
    
    //Check if the tuple is a cyclic permutation of (0,1,2) or (0, 2, 1), else it is 0

    int[] pos = {0, 1, 2};
    int[] neg = {0, 2, 1};

    if(isRot(abc, pos)){ 
      return 1;
    }else if(isRot(abc, neg)){
      return -1;
    }else{
      return 0;
    }
        
  }

  //For the output File
  public static List<String> linesOfFile = new ArrayList<>();
  public static StringBuilder output = new StringBuilder();
  
  //converts tuples to polynomials in the output string
  public static void dotstoPoly(int[] numTuple){
    int zeros = 0;
    for(int i = 0; i < numTuple.length; i++){
      if(numTuple[i] > 1){
        output.append(monomial[i] + numTuple[i] + "}");
      } else if(numTuple[i] == 1){
        output.append(monomial[i]);
        output.delete(output.length() - 2, output.length());
      } else{
        zeros++;  
      }
      
      if(zeros == numTuple.length){
        output.append("1");
      }
    }
    return;   
  }

  //the furnace evaluation. Returns t/f so it can be used in the base case of F3Images
  public static boolean furnaceEval(int[] dots){

    //Check if mlr is zero. If it is, then it returns true so F3Images skips to next combo. Otherwise it lets the calculator mode print the answer as zero.
    
    int[] mlr = new int[3];
    
    for(int i = 0; i < mlr.length; i++){
      mlr[i] = dots[i];
    }

    int[] zuv = {dots[DOTS - 2], 2, dots[DOTS - 1]}; //For theta evaluation
    
    //check if zuv or mlr is zero in each case in the evaluation formula

    int numZero = 0;
    int mlrSum = 0;
    int zeroCheck = 0;
    boolean isZero = false;

    for(int i = 0; i < zuv.length; i++){
      
      zuv[1] -= i;
      mlr[2] += i;
      
      mlrSum = 0;
      isZero = false;
      zeroCheck = 0;
        
      for(int j = 0; j < mlr.length; j++){
        mlrSum += mlr[j];
        if(mlr[j] > 2){
          isZero = true;  
          break;        
        } else if(mlr[j] < 1){
          zeroCheck++;
        }
      }

      if(isZero || (kTheta(zuv) == 0) || (mlrSum > 3) || (zeroCheck == 0)){
       numZero++;
      }

      zuv[1] += i;
      mlr[2] -= i;          
    }
      
    //if there are three zero terms then skip to next
     
    if(numZero > 2){
      output.append("0");
      return true;  //Lets F3Images skip
    }

    //write out answer in LaTeX/Reader friendly format
    
    int plusIndex = 0; //gets index of sign of first term when it is a plus     
    
    //Writes each term with sign
    
    for(int j = 0; j < 3; j++){
       
       zuv[1] -= j; //decrease degree of Z in each term 
       mlr[2] += j; //increase degree of Z in each term
        
      //Put the sign
      if(kTheta(zuv) < 0){
        output.append("+");
        plusIndex = output.length() - 1;
      } else if(kTheta(zuv) > 0){
        output.append("-");
      } else{
        //Reset zuv and mlr before skipping
        zuv[1] += j; 
        mlr[2] -= j;
        continue; //goes back and starts the next j. skips the current term
      }
     
      //Write the terms that appear
      dotstoPoly(mlr);
      
      //Reset the Z's for the next j
      
      zuv[1] += j; 
      mlr[2] -= j; 
    }
    
    //remove the leading plus sign
    if(output.charAt(plusIndex) == '+'){
      output.deleteCharAt(plusIndex);
    }
    
    return false;  
  }

  //prints an answer to the furnace evaluation in Calculator Mode
  public static void furnaceResult(int n, int m, int[] inputDots){
  
    //Writes result as user wishes
    boolean resultTape = false;

    if(m > 0){
      resultTape = true;
    }

    switch(n){

        case 0: //plug in polynomials
          output.append("F(");
          dotstoPoly(inputDots);
          break;

        case 1: //Original
          output.append("F_{" + inputDots[DOTS - 2] + inputDots[DOTS - 1] + "}(" + inputDots[0]);
          
          //write out coordinates in parentheses
          
          for(int i = 1; i < MLRLEN; i++){
            output.append(", " + inputDots[i]);
          }

          break;
                         
        default:
          output.append("F(");
          dotstoPoly(inputDots);
          break;                
      }
    
    if(resultTape){
      output.append(") &= ");
    } else{
      output.append(") = ");
    }

    furnaceEval(inputDots);
    System.out.println(output.toString());

    if(resultTape){
      output.append(" \\\\ "); // to help escape something weird
      //output.append(" & ");
      linesOfFile.add(output.toString()); //The Java 8 update changed append to add 
    }

    output.delete(0, output.length()); //clear the string for the next line
  }

  //Recursively generates all possible dot configurations and evaluates them 
  public static void F3Images(int n, int k){
    
    //base case for recursion

    if (k < 0){
                     
      //Change left column depending on user choice

      switch(n){
        case 0: //Original and Default
          output.append("F_{" + mlrab[DOTS - 2] + mlrab[DOTS - 1] + "}(" + mlrab[0]);
          
          //write out coordinates in parentheses
          
          for(int i = 1; i < MLRLEN; i++){
            output.append(", " + mlrab[i]);
          }
        
          //output.append(") = ");
          output.append(") &= "); //For alignment
          break;
        
        case 1: //plug in polynomials
          output.append("F(");
          dotstoPoly(mlrab);

          output.append(") &= "); //For alignment
          //output.append(") = ");
          break;
        
        case 2: //tabular form with polynomials
          dotstoPoly(mlrab);
          output.append(" & "); //For alignment in table or array 
          break;
                 
        default:
          output.append("F_{" + mlrab[DOTS - 2] + mlrab[DOTS - 1] + "}(" + mlrab[0]);
          
          //write out coordinates in parentheses
          
          for(int i = 1; i < MLRLEN; i++){
            output.append(", " + mlrab[i]);
          }
        
          //output.append(") = ");
          output.append(") &= "); //For alignment
          break;          
      }      
      
      if(furnaceEval(mlrab)){
        output.delete(0, output.length()); //clear the string for the next iteration
        return;
      } 
      
      //additional LaTeX formatting (optional)

      output.append(" \\\\ "); // to help escape something weird
      //output.append(" & ");   

      //generate the output in a list

      linesOfFile.add(output.toString()); //The Java 8 update changed append to add
      output.delete(0, output.length()); //clear the string for the next line
     
      return;
    }
    
    //Recursive generation and computation; generates all 5-tuples with 0, 1, or 2 only
    
    for(int i = 0; i < 3; i++){
      mlrab[k] += i;
      F3Images(n, k - 1);
      mlrab[k] -= i;
    }
    return;
    
  }

  public static void main(String[] args){
    Scanner scan = new Scanner(System.in);
    int n;

    System.out.print("Enter 0 for Calculator mode, 1 for a list of the nonzero evaluations: ");
    int c = scan.nextInt();
    
    if(c > 0){  
      System.out.print("Enter 0 for F_{ab}(mlr), 1 for F(X^mY^lZ^rU^aV^b), and 2 for tabular form with polynomial in H^{*}(J) on the left and image under Furnace on the right: ");
      n = scan.nextInt();

      F3Images(n, DOTS - 1);

      try{
        Files.write(Paths.get("HKJ.txt"), linesOfFile);
      } catch(Exception e) {
          e.printStackTrace();
      }
    } else{
          
      System.out.print("Enter 0 for no record of calculations, 1 for a text file record: ");   
      int m = scan.nextInt();
      
      System.out.println("In which form would you prefer your answer?");
      System.out.print("Enter 0 for F(X^mY^lZ^rU^aV^b), 1 for F_{ab}(mlr): ");
      n = scan.nextInt();
      
      while(true){
        int[] dots = new int[DOTS];
        System.out.println("Enter the number of dots (the power of):");
      
        for(int i = 0; i < monomial.length; i++){
          System.out.print(monomial[i]);
          dots[i] = scan.nextInt();
        }
       
        furnaceResult(n, m, dots);
      
        System.out.println("Would you like to do another calculation?");
        System.out.print("Enter 0 for yes 1 for no: ");     
        c = scan.nextInt();
       
        if(c > 0){
          break;
        }           
      }
       
      if(m > 0){
        try{
          Files.write(Paths.get("HKJRecord.txt"), linesOfFile);
        } catch(Exception e) {
            e.printStackTrace();
        }
      }
    }
  }
}
