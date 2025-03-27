/***************************************************************
 *                  Simple FCM Cluster                         *
 * *************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MaxIter 20

FILE *ReadFilePtr(FILE *f) {
  f =fopen("Data.txt","r");
  if(f == NULL) {
    printf("File name doesn't exist");
    exit(1);
  }
    
  return f;
}
void PrintDouArr(double **D, int NoDtpts, int NoDim) {
  int i, j;
  
  for(i = 0; i < NoDtpts; i++) {
    for(j = 0; j < NoDim; j++) {  
        printf("%6.4lf  ",D[i][j]);
    }
    printf("\n");
  }
}
double **initializeDouPtr(double **D, int NoDtpts, int NoDim) {
  int i, j;
  
  for(i = 0; i < NoDtpts; i++) {
    for(j = 0; j < NoDim; j++) {  
        D[i][j] = 0;
    }
  }
  return D;    
}
double **CreateDouPtr(double **D, int NoDtpts, int NoDim) {
  int i; 
  D = (double**)malloc(NoDtpts * sizeof(double*));
  if(D == NULL) {
    printf("\nError in Memory Allocation\n");
    exit(1);
  } else {
      for(i = 0; i < NoDtpts; i++) {
        D[i] = (double*)malloc(NoDim * sizeof(double));
        if(D[i] == NULL) {
          printf("\nError in Memory Allocation\n");
          exit(1);
        }
      }
  }
  
  D = initializeDouPtr(D, NoDtpts, NoDim);
  
  return D;
}
double **ReadData(double ** D, FILE* fp, int NoDtpts, int NoDim) {
  int i, j;
  for(i = 0; i < NoDtpts; i++) {  
    for(j = 0; j < NoDim; j++) {  
        fscanf(fp, "%lf  ", &D[i][j]);
    }
  }
  return D;    
}
double **RandonGenDegofMem(double ** D, int NoDtPts, int NoClust) {
  int i, j;
    
  for(i = 0; i < NoDtPts; i++) {
    for(j = 0; j < NoClust; j++) {
       D[i][j] = (double)rand()/RAND_MAX; 
    }
  }
  return D;
}
double **ModifiedDegofMem(double** D, int NoDtPts,int NoClust) {
  double row_sum, den;
  int i, j;
  
  for(i = 0; i < NoDtPts; i++) {
    row_sum = 0;
    for(j = 0; j < NoClust; j++) {  
      row_sum += D[i][j];   // Probability sum
    }
    // printf("\nrow_sum =%lf", row_sum);
    
    for(j = 0; j < NoClust; j++) {  
      den = D[i][j]/row_sum;   // Probability sum
      D[i][j] = den;
      //printf("\t\tden =%lf", den);
    }
  }
  return D;  
}
double Power(double res, double a, double b) {
  int n1, n2;
     
  res = pow(a,b);
  
  return res;  
}
double **ComputePowerofDegofMem(double **CV, int NoDtPts, int NoClust, double **DegofMem, double fuzziness) {
  int i, j;
  double p1;
  
  for(i = 0; i < NoDtPts; i++) {
    for(j = 0; j < NoClust; j++) {
      CV[i][j] = Power(CV[i][j], DegofMem[i][j], fuzziness);
    }
}
  return CV;  
}
double **CalculateCentVect(double **CC, int NoDtPts, int NoClust, int NoDim, double **DegofMem, double fuzziness, double **D, double **CV) {
  int i, j, k;
  double num, den, value, value1;

  for(j = 0; j < NoClust; j++) {
     for(k = 0; k < NoDim; k++) {
       num = 0.0;
       den = 0.0;
       for(i = 0; i < NoDtPts; i++) {
         // printf("\nD[%d][%d] = %lf", i, k, D[i][k]);
         value = 0.0; value1 = 0.0;
         value = CV[i][j]*D[i][k];
         value1 = CV[i][j];
         // printf("\t\tValue = %lf\tvalue1 = %lf", value, value1);
         num += value;
         den += value1;
        }
       CC[j][k] = num/den;
    }
  }  
  return CC;  
}

/**************** Inside function **********************/
double GetNorm(double res, int i, int j, double **D, double **CC, double fuzziness, int NoDim) {
  int k;
  double sum, term1, term2;
  sum = 0.0; term1 = 0; term2 = 0;
  for(k = 0; k < NoDim; k++) {
      term2 = D[i][k]-CC[j][k];
      term1 = pow(term2,fuzziness);
      sum = sum + term1;
  }
  res = sqrt(sum);
  return res;
}
double GetNewValue(double Newuij, int i, int j, double **D, double **CC, int NoDim, int NoClust, double fuzziness) {
  int k;
  double t, p, sum, num, den, p1;
  sum = 0.0;
  p1 = fuzziness -1;
  p = 2/p1;
  for(k = 0; k < NoClust; k++) {
    num = 0; den =0;
    num = GetNorm(num, i, j, D, CC, fuzziness, NoDim);
    den = GetNorm(den, i, k, D, CC, fuzziness, NoDim);
    t = num/den;
    t = pow(t,p);
    sum += t;
  }
  Newuij = 1.0/sum;
  return Newuij;
}
double UpdateDegreeofMembership(double MaxDiff, int NoClust, int NoDtPts, int NoDim, double fuzziness, double **Data, double **DegofMem, double **CentreV, double **CC) {
  int i, j;
  double Newuij;
  double diff, maxdiff, a;
  maxdiff = 0;
  
  for(j = 0; j < NoClust; j ++ ) {
     for(i = 0; i < NoDtPts; i ++ ) {
        Newuij = GetNewValue(Newuij,  i, j, Data, CC, NoDim, NoClust, fuzziness);
        diff = 0;
        diff = Newuij-DegofMem[i][j];
        // printf("\nNewuij = %lf\tDegofMem[%d][%d] = %lf\t diff = %lf", Newuij, i, j, DegofMem[i][j], diff);
        if(diff > maxdiff)
            maxdiff = diff;
        DegofMem[i][j] = Newuij;
     }
    }
    MaxDiff = maxdiff;
  return MaxDiff;
}


int main(int argc, char *argv[]) {
  
  int i;
  char fname[20];
  
  int NoDtpts, NoDim, NoClus, fuzziness, epsilon;
  
  NoDtpts = 5;
  NoDim = 2;
  
  NoClus = 2; fuzziness = 2.0; epsilon = 0.01;
    
  // ***** Pointer Creation ***** 
  
  FILE *fp = ReadFilePtr(fp);
  
  double **Dat = CreateDouPtr(Dat, NoDtpts, NoDim);
  // PrintDouArr(Dat, NoDtpts, NoDim);
  
  double ** DegMem = CreateDouPtr(DegMem, NoDtpts, NoClus);
  // PrintDouArr(Dat, NoDtpts, NoClus);
  
  double **CentreV = CreateDouPtr(CentreV, NoDtpts, NoClus); // Calculation of centre vector 
  // PrintDoublePtr(CentreV, NoDtpts, NoDim);
  
  double **CC = CreateDouPtr(CC, NoClus, NoDim);
  // PrintDouArr(Dat, NoClus, NoDim);
  
  // ***** Functions ***** 
  
  Dat = ReadData(Dat, fp, NoDtpts, NoDim);
  printf("\nInput Data\n\n"); PrintDouArr(Dat, NoDtpts, NoDim);
  
  DegMem = RandonGenDegofMem(DegMem, NoDtpts, NoClus);
   printf("\nActual generated random number\n"); PrintDouArr(DegMem, NoDtpts, NoClus);
  
  DegMem = ModifiedDegofMem(DegMem, NoDtpts, NoClus);
   printf("\nModified random number\n"); PrintDouArr(DegMem, NoDtpts, NoClus);
  
  double MaxDiff;
  int iter;
  iter = 0;
  MaxDiff = 0;
  do {
    CentreV = initializeDouPtr(CentreV, NoDtpts, NoClus);
    
    CentreV = ComputePowerofDegofMem(CentreV, NoDtpts, NoClus, DegMem, fuzziness);  
    // printf("\nPower of DegofMem\n\n"); PrintDouArr(CentreV, NoDtpts, NoClus);
    
    CC = initializeDouPtr(CC, NoClus, NoDim);
    
    CC = CalculateCentVect(CC, NoDtpts, NoClus, NoDim, DegMem, fuzziness, Dat, CentreV);
    // printf("\nCluster Center is as follows:\n\n");  PrintDoublePtr(CC, NoClus, NoDim);
    
    MaxDiff = UpdateDegreeofMembership(MaxDiff, NoClus, NoDtpts, NoDim, fuzziness, Dat, DegMem, CentreV, CC);
    printf("\nMaxDiff = %lf", MaxDiff);
    
    printf("\nUpdated DegMem @ iter = %d\n", iter); PrintDouArr(DegMem, NoDtpts, NoClus);
    
    if(MaxDiff <= epsilon)
        break;
    
    iter = iter + 1;
  } while(iter < 20);
  
  printf("\n\tNumer of iteration = %d\n\n", iter);
  
  

  
  
  // Free memory
  fclose(fp);
  for(i = 0; i < NoDtpts; i++) {
      free(Dat[i]); free(DegMem[i]); free(CentreV[i]);
  } free(Dat); free(DegMem); free(CentreV);
  for(i = 0; i < NoDim; i++) {
      free(CC[i]); } free(CC);
  return 0;  
}
