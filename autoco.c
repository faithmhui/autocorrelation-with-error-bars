// autocorrelation-with-error-bars \(o^O)/ 
// Runs autocorrelation on 1-D data and also returns the errors on each time lag.
// 2015 Faith-Mei-Hui 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#define MAX 1000   //max number of data

void autoco (int nbin, double bincenter[],double data [],double sdata[], FILE ** pntr){
  
  //find the mean;
  double sum=0;
  for (int i=0; i<nbin; i=i+1){
    sum=sum+data[i]*bincenter[i];
  }
  double mean= sum/nbin;
  
  //find the variance
  sum=0;
  for(int k=0;k<nbin;k++){
    sum=sum+(data[k]-mean)*(data[k]-mean);    
  }
  double sumdsqr=sum;
  printf("mean=%f sumdsqr =%f\n",mean, sumdsqr);
  
  
  /*autocorelation starts here*/
  fprintf(*pntr,"#tau, position of bin,autosumresult, autosumresulterr\n");
  for (int tau=0;tau<nbin;tau++){ 
    double autosum=0; 
    for (int bincount=0;bincount<nbin;bincount++){

      /*autocorelation equation*/
      autosum=(data[bincount-tau]-mean)*(data[bincount]-mean)+autosum;
    }//for
    
    /*autocorrelation error analysis*/
    double autosumerr=0;    
    for (int bincount=0;bincount<nbin;bincount++){
 

      if(bincount+tau>nbin)            //when autocorrelation exceed the first histogram set to 0
        data[bincount+tau]=0;
      //differentiate with tau and differentiate with bincount
      double numer1=(data[bincount-tau]-mean + data[bincount+tau]-mean);   
      double numer2=2.0*(data[bincount]- mean);  
      double denom=sumdsqr*sumdsqr;
      double combined=(numer1*sumdsqr-numer2*autosum)/denom * sdata[bincount];
      autosumerr=combined*combined+autosumerr;
  
    }//for
    double autosumresult=autosum/sumdsqr;
    double autosumresulterr=sqrt(autosumerr);//to make positive the error
   
    fprintf(*pntr,"%d %f %f %f\n",tau,bincenter[tau],autosumresult,autosumresulterr);//print autosumresult and its error
  
  }//for
  
}
//=======================================================================================
int main(){

  FILE *inp = fopen("auto.txt","r");           // retrieve three column spaced data eg. 1 2 3
  FILE *outp = fopen("resutlauto.txt","w");    // print result to here
  if (inp == NULL || outp == NULL)
    {  
      printf("Error! Could not open file\n");
      exit(-1);
    }

  int i=0;
  double data[MAX]={0};
  double bincenter[MAX]={0};
  double sdata[MAX]={0};
  double n=0;
  double m=0;
  double s=0;

  
  while( fscanf(inp, "%lf %lf %lf", &n, &m, &s) !=EOF){    //Unpack data into array
    bincenter[i]=n;
    data[i]=m;
    sdata[i]=sqrt(s);
    i=i+1;
  }
  int total=i;
  
  autoco (total, bincenter ,data ,sdata, &outp );


}




