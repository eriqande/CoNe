#define UN_EXTERN

#define MAX_ALLELES 500

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"




int main(int argc, char **argv)
{
	int n0,i,j,num_ts,k,test,n,reps=1000000;
	double *ts,t,y;  /* a list of the ts that we want to compute probs for */
	double mean;
	ival ***r;  /* to store all the average values. The r.v.'s themselves are 0/1 indicators
					subscripted by [ time ] [ i ] where i is the number of lineages */
	double *ctimes;  /* coalescent times= ctimes[i] is the time at which the genealogy
						having i lineages became one having i-1 */
	
	
	
	if(argc==1) {
		printf("usage: simCoNeprob n0  reps\n\n");
		return(0);
	}
	
	n0 = atoi(argv[1]);
	reps = atoi(argv[2]);
	num_ts = 100;
	
	ctimes = (double *)ECA_CALLOC(n0+5,sizeof(double));
	ts = (double *)ECA_CALLOC(num_ts,sizeof(double));
	r = (ival ***)ECA_CALLOC(num_ts,sizeof(ival **));
	for(j=0;j<num_ts;j++)  {
		r[j] = IvalVector(0,n0,0,0,0);
	}
	
	/* here I set the values of t that I want */
	y = 1000.0;
	ts[0] = 1/y;
	for(i=1;i<50;i++)  {
		y-=18.36734694;
		ts[i] = 1/y; 
	}
	y=100.0;
	ts[49]=1/y;
	for(i=50;i<100;i++)  {
		y-=1.66666666666;
		ts[i] = 1/y; 
	}

/*	for(i=0;i<num_ts;i++)  {
		printf("%d  %f\n",i,ts[i]);
	}
*/

	
	for(n=0;n<reps;n++) {
		/* simulate the coalescent times */
		for(i=n0,t=0.0;i>1;i--)  {
			mean = 2.0 / ((double)i*(i-1.0));
			t += (double)genexp(mean);
			ctimes[i] = t;
		}
		
		/* then record those */
		/* initialize everything to a 0 */
		for(j=0;j<num_ts;j++)  {
			for(k=0;k<=n0;k++)  {
				r[j][k]->v = 0;
			}
		}
		test = n0;
		j=0;
		while(test > 1 && j < num_ts)  {
			/* printf("test= %d   t= %f  j= %d\n", test,ctimes[test],j); */
			if(test==n0 && ts[j] < ctimes[test]) {
				r[j][test]->v = 1;
				/*printf("\nts[%d]= %f   Has: %d\n",j,ts[j],test);*/
				j++;
			}
			else if(ts[j] < ctimes[test] && ts[j] > ctimes[test+1]) {
				r[j][test]->v = 1;
				/*printf("ts[%d]= %f   Has: %d\n",j,ts[j],test+1);*/
				j++;
			}
			else if(ts[j] > ctimes[test]) {
				test--;
			}
		}
		/*printf("\n\nCatching the extra bits.  test= %d  j= %d\n\n",test,j); */
		while(j<num_ts) {  /* then take care of the remaining time points, if any */
			if(ts[j] < ctimes[2])  {
				r[j][test+1]->v = 1;
				/*printf("ts[%d]= %f   Has: %d\n",j,ts[j],test+1); */
			}
			else {
				r[j][test]->v = 1;
				/* printf("ts[%d]= %f   Has: %d\n",j,ts[j],test);*/
			}
			j++;
		}
		/* now increment all the ivals */
		for(j=0;j<num_ts;j++)  {
			for(k=0;k<=n0;k++)  {
				IncrementIval(r[j][k]);
			}
		}
		
		/*  printf("Completed Rep %d\n",n+1);  */
	}

	/*  then print all the results to standard out */
	for(j=num_ts-1;j>=0;j--)  {int start=0; int end=n0; int cnt=0; int alreadyprinted=0; double psum=0.0;
		printf("t: %.12f ",ts[j]);
		for(k=1;k<=n0;k++)  {
			if(r[j][k]->Ave > 0.0 && start==0) {
				start = k;
				printf("start: %d ",k);
			}
			if(start>0) {
				if(r[j][k]->Ave == 0.0) {
					cnt++;
				}
				psum += r[j][k]->Ave;
			}
			if(cnt>10 && psum > 1.0 - (1.0/reps) ) {  /* count it as an end if all (or about all) of  the probability mass is there
													    and there have been a lot of zeroes already. */
				end = k;
				printf("end: %d ",end);
				alreadyprinted=1;
				break;
			}
		}
		if(!alreadyprinted)
			printf("end: %d ",end);
		for(k=start;k<=end;k++)  {
			printf("%f ",r[j][k]->Ave);
		}
		printf("\n");
	}
	
}