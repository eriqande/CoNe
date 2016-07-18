#define UN_EXTERN



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"
#include "ECA_Opt3.h"
#include "cone.h"


int main(int argc, char **argv)
{
	int nf; /* number of founder lineages, assumed known here */
	double T; /* the number of generations */
	int i,NumReps;
	long seed1, seed2;
	char phrase[200];
	ldat *locarray;
	Nstruct *Ns;
	int NumLoc;
	int MaxN_index = 0;
	double MaxN_LogLike = -999999.9;
	double x,y,uconf,lconf,bd;
	double TotOverLociExtract;
	char FILENAME[10000];
	int max_nf=0,min_nf=99999;  /* for keeping track of where the probs will be zero over all loci */
	double all_loci_logl;
	
	gMax_n0 = 0;
	
	
	GetCoNeData(FILENAME,&locarray,argv,argc,&NumLoc,&T,&NumReps, phrase, &Ns);
	GetCoNeFileProbs(Ns, gNumNs, locarray,NumLoc,T); 

	
	EchoCoNeData(locarray,NumLoc);
		
	/* PrintNProbs(Ns, gNumNs, locarray, NumLoc, 1); */
/*	PrintNsLong(Ns);  */
	
	if(strcmp(phrase,"")!=0) {
		phrtsd(phrase,&seed1, &seed2);
		setall(seed1, seed2);
	}
	else {
		SeedFromFile("cone_seeds");
	}
	
	/* print out the run settings here */
	printf("SETTINGS : file-name : %s\n", FILENAME);
	printf("SETTINGS : pathway-to-probs-file : %s\n",gProbsPath);
	printf("SETTINGS : gens-between : %.4f\n", T);
	printf("SETTINGS : mc-reps : %d\n",NumReps);
	printf("SETTINGS : Ne\'s to evaluate likelihood at : ");
	for(i=0;i<gNumNs;i++)  printf("%.4f ",Ns[i].N);
	printf("\n");
	printf("SETTINGS : Scaled times to evaluate likelihood at : ");
	for(i=0;i<gNumNs;i++)  printf("%.4f ",Ns[i].t);
	printf("\n");
	printf("SETTINGS : prior : ",NumReps);
	if(gPrior==0.0) {
		printf("Using Uniform Prior for Allele Frequencies\n");
	}
	else {
		printf("Using Dirichet parameters %f/K for Allele Frequencies\n",gPrior);
	}
	
	if(strcmp(phrase,"")!=0) {
		printf("SETTINGS : seed-phrase : %s\n", phrase);
		printf("SETTINGS : Actual Ranlib Seeds :  %ld %ld",seed1,seed2);
	}
	else {
		printf("SETTINGS : No Seed Phrase Given\n");
	}
	printf("\n\n");
	

	for(i=0;i<NumLoc;i++)  {  double maxlogprob=-9999999.9;
		printf("LOCUS : %d : nf       Ave          Var               SD             NormoExtract      expNormoExtract     LogAvePlusNormoExtract  \n",i+1);
		if(locarray[i].Kin0 > max_nf) {  /* record the largest Kin0 of all the loci */
			max_nf = locarray[i].Kin0; 
		}
		if(locarray[i].n0 < min_nf) {  /* record the largest Kin0 of all the loci */
			min_nf = locarray[i].n0; 
		}
		for(nf=locarray[i].Kin0; nf<=locarray[i].n0; nf++) {
			AddReps(locarray[i].R,nf,NumReps, &(locarray[i]));
			printf("LOCUS : %d : %d  %e   %e     %e        %f        %e       %f\n",i+1,nf, locarray[i].R[nf]->Ave, locarray[i].R[nf]->Var, 
					sqrt(locarray[i].R[nf]->Var),locarray[i].NormoExtract[nf],
					exp(locarray[i].NormoExtract[nf]), log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf]);
			
			/* we are going to cut this off if the prob has dropped very low */
			if(  maxlogprob < log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf])
				maxlogprob = log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf];
			if(log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf] - maxlogprob < CALLITZERO) {
				printf("Truncating evaluation beyond this point.  Calling it zero...\n");
				if(nf<min_nf) {
					min_nf = nf;
				}
				while(nf<=locarray[i].n0) {
					locarray[i].R[nf]->Ave = 0.0;
					locarray[i].R[nf]->Var = 0.0;
					locarray[i].NormoExtract[nf] = -999999.9;
					nf++;
				}
				
			} 
			

		}
	}
	
	/* now, here we are going to print out the log-likelihood for each nf.  I don't want
	to deal with calculating the MC variance for them at this point, because that is a bit hairy 
	with all the renormalizing stuff.  I'll do that later. */

/**********  I 	CUT THIS SECTION OUT BECAUSE IT DIDN'T MAKE SENSE IF THE NUMBER
 OF LINEAGES AT TIME T IS NOT THE SAME FOR ALL LOCI **********/
 /*
	printf("ALL_LOCI :    nf     LogL\n");
	for(nf=max_nf;nf<=min_nf;nf++) { 
		all_loci_logl = 0.0;
		for(i=0;i<NumLoc;i++)  {
			all_loci_logl += log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf];
		}
		printf("ALL_LOCI :    %d     %f\n",nf,all_loci_logl);
	}
* END OF REMOVED SECTION. */
	
	TotOverLociExtract = ComputeNeLikes(locarray,Ns,NumLoc);
	
	printf("NE_LOGLIKE  :    Ne      t      LogLike    MC.CI.Lo  MC.CI.Hi      Like       SD   \n");
	for(i=0;i<gNumNs;i++)  {
		printf("NE_LOGLIKE  :  %.3f  %f  %f  %f  %f   %e     %e\n",Ns[i].N,Ns[i].t, 
						log(Ns[i].Like),
						log(Ns[i].Like-1.96*sqrt(Ns[i].LikeVar) ), log(Ns[i].Like+1.96*sqrt(Ns[i].LikeVar) ),
						Ns[i].Like,sqrt(Ns[i].LikeVar));
		if(MaxN_LogLike < log(Ns[i].Like) ) {
			MaxN_LogLike = log(Ns[i].Like);
			MaxN_index = i;
		}
	}
	
	
	/* now report the max (by parabolic interpolation) and the 1.96 units of 
		log likelihood "confidence intervals" (by linear interpolation).  */
	if(MaxN_index == 0  || MaxN_index == gNumNs-1) {
		printf("MaxByParabolicInterplotation :  -999.999  -999.999  -999.999\n");
	} 
	else {
		i=MaxN_index;
		y = ParabolicMaxInterp(Ns[i-1].N,Ns[i].N,Ns[i+1].N,
				log(Ns[i-1].Like),log(Ns[i].Like),log(Ns[i+1].Like), &x);
		printf("MaxByParabolicInterpolation :  %f   %f    %f  :  Ne  t  LogLike \n",x,(double)T/(2.0*x), y);
	}
	lconf = -999.999;
	uconf = -999.999;
	bd = y - 1.96;
	for(i=0;i<gNumNs-1;i++)  {
		if(log(Ns[i].Like) < bd && log(Ns[i+1].Like) > bd)  {
			lconf =  Ns[i].N + (bd - log(Ns[i].Like)) * 
						(Ns[i+1].N - Ns[i].N) /  ( log(Ns[i+1].Like) - log(Ns[i].Like) );
		}
		if(log(Ns[i].Like) > bd && log(Ns[i+1].Like) < bd)  {
			uconf = Ns[i].N + (bd - log(Ns[i].Like)) * 
						(Ns[i+1].N - Ns[i].N) /  ( log(Ns[i+1].Like) - log(Ns[i].Like) );
		}
	}
	printf("LowerSupportLimit  : %f   %f  :  N    t  \n",lconf, T/(2.0*lconf) );
	printf("UpperSupportLimit  : %f   %f  :  N    t \n",uconf, T/(2.0*uconf) );
	
	printf("ThisAmountTakenOutOfLogLikelihood :  %f\n", TotOverLociExtract);
	
	SeedToFile("cone_seeds");
	return(0);
}