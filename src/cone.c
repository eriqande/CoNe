#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"
#include "ECA_Opt3.h"
#include "cone.h"

	
/*	
void PrintUsage(char *msg) 
{
	printf("\n%s\n",msg);
	printf("Usage:\n\nCoNe pathway_to_probs_files -c  K    AT,1 ... AT,K    A0,1 A0,2 ... A0,K     T  NumReps  SeedPhrase Nlo Nhi Nstep  [Prior] \n\nor");
	printf("\n\nCoNe  pathway_to_probs_files -f filename T  NumReps  SeedPhrase Nlo Nhi Nstep [Prior] \n\n");
	
	printf("The likelihood will be computed for values of Ne:\n");
	printf("	Nlo, Nlo+step, Nlo+2*step,...NloTop.  NloTop is <= Nhi \n");
	printf("\npathway_to_probs_files is the path to the directory that holds\n");
	printf("the files that have all the small-t probabilities in them.\n");
	printf("For example, on my system it is:\n\t~/Documents/eca_code/CoNe/data/\n\n");
	printf("Note that you must have the final slash on the end.\n\n");
	printf("\n\n[Prior] represents an optional argument specifying \nthe prior distribution for the allele frequencies:\n");
	printf("\t0  --   The prior is a uniform Dirichlet distribution");
	printf("\n\tR, where R is a real number greater than 0.0:  \n\t\tThe parameters of the Dirichlet prior are R/K.");
	printf("\n\t\tHence, R=1 is the \"Unit information prior\".  This is the default.\n\n");
	exit(1);
}
*/	



void GetCoNeData(char *FILENAME, ldat **loc, char **argv, int argc, int *NumLoc, double *T,  int *NumReps, char *phrase, Nstruct **Ns)
{
	int i,j,l,K,n0,nT,temp;
	FILE *in;
	double Nlo,Nhi,Nstep;
	int NumSam;

	int File_f = 0,
		T_f = 0,
		Path_f = 0,
		N_f = 0,
		Seed_f = 0,
		Nlo_hi_step_f = 0,
		Prior_f = 0;
		
	
	DECLARE_ECA_OPT_VARS;
	
	/* some defaults */
	gPrior = 1.0;
	sprintf(phrase,"");
	
	
	SET_OPT_WIDTH(28);
	SET_ARG_WIDTH(17);
	SET_PROGRAM_NAME("cone");
	SET_PROGRAM_SHORT_DESCRIPTION("a program for estimating Ne");
	SET_PROGRAM_LONG_DESCRIPTION(
		CoNe
		computes the likelihood of Ne given data on two temporally spaced
		genetic samples.  The statistical model used is based on the
		coalescent of the gene copies drawn in the second sample\054 as
		described in Berthier et al. (2003) Genetics 2003. 160:741-51.
		The Monte Carlo computations to compute the likelihood\054
		however\054 were developed by Eric Anderson\054 and are orders of
		magnitude faster than previous implementations. 
 
		\n\nDetails of the algorithm are given in Anderson (2005) Genetics 170:955-967.
	);
	SET_VERSION("VERSION: 1.02\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 1 August 2007")
	SET_PROGRAM_AUTHOR_STRING("Eric C. Anderson (eric.anderson@noaa.gov)");
	SET_VERSION_HISTORY("\
				\n\nVERSION: 1.02\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 1 August 2007\
				\nCHANGES:\
				\n1. Modified options input to use ECA_Opt3, so it can ouput to guiLiner.\
				\n\nVERSION: 1.01\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 28 September 2005\
				\nCHANGES:\
				\n1. Added the LOC_SPECIFIC_LOGLS lines on the output.  These\
				\n   provide locus specific log-likelihood curves.\
				\n2. Removed the ALL_LOCI lines of output because they are meaningless\
				\n   when sample sizes differ between loci.\
				\n3. Fixed a bug in the companion utility \"simCoNeprob\" that \
				\n   caused spurious results if the number\
				\n   of lineages at time T was very large, and few replicates were\
				\n   run in simCoNeprob.  (Thanks to Stuart Barker for catching this.) \
				\nCOPYRIGHT: Federal Gov't Work.  No copyright.\n\n\
				\n\nVERSION: 1.0\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 7 March 2005\nCOPYRIGHT: Federal Gov't Work.  No copyright.\n\n\
				\n\nVERSION: 1.0 beta\nAUTHOR: Eric C. Anderson (eric.anderson@noaa.gov)\nDATE: 21 April 2004\nCOPYRIGHT: Federal Gov't Work.  No copyright. \n\n")
	

	BEGIN_OPT_LOOP
	
		OPEN_SUBSET(Data Analysis Options, Data Analysis Options, These options control the inputs and methods for the CoNe analysis);
	
		if ( REQUIRED_OPTION(Data File,
				File_f,
				f,
				file-name,
				F,
				pathname of the data file,
				F is the name of the file in which you have your data.  It is in the same format as
				data files for TM3 (by Pierre Berthier) 
				and TMVP (by Mark Beaumont)
				The data file should start  with a 0 (this is a strange vestige of some sort 
				from TM3 or TMVP) followed by the number of time periods
				***which in this case must
				always be 2*** 
				followed by the number of loci.  Then data for each locus consists of
				the number of alleles observed at the locus followed by a row of counts
				of the different alleles observed in the first sample and a row of counts
				of the different alleles observed at the second sample. (By first sample I mean the
				sample taken first in time going forward.  Hence\054 the second sample is the sample that
				was collected most recently.)  There must be only
				integers and whitespace in the file.  An example file is shown in the FILES section
				of the manual pages.  NOTE! It turns out that it is not essential that the counts of alleles
				observed in the first sample (the one further back in the past) be integers.  It turns out
				to be convenient to express them as real numbers in some cases\054 so I have recoded it so that they
				can be real numbers and not just integers.  Now when the data file is echoed to standard input\054 these counts are expressed
				as real numbers.  Do not let this alarm you.  Everything is OK\054 still.
				 ) ) {
			if( ARGS_EQ(1) ) {
				GET_STR(FILENAME);
			}
		}
		if ( REQUIRED_OPTION(Probs File Path,
					Path_f,
					p,
					path-to-probs-files,
					D, 
					directory path to XXXpr.txt files ,
					This is the pathway to files containing the precomputed probabilities of
					having j lineages remaining at scaled time t given that you started with
					i lineages.
					Note that the trailing slash is required.
					For example on my system 
					~/Documents/eca_code/CoNe/probs/ is the pathway.
					The probs files are a collection of files named
					XXXpr.txt
					where XXX is a number giving the number of gene copies in the 
					second sample.  These files have been precomputed using the program
					simCoNeprob and are described below in FILES.  The CoNe distribution
					includes precomputed XXXpr.txt files with XXX ranging from 10 to 400.  This
					represents samples of between 5 and 200 diploid organisms.  For different sample
					sizes it is necessary to create new XXXpr.txt files using simCoNeprob which is
					also included with the CoNe distribution.
		) ) {
			if( ARGS_EQ(1) ) {
				GET_STR(gProbsPath);
			}
		}
		if ( REQUIRED_OPTION(
				Number of Generations,
				T_f,
				T,
				gens-between,
				R,
				generations between samples, 
				R is the number
					of generations between samples.  It may be specified as a non-integer in order to allow
					for a non-integer number of generations.  
		) ) {
			if( ARGS_EQ(1) ) {
				*T = GET_DUB;
			}
		}
		if ( REQUIRED_OPTION(
				Monte Carlo Reps,
				N_f,
				m,
				mc-reps,
				J,
				number of Monte Carlo reps,
				J is the number
					of importance sampling reps to perform for each value of the number of 
					genes ancestral to the second sample on each locus.  I have found the importance
					sampling algorithm to be good enough that 100 J=100 gives reliable results and
					usually runs very quickly.  However; on a final run of your data set J should be
					much larger.  You can get an idea of whether J should be larger by the width of the
					Monte Carlo confidence intervals around the estimated likelihood curve.
		) ) {
			if( ARGS_EQ(1) ) {
				*NumReps = GET_INT;
			}
		}
		if ( REQUIRED_OPTION(
				Ne Values,
				Nlo_hi_step_f,
				n,
				ne-lo-hi-step,
				R1 R2 R3, 
				values of Ne to compute L(Ne) at ,
				sets the values of Ne for which the likelihood will be computed.
					R1 is the lowest value of Ne. R2 is the highest value of Ne. R3 is
					the step size between values of Ne.  For values of Ne such that T/(2Ne)
					is smaller than .06 (or so) the precomputed values of scaled time
					(stored in the appropriate XXXpr.txt file (see the description of the -p option)) will be
					used.  So the step size may not be that given by R3.  Arguments need 
					not be given as real numbers.  An integer (like 250) will work just fine.
					) ) {
			if( ARGS_EQ(3) ) {
				Nlo = GET_DUB;
				Nhi = GET_DUB;
				Nstep = GET_DUB;
				if(Nlo>Nhi) {
					fprintf(stderr,"Error Processing Option -n/--ne-lo-hi-step! The first argument is less than the second argument.\n");
					OPT_ERROR;
				}
			}
		}
		if (OPTION(
				Allele Freq Prior,
				Prior_f,
				q,
				prior,
				R,
				allele frequency prior parameter,
				R specifies the prior distribution  for  the
              allele frequencies.  R=0 makes the prior a uniform Dirichlet distribution.
              R>0 makes the parameters
              of the Dirichlet prior are R/K. Hence R=1 is  the so-called unit 
              information prior.  R=1 is the default.) ) {
			if( ARGS_EQ(1) ) {
				gPrior = GET_DUB;
				if(gPrior < 0.0) {
					fprintf(stderr,"Error Processing Option -q/--prior!  The argument of this option must be greater than zero\n");
					OPT_ERROR;
				}
			}
		}
		if ( OPTION(
				Random Seed Phrase,
				Seed_f,
				s,
				seed-phrase,
				S,
				random number seed-phrase ,
				S is a single string (no spaces) that will be used to seed the random number generator.  If 
			this option is not invoked then a seed is chosen based on the current time or---if the file
			cone_seeds is present---the seeds are taken from that file.  Upon completion of the 
			program the next random number seeds in series are printed to the file cone_seeds.
			) ) {
			if( ARGS_EQ(1) ) {
				 GET_STR(phrase);;
			}
		}
		CLOSE_SUBSET
		
	END_OPT_LOOP
		
/*		
	if(!Prior_f)  {  /* in this case, there is not  specification for PRIOR */
/*		printf("Using Default Unit-Information Prior for Allele Frequencies\n");
	}
	else if(gPrior==0.0) {
		printf("Using Uniform Prior for Allele Frequencies\n");
	}
	else {
		printf("Using Dirichet parameters %f/K for Allele Frequencies\n",gPrior);
	}
*/
	
	if( (in = fopen(FILENAME,"r")) == NULL ) {
		printf("\nCan't open file: %s\n\nExiting...\n",FILENAME);
		exit(1);
	}
	fscanf(in,"%d",&temp);
	if(temp != 0)  {
		printf("\nSorry, can't take transposed data.  Data file must start with a 0\n\nExiting...\n\n");
		exit(1);
	}
	fscanf(in,"%d",&NumSam);
	if(NumSam != 2)  {
		printf("\nSorry, this program allows for only two samples in time.\n\nExiting...\n\n");
		exit(1);
	}
	fscanf(in,"%d",NumLoc);
	*loc = (ldat *)ECA_CALLOC(*NumLoc, sizeof(ldat));
	for(l=0;l<*NumLoc;l++) {
		fscanf(in,"%d",&((*loc)[l].K) );
		K = (*loc)[l].K;
		(*loc)[l].A0 = (int *)ECA_CALLOC(K,sizeof(int));
		(*loc)[l].AT = (double *)ECA_CALLOC(K,sizeof(double));
		(*loc)[l].predpars = (double *)ECA_CALLOC(K,sizeof(double));
		for(j=0,nT=0;j<K;j++)  {
			fscanf(in,"%lf",&((*loc)[l].AT[j]));
			nT += (*loc)[l].AT[j];
			(*loc)[l].predpars[j] = (*loc)[l].AT[j] + 1.0;
		}	
		(*loc)[l].Kin0 = 0;
		for(n0=0,j=0;j<K;j++) {
			fscanf(in,"%d",&((*loc)[l].A0[j]) );
			n0 += (*loc)[l].A0[j];
			(*loc)[l].Kin0 += ((*loc)[l].A0[j] > 0);
		}
		(*loc)[l].P = NULL;
		(*loc)[l].R = DvalVector(0,n0,0.0,0.0,0.0);
		(*loc)[l].n0 = n0;
		(*loc)[l].nT = nT;
		(*loc)[l].NormoExtract = (double *)ECA_CALLOC(n0+1,sizeof(double));
		if(gMax_n0 < n0) gMax_n0 = n0;
		
	}
	fclose(in);
	
	/* Now, at the end here we can deal with allocating memory to, and setting up,
		the Ns struct where the likelihoods will be stored */
	*Ns = (Nstruct *)ECA_CALLOC( (int)((Nhi-Nlo)/Nstep) + FILE_LINES, sizeof(Nstruct));
		
	for(i=0;Nlo<=Nhi;Nlo+=Nstep)  {double t; int nn;
		t = (double)(*T)/(2.0*Nlo);
		if(t >= .06) {  /* if we can actually compute this thing we do, otherwise we rely on the probs in files */
			(*Ns)[i].N = Nlo;
			(*Ns)[i].t = t;
			(*Ns)[i].LocLikes = (double *)ECA_CALLOC(*NumLoc,sizeof(double));
			(*Ns)[i].LocLikeVars = (double *)ECA_CALLOC(*NumLoc,sizeof(double));
			(*Ns)[i].CP = (CoalProbs *)ECA_CALLOC(gMax_n0 + 1, sizeof(CoalProbs));
			for(j=0;j<=gMax_n0;j++)  {
				(*Ns)[i].CP[j].P = NULL;
			}
			
			/* here we allocate memory to the CoalProbs[n0].P (if necessary) and then also compute them  */
			for(l=0;l<*NumLoc;l++)  { double y;  int BackToZero; int GotBeyondZero;
				nn = (*loc)[l].n0;  
				if((*Ns)[i].CP[nn].P == NULL)  {  /* in this case we have to allocate memory and do the computations */
					(*Ns)[i].CP[nn].P = (double *)ECA_CALLOC(nn+1,sizeof(double));
					
					/* cycle over the values of j and record the probs */
					BackToZero = 0;
					GotBeyondZero = 0;
					for(j=1;j<=nn;j++)  {
						if(BackToZero==0)  {
							y = Tavare1984by_recurs(nn,j,t);
							if(y < 1.0e-10)
								y = 0.0; 
						}
						if(y > .1 / (double)i )
							GotBeyondZero = 1;
						if(GotBeyondZero==1 && y == 0.0) {
							BackToZero=1;
							y = 0.0;
						}
						(*Ns)[i].CP[nn].P[j] = y;
						
					}
					/* then cycle again to record the start and end */
					(*Ns)[i].CP[nn].start = 0;
					for(j=1;j<=nn;j++)  {
						if((*Ns)[i].CP[nn].P[j] > 0.0) {
							(*Ns)[i].CP[nn].start = j;
							break;
						}
					}
					(*Ns)[i].CP[nn].end = nn;
					for(j=nn;j>=1;j--)  {
						if((*Ns)[i].CP[nn].P[j] > 0.0) {
							(*Ns)[i].CP[nn].end = j;
							break;
						}
					}
				}  /* ends the if(...P==NULL) conditional */
			
			} /* ends the for loop over l */
			i++;  
		}  /* ends the if(t >= .06) */
	}
	gNumNs = i;
		
}

void PrintNsLong(Nstruct *Ns)
{
	int i,l,j;
	
	for(l=0;l<=gMax_n0;l++)  {
		for(i=0;i<gNumNs;i++)  {
			if(Ns[i].CP[l].P != NULL) {
				for(j=Ns[i].CP[l].start; j<=Ns[i].CP[l].end;j++)  
					printf("n0= %d  j= %d  prob= %e  t= %f  N = %f\n",l,j,Ns[i].CP[l].P[j],Ns[i].t,Ns[i].N);
			}
		}
	}

}

/* for getting the file probs *after* computing probs from the small N's */
void GetCoNeFileProbs(Nstruct *Ns, int Num, ldat *loca, int NumLoc, int T)
{
	int j,i,nn,k;
	double dT = (double)T;
	int *DoneIt = (int *)ECA_CALLOC(100000,sizeof(int));
	char FNAME[1000];
	FILE *in;
	char temp[100], temp2[100];
	
	
	gNumNs = Num + FILE_LINES;
	
	/* first allocate some space to CoalProbs */
	for(i=Num;i<gNumNs;i++)  {
		Ns[i].CP = (CoalProbs *)ECA_CALLOC(gMax_n0 + 1, sizeof(CoalProbs));
		Ns[i].LocLikes = (double *)ECA_CALLOC(NumLoc,sizeof(double));
		Ns[i].LocLikeVars = (double *)ECA_CALLOC(NumLoc,sizeof(double));
		for(j=0;j<=gMax_n0;j++)  {
			Ns[i].CP[j].P = NULL;
		}
	}
	
	
	/* first, allocate the memory */
	for(i=Num;i<gNumNs;i++) {
		for(j=0;j<NumLoc;j++)  {
			nn = loca[j].n0;
			if(Ns[i].CP[nn].P == NULL) {
				Ns[i].CP[nn].P = (double *)ECA_CALLOC(nn+1,sizeof(double));
			}
		}
	}
	
	/* then go get stuff out of the files */
	for(j=0;j<NumLoc;j++)  {
		nn = loca[j].n0;
		if(DoneIt[nn] == 0) {
			DoneIt[nn] = 1;
			sprintf(FNAME,"%s%dpr.txt",gProbsPath,nn);
			if( (in = fopen(FNAME, "r"))==NULL) {
				printf("\n\nError! Unable to find probs file:\n\t%s\nfor locus %d\n",FNAME,j+1);
				printf("You may have incorrectly specified the pathway-to-probs-files\n");
				printf("or you may not have the appropriate XXXpr.txt file.  See the documentation\n");
				printf("for instructions on using simCoNeprob to make new files as needed.\n\nExiting to system...\n\n");
				exit(1);
			}
			for(i=Num;i<gNumNs;i++) {
				/* eat the t:  and get the time t */
				fscanf(in," %s %lf",temp, &(Ns[i].t));
				
				/* then set the N value in there too */
				Ns[i].N = T / (2.0 * Ns[i].t);
				
				/* then eat and get the starts and ends */
				fscanf(in," %s %d %s %d", temp, &(Ns[i].CP[nn].start),temp2, &(Ns[i].CP[nn].end));
				
				/* then get all the rest of the goodies */
				for(k=Ns[i].CP[nn].start;k<=Ns[i].CP[nn].end;k++)  {
					fscanf(in, " %lf",&(Ns[i].CP[nn].P[k]));
				}
				
			}
			fclose(in);
			
		}
	}
	
	free(DoneIt);
}

/*
	This echoes the data back out in the format it was read in.
*/
void EchoCoNeData(ldat *loc, int NumLoc)
{
	int i,j;
	printf("DATFILE : 0\nDATFILE : 2\nDATFILE : %d\n",NumLoc);
	for(i=0;i<NumLoc;i++)  {
		printf("DATFILE : %d\nDATFILE : ",loc[i].K);
		for(j=0;j<loc[i].K;j++)  
			printf("%f ",loc[i].AT[j]);
		printf("\nDATFILE : ");
		for(j=0;j<loc[i].K;j++)  
			printf("%d ",loc[i].A0[j]);
		printf("\n");
	}
}





/* 
	Prints out the probs for the different Ns and the different loci
*/
void PrintNProbs(Nstruct *Ns, int NumNs, ldat *loc, int NumLoc, int Verbosity)
{
	int i,j,nf;
	double temp;
	
	/*  first summarize the stuff:  */
	printf("NeValues: ");
	for(i=0;i<NumNs;i++)  {
		printf("%f ",Ns[i].N);
	}
	printf("\n\n");
	
	if(Verbosity > 0)
	{
		/* cycle over the Ne values */
		for(i=0;i<NumNs;i++)  {
			printf("NeValue: %f\n",Ns[i].N);
			/* print a header line */
			printf("nf ");
			for(j=0;j<NumLoc;j++)  
				printf("loc%d ",j+1);
			printf("\n");
			
			for(nf=1;nf<gMax_n0;nf++)  {
				printf("%d ",nf);
				for(j=0;j<NumLoc;j++)  {
					if(nf > loc[j].n0)
						temp = 0.0;
					printf("%e ",temp);
				}
				printf("\n");
			}
			printf("\n");
		}
	}
	
}




/*
	This function will add NReps reps of importance sampling to the
	average held in the dval probs[nf] for locus data held in loc.
	Note that loc also holds the probabilities that will get stored for it.
*/
void AddReps(dval **probs, int nf, int NReps, ldat *loc)
{
	int n0,k,j,i;
	int K = loc->K;
	int Af[MAX_ALLELES];
	double predpars[MAX_ALLELES];
	int nonzero=0;
	double temp;
	
	/* count up n0  and set the values of predpars and count up nonzeros in A0*/
	for(n0=0,k=0;k<K;k++)  {
		n0 += loc->A0[k];
		nonzero += (loc->A0[k] > 0);
		if(gPrior==0.0) /* here is where we set the prior on allele freqs */
			predpars[k] = (loc->AT[k]) + 1.0; 
		else 
			predpars[k] = (loc->AT[k]) + gPrior / (double)K;    
	} 
		
	/* if loc->P is not allocated to, then take care of that */
	if(loc->P==NULL) 
		loc->P = (double ***)ECA_CALLOC( K, sizeof(double **));
	for(j=0;j<K;j++)  {
		if(loc->P[j]==NULL)
			loc->P[j] = (double **)ECA_CALLOC( n0+1, sizeof(double *));
	}
	
	/* if nf is too small, or too large, send back a 0.0 with zero variance */
	if( nf < nonzero || nf > n0) {
		probs[nf]->v = 0.0;
		probs[nf]->Ave = 0.0;
		probs[nf]->Var = 0.0;
		probs[nf]->NumAved = 1;
		return;
	}
	
	if(loc->NormoExtract[nf]==0.0) {  /* if it's zero, it hasn't been set yet, so we need to set it.  */
		/* we set it to be just the log of the value of the first importance sampling ratio */
		temp = SimConstrainedCMD(Af,nf,K, predpars,loc->A0,loc->P);
		loc->NormoExtract[nf] = LogPrA0givenAf(loc->A0,Af,K) + LogCMDPMF(Af,predpars,K) - temp;
	}
	
	/*  cycle over the reps */
	for(i=0;i<NReps;i++)  {				
		
		/* simulate the value of Af (keeping track of the prob of simulating it) */
		temp = SimConstrainedCMD(Af,nf,K, predpars,loc->A0,loc->P);

#ifdef DO_RENORMO_OVER_LOCI		
		/* then compute the importance sampling ratios....WITH A LITTLE CONSTANT TAKEN OUT!! */
		probs[nf]->v = exp(LogPrA0givenAf(loc->A0,Af,K) + LogCMDPMF(Af,predpars,K) - 
						temp  - loc->NormoExtract[nf]  );
#endif

#ifndef DO_RENORMO_OVER_LOCI 
	/* then compute the importance sampling ratios and don't take the constant out */
		probs[nf]->v = exp(LogPrA0givenAf(loc->A0,Af,K) + LogCMDPMF(Af,predpars,K) - 
						temp);
#endif

		
/* printf("DoingReps.  nf= %d  ImpRat= %f\n",nf,probs[nf]->v);  */
		
		/* then record it */
		IncrementDval(probs[nf]);
		
		/* if the number of alleles is 2, then one iteration is exact, so record 0 for the variance and
			get out */
		if(K==2) {
			probs[nf]->Var = 0.0;
			break;
		}
	}
}






double ComputeNeLikes(ldat *loc, Nstruct *Ns, int NumLoc)
{
	int i,l,nf,lo,hi,n0;
	double c,u,v;
	double MaxNormoExtract,TotOverLociExtract=0.0;  /* for dealing with renormalizations */
	double term1, term2;  /* for computing the variance of likelihood over all loci...
							these refer to the two terms in equation 13 in Anderson and Thompson */
	
	for(l=0;l<NumLoc;l++)  {  /* cycle over loci */
	
		n0 = loc[l].n0;


		/* *******************************   */
		/* the first thing we do is find the maximum NormoExtract */
		lo = loc[l].Kin0;
		hi = n0;
		MaxNormoExtract=loc[l].NormoExtract[lo];
		for(nf=lo;nf<=hi;nf++)  {
			if(MaxNormoExtract < loc[l].NormoExtract[nf])	
				MaxNormoExtract = loc[l].NormoExtract[nf];
		}
		TotOverLociExtract +=  MaxNormoExtract;

		for(i=0;i<gNumNs;i++)  {  /* cycle over the NumNs */
			/* initialize to accumulate sums */
			Ns[i].LocLikes[l] = 0.0;
			Ns[i].LocLikeVars[l] = 0.0;
			
			/* then cycle over possible (prob > 0) values of nf.  We cycle from lo to hi which are determined so
			   that we don't sum over a lot of zeros */
			lo = ECA_MAX(loc[l].Kin0, Ns[i].CP[n0].start);
			hi = ECA_MIN(n0,Ns[i].CP[n0].end);
#ifdef DO_RENORMO_OVER_LOCI
			for(nf=lo;nf<=hi;nf++)  {
				
				/* here we deal with renormalizations */
				c = Ns[i].CP[n0].P[nf] ;
				u = log(loc[l].R[nf]->Ave)  + loc[l].NormoExtract[nf] - MaxNormoExtract;
				v = log(loc[l].R[nf]->Var) + 2.0 * (loc[l].NormoExtract[nf] - MaxNormoExtract);
				
/*printf("BEFOR: Locus%d i=%d  N[i]=%.3f    nf=%d   MaxNormo=%f Ave=%f  Var=%f    u=%f   v=%f  loc[l].R[nf]->Ave=%e\n",
	l,i,Ns[i].N,nf,MaxNormoExtract,log(loc[l].R[nf]->Ave),log(loc[l].R[nf]->Var),u,v,loc[l].R[nf]->Ave );	*/
			
				/* here we call any probability less than e^-345 or variance less than e^-690 a 0.0 
					except for the case in which the variance itself is exactly 0.0 (because there is
					only one configuration of a_f that works  */
				if(u < -345.0  || (v < -690.0 && loc[l].R[nf]->Var != 0.0) ) {
					u = 0.0;
					v = 0.0;
				}
				else if (loc[l].R[nf]->Var == 0.0) {
					u = exp(u);
					v = 0.0;
				}
				else {
					u = exp(u);
					v = exp(v);
				}
				
/* printf("AFTER: Locus%d i=%d  N[i]=%.3f    nf=%d   MaxNormo=%f Ave=%f  Var=%f    u=%e   v=%e \n",
	l,i,Ns[i].N,nf,MaxNormoExtract,log(loc[l].R[nf]->Ave),log(loc[l].R[nf]->Var),u,v );  */
				
				Ns[i].LocLikes[l] += c * u;
				Ns[i].LocLikeVars[l] += c * c * v;
				
			}
#endif
#ifndef DO_RENORMO_OVER_LOCI
			/*   UN-RENORMALIZING VERSION: */
			for(nf=lo;nf<=hi;nf++)  {
				c = Ns[i].CP[n0].P[nf];
				u = loc[l].R[nf]->Ave;
				v = loc[l].R[nf]->Var;
				
				Ns[i].LocLikes[l] += c * u;
				Ns[i].LocLikeVars[l] += c * c * v;
				
			} 
#endif
		}
		
	}
	
	
	/* now we print out some intermediate output about the log-like from each locus */
	printf("LOC_SPECIFIC_LOGLS  :   Ne     t    ");
	for(l=0;l<NumLoc;l++)
		printf("Loc%d  ",l+1);
	printf("\n");
	for(i=0;i<gNumNs;i++)  {
		printf("LOC_SPECIFIC_LOGLS   :    %.3f    %f   ",Ns[i].N, Ns[i].t);
		for(l=0;l<NumLoc;l++) {	
			printf("%.3f   ",log(Ns[i].LocLikes[l]));
		}	
		printf("\n");
	}
	
	
	
	/* now we have all the likelihoods for the various Ne's at different loci, and 
	   we have to combine those all into a single overall Ne.  This is not so hard...
	   it is just the product over the loci.
	 */
	 for(i=0;i<gNumNs;i++)  {
	 	Ns[i].Like = 1.0;
	 	
	 	term1 = 1.0;
	 	term2 = 1.0;
	 	for(l=0;l<NumLoc;l++)  {
	 		c = Ns[i].LocLikes[l];
	 		v = Ns[i].LocLikeVars[l];
	 		
	 		Ns[i].Like *= c;
	 		
	 		term1 *= (c*c);
	 		term2 *= ( c*c - v ); 
	 	}
	 	Ns[i].LikeVar = term1 - term2;
	 
	 }
	 return(TotOverLociExtract);
}
	



/* this ignores alleles in which the corresponding elements in 
both a0 and af are zero */
double LogPrA0givenAf(int *a0, int *af, int K)
{
	double res = 0.0;
	int k,n0=0,nf=0;
	
	for(k=0;k<K;k++)  {
		if( !(a0[k]==0 && af[k]==0) ) {
			res += LogFactorial(a0[k]-1) - LogFactorial(af[k]-1) - LogFactorial(a0[k]-af[k]);
			n0 += a0[k];
			nf += af[k];
		}
	}

	res -= 	LogFactorial(n0-1) - LogFactorial(nf-1) - LogFactorial(n0-nf);
	
	return(res);
}


/*
	This is a function for simulating Af from a distribution that is proportional to the
	product of a compound multinomial dsn with parameters beta with the distribution of
	A0 given af,  and with the total
	number of balls nf, but subject to the constraint that it is compatible with A0 (so, no element
	of Af can be greater than the corresponding element in A0, and an element of Af cannot be zero if
	it is non-zero in A0).  
	
	The array P is a stored array of cumulative probability functions for this distribution.
	
	This function returns the log probability of simulating the value of Af.  Af is an output variable.
	
	Currently, the max value that K can take is 500.
*/ 
double SimConstrainedCMD(int * Af, int nf, int K, double *beta, int *A0, double ***P)
{
	int i,j, lo, hi;
	int A0r[500], A0i[500];  /* reverse cumulative counts of the elements of A0, and the reverse
									cumulative count of the indicators that each element of A0>0.  */ 
	double b1,b2; /* To store the parameters of the component beta-binomials */
	double rando,div,simprob=0.0, prob;
	int n0,a0;  /* temp variables for dealing with the marginal beta-binomials */
	
	/* store the reverse cumulative arrays */								
	A0r[K-1] = A0[K-1];
	A0i[K-1] = (A0[K-1] > 0);
	for(i=K-2;i>=0;i--) {
		A0r[i] = A0[i] + A0r[i+1];
		A0i[i] = (A0[i] > 0) + A0i[i+1];
	}
			
			
	/* Cycle over the components, and determine the constraints, and then simulate them */
	for(i=0;i<K-1;i++)  {
		/* set single component constraints first */
		lo = (A0[i] > 0);
		hi = (A0[i]);
		
		/* then set the "multi-component" constraints, if they are more extreme
		  (see Notes: CoNe 10/29/03--10/30/03 page 2 */
		if( lo < nf - A0r[i+1]) lo = nf - A0r[i+1];
		if( hi > nf - A0i[i+1]) hi = nf - A0i[i+1];
		
		/* here, if lo > hi there is a problem.  Exit */
		if(lo>hi) {
			printf("\nlo>hi in SimConstrainedCMD. lo=%d hi=%d, i=%d, nf=%d\n\n",lo,hi,i,nf);
			printf("\nExiting\n\n");
			exit(1);
		}
		else if(lo == hi) {  /* if hi==lo, there is only one thing Af[i] can be and the prob of that is 1.0 */
			Af[i] = lo;
			prob = 1.0;
		}
		else {  
			
			/* check to see if the corresponding array in P is filled out, if not
				then fill it */
			if(P[i][nf]==NULL) {
				/* compute the corresponding parameters for the marginalized beta-binomial */
				b1 = beta[i];
				b2 = 0.0;
				for(j=i+1;j<K;j++)  {
					b2 += beta[j];
				}
				
				/* compute the corresponding marginalized A0's too */
				a0 = A0[i];
				n0 = a0;
				for(j=i+1;j<K;j++)  {
					n0 += A0[j];
				}
				/* then allocate to and fill the array */
				P[i][nf] = (double *)ECA_CALLOC((size_t)nf+1, sizeof(double));
				ProbAfGivenA0_CDF_Array(b1,b2,nf,P[i][nf],a0,n0);
			}
			
			/* now we have to draw a random number and scale it so it will pick out an
				element in P that is between lo and hi inclusive. */
			rando = (double)ranf();
			div = P[i][nf][hi];
			if(lo>0) div -= P[i][nf][lo-1];
			rando /= div;
			if(lo>0) rando += P[i][nf][lo-1];
			
			/* then call the binary search routine which returns the index of the 
				least element in an array of sorted doubles which is greater than rando.
				However, if hi-lo is less than 50, then just cycle over the probs... */
			if(hi-lo < 50) {
				for(j=lo;j<hi;j++) {
					if(P[i][nf][j] > rando)
						break;
				}
				Af[i] = j;
			}
			else {
				Af[i] = IntFromArray_BinarySearch(rando,P[i][nf],lo,hi);
				if(Af[i] < 0) {
					printf("\n\nIntFromArray_BinarySearch() returned %d\n\n",Af[i]);
					exit(1);
				}
			}
			
			/* then, compute the probability of realizing that and add it to the running total */
			if(Af[i]==0) 
				prob = P[i][nf][0];
			else
				prob = P[i][nf][Af[i]] - P[i][nf][Af[i]-1];
			prob /= div;
		}
		
		simprob += log(prob);
		
		/* now we must decrement nf by the value of Af[i] just realized */
		nf -= Af[i];
		
	}				
	
	/* the remainders go into Af[K-1] */
	Af[K-1] = nf;		
									
	return(simprob);									
									
									
}


/*
	This puts the Cumulative distribution function for P(af|a0,nf,aT)	
	where af is a beta binomial r.v. with
	parameters a and b (both positive), and n balls drawn total,  a0|af is like a hypergeometric. (see notes 
	and paper), into the array P.
	
	It does not allocate memory to P.  That must be done ahead of time.  It must have
	space allocated from subscript 0 to subscript n.  (that is it needs n+1 elements).
	
	At the end P[0] is Prob(X=0)
			   P[1] is Prob(X<=1)
			   .
			   .
			   .
			   P[j] is Prob(X <= j).
			  
	This computes things a little recursively.
	
	The usual normalization factors are left out of this function, because it
	all gets normalized in the end anyway.
	
	n plays the role of nf here.
	
*/
void ProbAfGivenA0_CDF_Array(double a, double b, int n, double *P, int a0, int n0)
{
	double y=1.0,q,j=0.0;
		int i=0;
	int bot, top;
	double hg;  /* product of binom coefficients for the hypergeometric part */
	int tempa0[2], tempaf[2];
	int r;
	double renormo;
	double pre_y;
	
	/* NOTE:  
		y is the "left hand" factor in the beta-binomial 
		q is the "right hand" factor in the beta-binomial
		
		bc1 is the left side binomial coefficient in the hypergeometric
		bc2 is the right side binom coeff in the hypergeometric
	*/
	
	/* the limits on the value of af: */
	bot = ECA_MAX(1,n-(n0-a0));
	top = ECA_MIN(a0,n-1);
		 


/* this is starting point for working backward in the beta-binomial part:
	
	q = exp( LogGamma(n+b) - LogFactorial(n) - LogGamma(b));
	
	HOWEVER BECAUSE WE ARE GOING TO NORMALIZE EVERYTHING IN THE END,
	WE CAN JUST START IT AT 1.  THE CONSTANTS ALL GET PULLED OUT!!
*/
	q = 1.0;
	P[0] = 0.0;  /* in this situation, af must be greater than 0 */  
	y=1.0;
	j=0.0;
	while(i<top) {
		j += 1.0;
		i++;
		
		/* this is for the beta-binomial part */
		pre_y = y;
		y *= (a + j - 1.0)/j;
		q *= ( 1.0*n - j + 1.0 ) / (b + 1.0*n - j);
		
		
		/* this is for the hypergeometric part */
		if(i<bot) {
			P[i] = 0.0;
		}
		else if(i==bot) {
			/* the following is the beginning for the hg recursion:
			
			hg = exp( LogFactorial(a0 - 1) - LogFactorial(i-1) - LogFactorial(a0-i)  + 
					LogFactorial(n0 - a0 - 1) - LogFactorial(n - i - 1) - LogFactorial(n0 - a0 - n + i)   );
				
				BUT, ONCE AGAIN, WE CAN JUST START IT AT 1.0, AND THEN WE ARE LESS LIKELY TO GET OVERFLOW PROBLEMS	
			*/
			hg = 1.0;
			tempa0[0]=a0;
			tempa0[1]=n0-a0;
			tempaf[0]=i;
			tempaf[1]=n-i;
			/*hg = exp(LogPrA0givenAf(tempa0,tempaf,2));  I had this in there to test that it was working OK
				It was...but the recursion method will be a little faster */
		}
		else {  /* when i > bot */
			hg *= ((a0-j+1.0)*(n-j)) / ( (j-1.0)*(n0-a0-n+j) );
			tempaf[0]=i;
			tempaf[1]=n-i;
			/* hg = exp(LogPrA0givenAf(tempa0,tempaf,2));  I had this in there to test that it was working OK
				It was...but the recursion method will be a little faster */
		}
		
		if(i>=bot) {  /* here we do the cumulative probability calculation */
			P[i] = P[i-1] + (y * q * hg);
			if(!isfinite(P[i])  || !isfinite(y*q*hg)  || !isfinite(y) || !isfinite(q) || !isfinite(hg) ) {
				printf("\n\nOverflow in ProbAfGivenA0_CDF_Array()   bot= %d   top= %d \n", bot,top);
				printf("Parameters:  a= %e    b= %e        n= %d    a0= %d    n0= %d\n", a,  b,  n,  a0,  n0);
				printf("\ni is %d  P[i-1] is %e  P[i-2] is %e (y*q*hg) = %e  y= %e  q= %e  hg= %e\n",i, P[i-1],P[i-2],y*q*hg,y,q,hg );
		if(P[i-1] == 0.0) {
			for(r=0;r<=i;r++)  {
				printf("\n\tP[%d]= %e    i= %d",r,P[r],i);
			}
			printf("\n");
		}
				printf("Attempting Renormalization");
				renormo = -log(P[i-1]);  /* this is what we will multiply each cell by */
				for(r=1;r<i;r++) {
					if( P[r]/P[i-1] < 1.0e-300 ) 
						P[r] = 0.0;
					else
						P[r] /= P[i-1];
				}
				renormo += log(pre_y) + log((a + j - 1.0)/j) + log(q) + log(hg);  /* now renormo holds the log of the renormalized value of the non-cumulative part of P[i] */
				renormo = pow(exp(renormo), 1.0/3.0);  /* divvy up the difference between y,q, and hg */
				y=q=hg = renormo;
				P[i] = P[i-1] + (y * q * hg);
				printf("\nRenormalized: P[i-1] = %e  (y*q*hg) = %e  y= %e  q= %e  hg= %e  P[i]= %e\n", P[i-1],y*q*hg,y,q,hg,P[i] ); 
			}
		}
		
/*		printf("bot=%d  top=%d || i=%d   hg=%e   y=%e   q=%e\n",bot,top,i,hg,y,q);   */
		
	}
	/* fill the rest, from top+1 to n in with P[top] */
	for(i=top+1;i<=n;i++)  {
		P[i] = P[top];
	}
	
	/* now, normalize them all */
	for(i=0;i<=n;i++)  {
		P[i] /= P[n];
	}
	
}

double Rrecurs(int ii,int ji,int ki, double t) 
{
	double i,j,k;
	double prod = -1.0;
	
	i = (double)ii;
	k = (double)ki;
	j = (double)ji;
	
	
	
	prod *= (2*k-1.) / (2*k-3.);
	
	prod *= (i-k+1.) / (i+k-1.);
	
	prod *= (j+k-2.) / (k-j);
	
	prod *= exp( -t * (k-1.));
	
	return(prod);
}

double Wrecurs(int ii, int ji, double t)
{
	double i,j;
	double prod = 1.0;
	
	i = (double)ii;
	j = (double)ji;
	
	prod *= exp(-t * j);
	
	prod *= (2.0*j-3.0) / (2.0*j-1.0);
	
	prod *= (4.0*j-2.0) / (j+1.0);
	
	prod *= (i-j) / (i+j);
	
	return(prod);
}

double ExpectedNumberOfLineages(int i, double t)
{
	double k,j,di;
	double w = 1.0, sum = 0.0;
	
	
	di = (double)i;
	sum += w;
	for(j=2,k=2.0 ; j<=i ; j++, k+=1.0)  {
		w *= exp(-t * (k-1)) * (di-k+1.0)/(di+k-1.0)  * (2.0*k-1.0)/(2.0*k-3.0);
		sum += w;
	}
	return(sum);
}


/* 
	This computes the an unnormalized probability recursively
	for any j between 2 and i, inclusive.  The variable Const
	is the value given to the first term in the sum (the anchor
	of the recursion).
	
	This takes account of numerical precision issues by returning
	0.0 any time the max(abs(summand)) is 13 orders of magnitude greater
	than the sum. 
*/
double Tavare1984recurs(int i, int j, double t, double Const)
{
	int k;
	double sum=Const;
	double w = Const;
	double the_max = Const;
	double abs_of_it, sum_abs;
	
	if(j<2 || j>i) {
		printf("\n\nj out of range in Tavare1984recurs().  \nExiting...\n");
		exit(1);
	}
	
	for(k=j+1;k<=i;k++)  {
		w *= Rrecurs(i,j,k,t);
		abs_of_it = ECA_ABS(w);
		if(abs_of_it > the_max) the_max = abs_of_it;
		sum += w;
if(0)		printf("\tj= %d  w= %e  sum= %e\n",j,w,sum); 
	}	
	
	sum_abs = ECA_ABS(sum);
	
	if( log10( the_max / sum_abs)  > 15 )  sum=0.0;
	return(sum);
}


/* this is a version that doesn't ask for a starting point.  It goes through and tries to find 
 one */
double Tavare1984by_recurs(int i, int j, double t)
{
	int k=j;
	double sum=0.0,w;
	
	if(j==1) {
		sum = Tavare1984prob(i,j,t);
	}
	
	else {
		w = Tavare1984summand(t,i,j,j);
		while( (w == 0.0  || w == -0.0)  && k<i) {  /* work up to a non-zero starting point */
			k++;
			w = Tavare1984summand(t,i,j,k);
		}
		sum+=w;
		while( k < i) {
			k++;
			w *= Rrecurs(i,j,k,t);
			sum += w;
		}
	}
	
	if(sum < 0 && ECA_ABS(sum) < .000001)
	sum *= -1.0;
	return(sum);
}