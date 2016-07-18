#define UN_EXTERN

/* comment out this #define line to avoid the renormalizing involving the
ExtractNormo */ 
#define DO_RENORMO_OVER_LOCI

#define MAX_ALLELES 500
#define FILE_LINES 100   /* the number of lines in the XXp.txt files. */

/* if the log of the probability ratios reaches this point just call it zero */
#define CALLITZERO -40

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "MathStatRand.h"
#include "ECA_MemAlloc.h"
#include "ranlib.h"
#include "MCTypesEtc.h"

typedef struct {
	int K; /* number of alleles */
	int Kin0;  /* number of alleles found in the time 0 sample */
	int n0;
	int nT;
	int *A0; 
	int *AT;
	double *predpars;
	dval **R; /* an array of the MC estimates of P(A0|nf)  */
	double ***P;  /* the stored probabilities, indexed by [allele][nf][af] */
	double *NormoExtract;   /* for a normalization to pull out of each value so they sum to one */
	double LocOverallNormo;
} ldat;  /* locus data */

typedef struct {
	int start;
	int end;
	double *P;
	
} CoalProbs;


typedef struct {
	double N;
	double t;    /* the scaled time that this corresponds to */
/*	double **P;  /* for pointing to the probabilities, one array fo */
	double *LocLikes;  /* log likelihood for the individual loci */
	double *LocLikeVars;
	double Like;
	double LikeVar;
	CoalProbs *CP;
} Nstruct;



void GetCoNeFileProbs(Nstruct *Ns, int Num, ldat *loca, int NumLoc, int T);
void PrintNsLong(Nstruct *Ns);
void GetCoNeData(ldat **loc, char **argv, int argc, int *NumLoc, double *T,  int *NumReps, char *phrase, Nstruct **Ns);
void EchoCoNeData(ldat *loc, int NumLoc);
void PrintNProbs(Nstruct *Ns, int NumNs, ldat *loc, int NumLoc, int Verbosity);
double LogPrA0givenAf(int *a0, int *af, int K);
double ComputeNeLikes(ldat *loc, Nstruct *Ns, int NumLoc);
double SimConstrainedCMD(int * Af, int nf, int K, double *beta, int *A0, double ***P);
void ProbAfGivenA0_CDF_Array(double a, double b, int n, double *P, int a0, int n0);
void AddReps(dval **probs, int nf, int NReps, ldat *loc);
void PrintUsage(char *msg);
double Rrecurs(int ii,int ji,int ki, double t);
double Wrecurs(int ii, int ji, double t);
double Tavare1984recurs(int i, int j, double t, double Const);
double ExpectedNumberOfLineages(int i, double t);
double Tavare1984by_recurs(int i, int j, double t);

/* some global variables */
int gMax_n0;  /* for the largest n0 out of all the loci */
int gNumNs;  /* the nuumber of N's we are keeping track of */
char gProbsPath[10000];	
double gPrior;
	
	
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
	
	gMax_n0 = 0;
	
	if(argc == 1) {
		PrintUsage(NULL);
	}
	else {
		
		GetCoNeData(&locarray,argv,argc,&NumLoc,&T,&NumReps, phrase, &Ns);
		GetCoNeFileProbs(Ns, gNumNs, locarray,NumLoc,T); 

	}
	
	EchoCoNeData(locarray,NumLoc);
		
	/* PrintNProbs(Ns, gNumNs, locarray, NumLoc, 1); */
/*	PrintNsLong(Ns);  */
	
	
	printf("\n\tNumReps: %d",NumReps);
	phrtsd(phrase,&seed1, &seed2);
	setall(seed1, seed2);
	printf("\n\tSeedPhrase: %s", phrase);
	printf("\n\tSeeds:  %ld %ld",seed1,seed2);
	printf("\n\n");
	

	for(i=0;i<NumLoc;i++)  {  double maxlogprob=-9999999.9;
		printf("LOCUS: %d\n",i+1);
		printf("nf Ave  Var  SD  NormoExtract   expNormoExtract  LogAvePlusNormoExtract  \n");
		for(nf=locarray[i].Kin0; nf<=locarray[i].n0; nf++) {
			AddReps(locarray[i].R,nf,NumReps, &(locarray[i]));
			printf("%d  %e   %e     %e   %f   %e    %f\n",nf, locarray[i].R[nf]->Ave, locarray[i].R[nf]->Var, 
					sqrt(locarray[i].R[nf]->Var),locarray[i].NormoExtract[nf],
					exp(locarray[i].NormoExtract[nf]), log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf]);
			
			/* we are going to cut this off if the prob has dropped very low */
			if(  maxlogprob < log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf])
				maxlogprob = log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf];
			if(log(locarray[i].R[nf]->Ave)+locarray[i].NormoExtract[nf] - maxlogprob < CALLITZERO) {
				printf("Truncating evaluation beyond this point.  Calling it zero...\n");
				while(nf<=locarray[i].n0) {
					locarray[i].R[nf]->Ave = 0.0;
					locarray[i].R[nf]->Var = 0.0;
					locarray[i].NormoExtract[nf] = -999999.9;
					nf++;
				}
			} 
			

		}
	}
	
	
	TotOverLociExtract = ComputeNeLikes(locarray,Ns,NumLoc);
	
	printf("\nOverall Results:\n");
	printf("Ne  t  LogLike  Like  SD\n");
	for(i=0;i<gNumNs;i++)  {
		printf("%.3f  %f  %f   %e     %e\n",Ns[i].N,Ns[i].t, 
						log(Ns[i].Like),Ns[i].Like,sqrt(Ns[i].LikeVar));
		if(MaxN_LogLike < log(Ns[i].Like) ) {
			MaxN_LogLike = log(Ns[i].Like);
			MaxN_index = i;
		}
	}
	
	
	/* now report the max (by parabolic interpolation) and the 1.96 units of 
		log likelihood "confidence intervals" (by linear interpolation).  */
	if(MaxN_index == 0  || MaxN_index == gNumNs-1) {
		printf("MaxByParabolicInterplotation:  -999.999  -999.999  -999.999\n");
	} 
	else {
		i=MaxN_index;
		y = ParabolicMaxInterp(Ns[i-1].N,Ns[i].N,Ns[i+1].N,
				log(Ns[i-1].Like),log(Ns[i].Like),log(Ns[i+1].Like), &x);
		printf("MaxByParabolicInterpolation:  %f   %f    %f  ( N t LogLike )\n",x,(double)T/(2.0*x), y);
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
	printf("LowerSupportLimit: %f   %f  ( N  t )\n",lconf, T/(2.0*lconf) );
	printf("UpperSupportLimit: %f   %f  ( N  t )\n",uconf, T/(2.0*uconf) );
	
	printf("ThisAmountTakenOutOfLogLikelihood:  %f\n", TotOverLociExtract);
	printf("PriorUsedForAlleleFrequencies_[PRIOR]:  %f\n",gPrior);
	return(0);
}


void GetCoNeData(ldat **loc, char **argv, int argc, int *NumLoc, double *T,  int *NumReps, char *phrase, Nstruct **Ns)
{
	int i,j,l,K,n0,nT,temp;
	FILE *in;
	char INPUT, FILENAME[500];
	int NumSam;
	double Nlo,Nhi,Nstep;
	
	/* test to see if it is command line for one locus, or file for multiple loci */
	i=1;
	sprintf(gProbsPath,"%s",argv[i++]);
	if(argv[i][0] == '-') {
		INPUT = argv[i][1];
		i++;
	}
	else {
		PrintUsage("You've gotta give us the -f or the -c option!");
	}
	
	if(INPUT=='c') {
		*NumLoc=1;
		*loc = (ldat *)ECA_CALLOC(*NumLoc, sizeof(ldat));
		(*loc)[0].K = atoi(argv[i++]);
		K = (*loc)[0].K;
		(*loc)[0].A0 = (int *)ECA_CALLOC(K,sizeof(int));
		(*loc)[0].AT = (int *)ECA_CALLOC(K,sizeof(int));
		(*loc)[0].predpars = (double *)ECA_CALLOC(K,sizeof(double));
		for(nT=0,j=0;j<K;j++)  {
			(*loc)[0].AT[j] = atoi(argv[i++]);
			nT += (*loc)[0].AT[j];
			(*loc)[0].predpars[j] = (*loc)[0].AT[j] + 1.0;
		}
		(*loc)[0].Kin0 = 0;
		for(n0=0,j=0;j<K;j++) {
			(*loc)[0].A0[j] = atoi(argv[i++]);
			n0 += (*loc)[0].A0[j];
			(*loc)[0].Kin0 += ((*loc)[0].A0[j] > 0);
		}
		(*loc)[0].P = NULL;
		(*loc)[0].R = DvalVector(0,n0,0.0,0.0,0.0);
		(*loc)[0].n0 = n0;
		(*loc)[0].nT = nT;
		(*loc)[0].NormoExtract = (double *)ECA_CALLOC(n0+1,sizeof(double));
		if(gMax_n0 < n0) gMax_n0 = n0;
		
		*T = atof(argv[i++]);
		*NumReps = atoi(argv[i++]);
		sprintf(phrase,"%s",argv[i++]);
		Nlo = atof(argv[i++]);
		Nhi = atof(argv[i++]);
		Nstep = atof(argv[i++]);
		if(i==argc)  {  /* in this case, there is not  specification for PRIOR */
			gPrior = 1.0;
			printf("Using Default Unit-Information Prior for Allele Frequencies\n");
		}
		else if(i<argc) {
			gPrior = (double)atof(argv[i++]);
			if(gPrior==0.0) {
				printf("Using Uniform Prior for Allele Frequencies\n");
			}
			else {
				printf("Using Dirichet parameters %f/K for Allele Frequencies\n",gPrior);
			}
		}	
		
		
		
	}
	if(INPUT=='f')  {
		sprintf(FILENAME,"%s",argv[i++]);
		*T = atof(argv[i++]);
		*NumReps = atoi(argv[i++]);
		sprintf(phrase,"%s",argv[i++]);
		Nlo = atof(argv[i++]);
		Nhi = atof(argv[i++]);
		Nstep = atof(argv[i++]);
		if(i==argc)  {  /* in this case, there is not  specification for PRIOR */
			gPrior = 1.0;
			printf("Using Default Unit-Information Prior for Allele Frequencies\n");
		}
		else if(i<argc) {
			gPrior = (double)atof(argv[i++]);
			if(gPrior==0.0) {
				printf("Using Uniform Prior for Allele Frequencies\n");
			}
			else {
				printf("Using Dirichet parameters %f/K for Allele Frequencies\n",gPrior);
			}
		}	
		
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
			(*loc)[l].AT = (int *)ECA_CALLOC(K,sizeof(int));
			(*loc)[l].predpars = (double *)ECA_CALLOC(K,sizeof(double));
			for(j=0,nT=0;j<K;j++)  {
				fscanf(in,"%d",&((*loc)[l].AT[j]));
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
				printf("\n\nError! Unable to find probs file:\n\t%s\nfor locus %d\nExiting to system!\n",FNAME,j+1);
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
	printf("0\n2\n%d\n",NumLoc);
	for(i=0;i<NumLoc;i++)  {
		printf("%d\n",loc[i].K);
		for(j=0;j<loc[i].K;j++)  
			printf("%d ",loc[i].AT[j]);
		printf("\n");
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
			predpars[k] = (double)(loc->AT[k]) + 1.0; 
		else 
			predpars[k] = (double)(loc->AT[k]) + gPrior / (double)K;    
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