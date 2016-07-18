/* comment out this #define line to avoid the renormalizing involving the
ExtractNormo */ 
#define DO_RENORMO_OVER_LOCI

#define MAX_ALLELES 500
#define FILE_LINES 100   /* the number of lines in the XXp.txt files. */

/* if the log of the probability ratios reaches this point just call it zero */
#define CALLITZERO -40


typedef struct {
	int K; /* number of alleles */
	int Kin0;  /* number of alleles found in the time 0 sample */
	int n0;
	double nT;
	int *A0; 
	double *AT;
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
void GetCoNeData(char *FILENAME, ldat **loc, char **argv, int argc, int *NumLoc, double *T,  int *NumReps, char *phrase, Nstruct **Ns);
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
GLOB int gMax_n0;  /* for the largest n0 out of all the loci */
GLOB int gNumNs;  /* the nuumber of N's we are keeping track of */
GLOB char gProbsPath[10000];	
GLOB double gPrior;