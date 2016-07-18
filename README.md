# cone  --  a program for estimating Ne

To get and compile this do like so:
```sh
  # first clone it from GitHub
  git clone https://github.com/eriqande/CoNe

  # then get the submodules
  cd CoNe
  git submodule init
  git submodule update

  # then compile it.  The result is the executable with the 
  # operating system appended. For example on my Mac I get
  # the binary ms2geno-Darwin
  ./Compile_CoNe.sh 


  # You can do this to get the full built-in help
  CoNe-Darwin --help-full
```

Here is what the built in help tells you, which should be sufficient to get going
with the program:

```
cone  --  a program for estimating Ne

Author(s):
	Eric C. Anderson (eric.anderson@noaa.gov)

About the Program:
    CoNe computes the likelihood of Ne given data on two temporally spaced genetic
    samples. The statistical model used is based on the coalescent of the gene
    copies drawn in the second sample, as described in Berthier et al. (2003)
    Genetics 2003. 160:741-51. The Monte Carlo computations to compute the
    likelihood, however, were developed by Eric Anderson, and are orders of
    magnitude faster than previous implementations. 
    
    Details of the algorithm are given in Anderson (2005) Genetics 170:955-967.

In the following:
	"J" refers to an integer argument to an option

	"R" refers to a real number argument to an option

	"S" refers to a string argument to an option

	"F" refers to a file path argument to an option. For example,
		"datfile.txt" if the file is in the current working directory, or
		something like "~/eriq/Documents/data/datfile.txt" if you want to
		provide a complete file path.  (Beware of spaces in file paths!)

	"D" refers to a directory path argument to an option. For example,
		"data_direcory/" if the directory is in the current working directory, or
		something like "~/eriq/Documents/data_directory/" if you want to
		provide a complete directory path.  Note that the trailing slash should be
		optional, but currently is not.  (ERIC ADD MACROS FOR GETTING FILES AND DIRECTORIES

	"G" refers to a string that gives a (possibly) discontinous range of
		nonnegative integers.  For example:  "1-5,7,9,10-15" specifies
		the integers 1 through 5, 7, 9, and 10 through 15.  There can be no
		whitespace in the string specifying the range, and the numbers must all
		be increasing.  Also, the string cannot start nor end with a comma or a dash.
		Finally, you should not use "-" to denote two ranges without placing any
		commas in between.

	"C" refers to a "constrained" string argument to an option,
		i.e., the argument is a string that may only be drawn from a small
		set of alternatives, as specified in the help-full description.

   ****  Program-description and documentation  parameters  ****

     --help                 
        this returns a short list of all program options and associated
        arguments

     --help-full            
        this returns a full list of all program options and associated
        arguments

     --help-nroff           
        this returns a full list of all program options and associated
        arguments using the formatting styles in nroff that give you the look
        of a man page. View the formatted ouput by doing: 'prog --help-nroff
        | nroff -man | more' where prog is the name of the program.

     --help-xml             
        This returns a list of all options in a simple XML format which is
        suitable for input to the guiLiner front end.

     --version              
        prints program version information

     --version-history      
        prints history of different versions

     --command-file         F
        By using this option you can store command line arguments in the file
        named in F. You may have any kind of white space (including line
        endings) in the file. The line endings are treated as just another
        space. Comments may be included in the file by enclosing them within
        a pair of ampersands (the & character). Note that you must have a &
        at the beginning and at the end of the comment. You cannot use just a
        single & to comment to the end of a line. Your comments may spread
        over several lines---they will still be stripped from the resulting
        command line so long as the are enclosed in ampersands. This feature
        is helpful if you have long and complex command lines that you wish
        to store if it makes it easier to read the command line by breaking
        it across multiple lines or if you have certain clusters of options
        that you like to store together as a module. This option may be used
        up to 10000 times. Optional.


   ****  Data Analysis Options  ****

-f , --file-name            F
        F is the name of the file in which you have your data. It is in the
        same format as data files for TM3 (by Pierre Berthier) and TMVP (by
        Mark Beaumont) The data file should start with a 0 (this is a strange
        vestige of some sort from TM3 or TMVP) followed by the number of time
        periods ***which in this case must always be 2*** followed by the
        number of loci. Then data for each locus consists of the number of
        alleles observed at the locus followed by a row of counts of the
        different alleles observed in the first sample and a row of counts of
        the different alleles observed at the second sample. (By first sample
        I mean the sample taken first in time going forward. Hence, the
        second sample is the sample that was collected most recently.) There
        must be only integers and whitespace in the file. An example file is
        shown in the FILES section of the manual pages. NOTE! It turns out
        that it is not essential that the counts of alleles observed in the
        first sample (the one further back in the past) be integers. It turns
        out to be convenient to express them as real numbers in some cases,
        so I have recoded it so that they can be real numbers and not just
        integers. Now when the data file is echoed to standard input, these
        counts are expressed as real numbers. Do not let this alarm you.
        Everything is OK, still.

-p , --path-to-probs-files  D
        This is the pathway to files containing the precomputed probabilities
        of having j lineages remaining at scaled time t given that you
        started with i lineages. Note that the trailing slash is required.
        For example on my system ~/Documents/eca_code/CoNe/probs/ is the
        pathway. The probs files are a collection of files named XXXpr.txt
        where XXX is a number giving the number of gene copies in the second
        sample. These files have been precomputed using the program
        simCoNeprob and are described below in FILES. The CoNe distribution
        includes precomputed XXXpr.txt files with XXX ranging from 10 to 400.
        This represents samples of between 5 and 200 diploid organisms. For
        different sample sizes it is necessary to create new XXXpr.txt files
        using simCoNeprob which is also included with the CoNe distribution.

-T , --gens-between         R
        R is the number of generations between samples. It may be specified as
        a non-integer in order to allow for a non-integer number of
        generations.

-m , --mc-reps              J
        J is the number of importance sampling reps to perform for each value
        of the number of genes ancestral to the second sample on each locus.
        I have found the importance sampling algorithm to be good enough that
        100 J=100 gives reliable results and usually runs very quickly.
        However; on a final run of your data set J should be much larger. You
        can get an idea of whether J should be larger by the width of the
        Monte Carlo confidence intervals around the estimated likelihood
        curve.

-n , --ne-lo-hi-step        R1 R2 R3
        sets the values of Ne for which the likelihood will be computed. R1 is
        the lowest value of Ne. R2 is the highest value of Ne. R3 is the step
        size between values of Ne. For values of Ne such that T/(2Ne) is
        smaller than .06 (or so) the precomputed values of scaled time
        (stored in the appropriate XXXpr.txt file (see the description of the
        -p option)) will be used. So the step size may not be that given by
        R3. Arguments need not be given as real numbers. An integer (like
        250) will work just fine.

-q , --prior                R
        R specifies the prior distribution for the allele frequencies. R=0
        makes the prior a uniform Dirichlet distribution. R>0 makes the
        parameters of the Dirichlet prior are R/K. Hence R=1 is the so-called
        unit information prior. R=1 is the default.

-s , --seed-phrase          S
        S is a single string (no spaces) that will be used to seed the random
        number generator. If this option is not invoked then a seed is chosen
        based on the current time or---if the file cone_seeds is
        present---the seeds are taken from that file. Upon completion of the
        program the next random number seeds in series are printed to the
        file cone_seeds.

```