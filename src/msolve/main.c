/* This file is part of msolve.
 *
 * msolve is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * msolve is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with msolve.  If not, see <https://www.gnu.org/licenses/>
 *
 * Authors:
 * Jérémy Berthomieu
 * Christian Eder
 * Mohab Safey El Din */

#include "libmsolve.c"

#define DEBUGGB 0
#define DEBUGBUILDMATRIX 0
#define IO_DEBUG 0

#define LONG_OPT_LENGTH 15
#define ARG_OPT_LENGTH 4

static inline void display_option_help(char short_opt, char *long_opt,
				       char *arg_opt, char* str) {
  int long_opt_non_empty= strcmp (long_opt, "");
  
  if (short_opt == '\0') {
    fprintf (stdout, "    ");
  } else {
    fprintf (stdout, "-%c", short_opt);
    if (long_opt_non_empty) {
      fprintf (stdout, ", ");
    } else {
      fprintf (stdout, "  ");
    }
  }

  if (long_opt_non_empty) {
    fprintf (stdout, "--");
  } else {
    fprintf (stdout, "  ");
  }
  fprintf (stdout, "%-*s ",
	   LONG_OPT_LENGTH, long_opt);
  fprintf (stdout, "%-*s  ",
	   ARG_OPT_LENGTH, arg_opt);
  fprintf (stdout, "%s", str);
}

static inline void display_option_help_noopt(char* str) {
  display_option_help ('\0', "", "", str);
}

static inline void display_help(char *str){
  fprintf(stdout, "\nmsolve library for polynomial system solving, version %s\n", VERSION);
  fprintf(stdout, "implemented by J. Berthomieu, C. Eder, M. Safey El Din\n");
  fprintf(stdout, "\n");
  /* fprintf(stdout, "commit hash: %s\n\n", GIT_COMMIT_HASH); */

  fprintf(stdout, "Basic call:\n");
  fprintf(stdout, "\t ./msolve -f [FILE1] -o [FILE2]\n\n");
  fprintf(stdout, "FILE1 and FILE2 are respectively the input and output files\n\n");

  fprintf(stdout, "Standard options\n\n");
  display_option_help('f', "file", "FILE", "File name (mandatory).\n\n");
  display_option_help('h', "help", "", "Prints this help.\n");
  display_option_help('o', "output-file", "FILE",  "Name of output file.\n");
  display_option_help('t', "threads", "THR", "Number of threads to be used.\n");
  display_option_help_noopt("1       - one thread (default).\n");
  display_option_help_noopt("THR > 1 - THR threads.\n");
  display_option_help('v', "verbose", "VERB", "Level of verbosity, 0 - 2\n");
  display_option_help_noopt("0 - no output (default).\n");
  display_option_help_noopt("1 - global information at the start and\n");
  display_option_help_noopt("    end of the computation.\n");
  display_option_help_noopt("2 - detailed output for each step of the\n");
  display_option_help_noopt("    algorithm, e.g. matrix sizes, #pairs, ...\n");

  fprintf(stdout, "Input file format:\n");
  fprintf(stdout, "\t - first line: variables separated by a comma\n");
  fprintf(stdout, "\t                (no comma at end of line)\n");
  fprintf(stdout, "\t - second line: characteristic of the field\n");
  fprintf(stdout, "\t - next lines provide the polynomials (one per line),\n");
  fprintf(stdout, "\t   separated by a comma\n");
  fprintf(stdout, "\t   (no comma after the last polynomial)\n\n");
  fprintf(stdout, "Output file format: \n");
  fprintf(stdout, "When there is no solution in an algebraic closure of the base field\n[-1]:\n");
  fprintf(stdout, "Where there are infinitely many solutions in \nan algebraic closure of the base field: \n[1, nvars, -1,[]]:\n");
  fprintf(stdout,"Else:\n");
  fprintf(stdout,"Over prime fields: a rational parametrization of the solutions\n");
  fprintf(stdout,"When input coefficients are rational numbers: \nreal solutions to the input system (see the -P flag to recover a parametrization of the solutions)\n");
  fprintf(stdout, "See the msolve tutorial for more details (https://msolve.lip6.fr)\n");


  fprintf(stdout, "\nAdvanced options:\n\n");
  display_option_help('F', "", "FILE", "File name encoding parametrizations in binary format.\n\n");
  display_option_help('g', "groebner-basis", "GB", "Prints reduced Groebner bases of input system for\n");
  display_option_help_noopt("first prime characteristic w.r.t. grevlex ordering.\n");
  display_option_help_noopt("One element per line is printed, commata separated.\n");
  display_option_help_noopt("0 - Nothing is printed. (default)\n");
  display_option_help_noopt("1 - Leading ideal is printed.\n");
  display_option_help_noopt("2 - Full reduced Groebner basis is printed.\n");
  display_option_help('c',"", "GEN", "Handling genericity: If the staircase is not generic\n");
  display_option_help_noopt("enough, msolve can automatically try to fix this\n");
  display_option_help_noopt("situation via first trying a change of the order of\n");
  display_option_help_noopt("variables and finally adding a random linear form\n");
  display_option_help_noopt("with a new variable (smallest w.r.t. DRL)\n");
  display_option_help_noopt("0 - Nothing is done, msolve quits.\n");
  display_option_help_noopt("1 - Change order of variables.\n");
  display_option_help_noopt("2 - Change order of variables, then try adding a\n");
  display_option_help_noopt("    random linear form. (default)\n");
  display_option_help('d', "", "GEN", "Handling genericity further: If the staircase is not generic\n");
  display_option_help_noopt("enough, msolve can still try to perform the full computation\n");
  display_option_help_noopt("by computing some normal forms and build the multiplication matrix,\n");
  display_option_help_noopt("before fixing the situation via option -c\n");
  display_option_help_noopt("0 - No normal forms are computed.\n");
  display_option_help_noopt("1 - Few normal forms are computed.\n");
  display_option_help_noopt("2 - Some normal forms are computed. (default)\n");
  display_option_help_noopt("3 - Lots of normal forms are computed.\n");
  display_option_help_noopt("4 - All the normal forms are computed.\n");
  display_option_help('C', "", "", "Use sparse-FGLM-col algorithm:\n");
  display_option_help_noopt("Given an input file with k polynomials\n");
  display_option_help_noopt("compute the quotient of the ideal\n");
  display_option_help_noopt("generated by the first k-1 polynomials\n");
  display_option_help_noopt("with respect to the kth polynomial.\n");
  display_option_help('e', "elimination", "ELIM", "Define an elimination order: msolve supports two\n");
  display_option_help_noopt("blocks of variables, each block using the degree reverse\n");
  display_option_help_noopt("lexicographical monomial order. ELIM has to be a number\n");
  display_option_help_noopt("between 1 and #variables-1, and gives the number of\n");
  display_option_help_noopt("eliminated variables. The basis with the first block of\n");
  display_option_help_noopt("ELIM variables eliminated is then computed.\n");
  display_option_help('I', "isolate", "ISOL", "Isolates the real roots (provided some univariate data)\n");
  display_option_help_noopt("without re-computing a Gröbner basis\n");
  display_option_help_noopt("0 - no (default).\n");
  display_option_help_noopt("1 - yes.\n");
  display_option_help('l', "linear-algebra", "LIN", "Linear algebra variant to be applied:\n");
  display_option_help_noopt(" 1 - exact sparse / dense\n");
  display_option_help_noopt(" 2 - exact sparse (default)\n");
  display_option_help_noopt("42 - sparse / dense linearization (probabilistic)\n");
  display_option_help_noopt("44 - sparse linearization (probabilistic)\n");
  display_option_help('m', "", "MPR", "Maximal number of pairs used per matrix.\n");
  display_option_help_noopt("0 - unlimited (default).\n");
  display_option_help('n', "normal-form", "NF", "Given n input generators compute normal form of the last NF\n");
  display_option_help_noopt("elements of the input w.r.t. a degree reverse lexicographical\n");
  display_option_help_noopt("Gröbner basis of the first (n - NF) input elements.\n");
  display_option_help_noopt("At the moment this only works for prime field computations.\n");
  display_option_help_noopt("Combining this option with the \"-i\" option assumes that the\n");
  display_option_help_noopt("first (n - NF) elements generate already a degree reverse\n");
  display_option_help_noopt("lexicographical Gröbner basis.\n");
  display_option_help('p', "precision", "PRE", "Precision (in bits) on the output of real root isolation.\n");
  display_option_help_noopt("128 (default).\n");
  display_option_help('P', "parametrization", "PAR", "Get also rational parametrization of solution set.\n");
  display_option_help_noopt("0 (default). For a detailed description of the output\n");
  display_option_help_noopt("format please see the general output data format section\n");
  display_option_help_noopt("above.\n");
  display_option_help('L', "lifting-mulmat", "LIF", "Controls lifting of multplication matrices over the rationals.\n");
  display_option_help_noopt("0 - no lifting (default). \n");
  display_option_help_noopt("1 - matrices are lifted.\n");
  display_option_help_noopt("Warning: when activated, this option may cause higher memory consumption.\n");
  display_option_help('q', "", "Q", "Uses signature-based algorithms.\n");
  display_option_help_noopt("0 - no (default).\n");
  display_option_help_noopt("1 - yes.\n");
  display_option_help('\0', "random_seed", "SEED", "Random seed to initialize the pseudo\n");
  display_option_help_noopt("random generator\n");
  display_option_help_noopt("0        - time(0) will be used (default)\n");
  display_option_help_noopt("SEED > 0 - use at your own risks;\n");
  display_option_help_noopt("           this is intended for developers and\n");
  display_option_help_noopt("           debug purposes only\n");
  display_option_help('r', "reduce-gb", "RED", "Reduce Groebner basis.\n");
  display_option_help_noopt("0 - no.\n");
  display_option_help_noopt("1 - yes (default).\n");
  /* display_option_help('R', "", "REF", "Refinement fo real roots.\n"); */
  /* display_option_help_noopt("(not implemented yet).\n"); */
  display_option_help('s', "", "HTS", "Initial hash table size given\n");
  display_option_help_noopt("as power of two.\n");
  display_option_help_noopt("17 (default).\n");
  display_option_help('S', "", "", "Use f4sat saturation algorithm:\n");
  display_option_help_noopt("Given an input file with k polynomials\n");
  display_option_help_noopt("compute the saturation of the ideal\n");
  display_option_help_noopt("generated by the first k-1 polynomials\n");
  display_option_help_noopt("with respect to the kth polynomial.\n");
  display_option_help_noopt("Note: At the moment restricted to 32 bit\n");
  display_option_help_noopt("prime fields.\n");
  display_option_help('u', "", "UHT", "Number of steps after which the\n");
  display_option_help_noopt("hash table is newly generated.\n");
  display_option_help_noopt("0 - no update (default).\n");
  display_option_help('V', "version", "", "Prints msolve's version\n");
}

static void getoptions(
        int argc,
        char **argv,
        int32_t *initial_hts,
        int32_t *nthreads,
        int32_t *max_pairs,
        int32_t *elim_block_len,
        int32_t *linear_algebra,
        int32_t *use_signatures,
        int32_t *update_ht,
        int32_t *reduce_gb,
        int32_t *print_gb,
        int32_t *truncate_lifting,
        int32_t *genericity_handling,
        int32_t *unstable_staircase,
        int32_t *saturate,
        int32_t *colon,
        int32_t *normal_form,
        int32_t *normal_form_matrix,
        int32_t *is_gb,
        int32_t *lift_matrix,
        int32_t *get_param,
        int32_t *precision,
        int32_t *refine,
        int32_t *isolate,
        int32_t *generate_pbm_files,
	int32_t *seed,
        int32_t *info_level,
        files_gb *files){
  int opt, errflag = 0, fflag = 1;
  char *filename = NULL;
  char *bin_filename = NULL;
  char *out_fname = NULL;
  char *bin_out_fname = NULL;
  opterr = 1;
  /* char short_options[] = "hf:N:F:v:l:t:e:o:O:u:iI:p:P:L:q:g:c:s:SCr:R:m:M:n:d:Vf:"; */
  char short_options[] = "c:Cd:e:f:F:g:hiI:l:L:m:M:n:N:o:O:p:P:q:r:R:s:St:u:v:V";

  int random_seed_flag = 0;
  struct option long_options[] = {
    {"elimination", required_argument, NULL, 'e'},
    {"file", required_argument, NULL, 'f'},
    {"groebner-basis", required_argument, NULL, 'g'},
    {"help", no_argument, NULL, 'h'},
    {"isolate", required_argument, NULL, 'I'},
    {"linear-algebra", required_argument, NULL, 'l'},
    {"lifting-mulmat", required_argument, NULL, 'L'},
    {"normal-form", required_argument, NULL, 'n'},
    {"output-file", required_argument, NULL, 'o'},
    {"precision", required_argument, NULL, 'p'},
    {"parametrization", required_argument, NULL, 'P'},
    {"random-seed", required_argument, &random_seed_flag, 1},
    {"reduce-gb", required_argument, NULL, 'r'},
    {"threads", required_argument, NULL, 't'},
    {"verbose", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, 'V'},
    {NULL,0,NULL,0}
  };
  
  while(1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    opt = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (opt == -1) {
      /* processed all command-line options */
      break;
    }
    
    switch(opt) {
    case 0: /* no short option equivalent */
      if (random_seed_flag == 1) {
	*seed = strtol(optarg, NULL, 10);
      }
      break;
    case 'N':
      *truncate_lifting = strtol(optarg, NULL, 10);
      break;
    case 'h':
      display_help(argv[0]);
      exit(0);
    case 'V':
      fprintf(stdout, "%s\n", VERSION);
      exit(0);
    case 'e':
      *elim_block_len = strtol(optarg, NULL, 10);
      if (*elim_block_len < 0) {
          *elim_block_len = 0;
      }
      break;
    case 'u':
      *update_ht = strtol(optarg, NULL, 10);
      break;
    case 'p':
      *precision = strtol(optarg, NULL, 10);
      if (*precision < 0) {
          *precision = 64;
      }
      /* if (*precision > 100) { */
      /*     *precision = 100; */
      /* } */
      break;
    case 'q':
      *use_signatures = strtol(optarg, NULL, 10);
      if (*use_signatures < 0) {
          *use_signatures = 0;
      }
      if (*use_signatures > 3) {
          *use_signatures = 0;
      }
      break;
    case 'R':
      *refine = 1;
      break;
    case 'i':
      *is_gb = 1;
      break;
    case 'I':
      *isolate = strtol(optarg, NULL, 10);
      break;
    case 's':
      *initial_hts = strtol(optarg, NULL, 10);
      break;
    case 'S':
      *saturate = 1;
      break;
    case 'C':
      *colon = 1;
      break;
    case 'r':
      *reduce_gb = strtol(optarg, NULL, 10);
      break;
    case 'm':
      *max_pairs = strtol(optarg, NULL, 10);
      break;
    case 'v':
      *info_level = strtol(optarg, NULL, 10);
      break;
    case 'L':
      *lift_matrix= strtol(optarg, NULL, 10);
      break;
    case 't':
      *nthreads = strtol(optarg, NULL, 10);
      break;
    case 'l':
      *linear_algebra = strtol(optarg, NULL, 10);
      break;
    case 'f':
      fflag = 0;
      filename = optarg;
      break;
    case 'F':
      fflag = 0;
      bin_filename = optarg;
      break;
    case 'o':
      out_fname = optarg;
      break;
    case 'O':
      bin_out_fname = optarg;
      break;
    case 'P':
      *get_param = strtol(optarg, NULL, 10);
      if (*get_param <= 0) {
          *get_param = 0;
      }
      /* if (*get_param > 1) { */
      /*     *get_param = 1; */
      /* } */
      break;
    case 'g':
      *print_gb = strtol(optarg, NULL, 10);
      if (*print_gb < 0) {
          *print_gb = 0;
      }
      if (*print_gb > 2) {
          *print_gb = 2;
      }
      break;
    case 'c':
      *genericity_handling = strtol(optarg, NULL, 10);
      if (*genericity_handling < 0) {
          *genericity_handling = 0;
      }
      if (*genericity_handling > 2) {
          *genericity_handling = 2;
      }
      break;
    case 'd':
      *unstable_staircase = strtol(optarg, NULL, 10);
      if (*unstable_staircase < 0) {
          *unstable_staircase = 0;
      }
      if (*unstable_staircase > 4) {
          *unstable_staircase = 4;
      }
      break;
    case 'n':
      *normal_form = strtol(optarg, NULL, 10);
      if (*normal_form < 0) {
          *normal_form  = 0;
      }
      break;
    case 'M':
      *normal_form_matrix = strtol(optarg, NULL, 10);
      if (*normal_form_matrix < 0) {
          *normal_form_matrix  = 0;
      }
      break;
    default:
      errflag++;
      break;
    }
  }
  if(fflag){
    fprintf(stderr,"No given file\n");
    display_help(argv[0]);
    exit(1);
  }
  if(errflag){
    fprintf(stderr, "Invalid usage\n");
    display_help(argv[0]);
    exit(1);
  }
  files->in_file = filename;
  files->bin_file = bin_filename;
  files->out_file = out_fname;
  files->bin_out_file = bin_out_fname;
}


int main(int argc, char **argv){

    /* timinigs */
    double st0 = cputime();
    double rt0 = realtime();

    /**
      We get values from the command line.
     **/
    int32_t la_option             = 2; // by default
    int32_t use_signatures        = 0;
    int32_t nr_threads            = 1;
    int32_t info_level            = 0;
    int32_t initial_hts           = 17;
    int32_t max_pairs             = 0;
    int32_t elim_block_len        = 0;
    int32_t update_ht             = 0;
    int32_t generate_pbm          = 0;
    int32_t reduce_gb             = 1;
    int32_t print_gb              = 0;
    int32_t truncate_lifting      = 0;
    int32_t genericity_handling   = 2;
    int32_t unstable_staircase    = 2;
    int32_t saturate              = 0;
    int32_t colon                 = 0;
    int32_t normal_form           = 0;
    int32_t normal_form_matrix    = 0;
    int32_t is_gb                 = 0;
    int32_t lift_matrix           = 0;
    int32_t get_param             = 0;
    int32_t precision             = 64;
    int32_t refine                = 0; /* not used at the moment */
    int32_t isolate               = 0; /* not used at the moment */
    int32_t seed                  = 0;

    files_gb *files = malloc(sizeof(files_gb));
    if(files == NULL) exit(1);
    files->in_file = NULL;
    files->bin_file = NULL;
    files->out_file = NULL;
    files->bin_out_file = NULL;
    getoptions(argc, argv, &initial_hts, &nr_threads, &max_pairs,
               &elim_block_len, &la_option, &use_signatures, &update_ht,
               &reduce_gb, &print_gb, &truncate_lifting, &genericity_handling, &unstable_staircase, &saturate, &colon,
               &normal_form, &normal_form_matrix, &is_gb, &lift_matrix, &get_param,
               &precision, &refine, &isolate, &generate_pbm,
	       &seed, &info_level, files);

    /* srand initialization */
    if (seed == 0) {
      seed = time(0);
    }
    srand(seed);
    if (info_level) {
      fprintf (stdout,"Initial seed for pseudo-random number generator ");
      fprintf (stdout,"is %d\n",seed);
    }

    FILE *fh  = fopen(files->in_file, "r");
    FILE *bfh  = fopen(files->bin_file, "r");

    if (fh == NULL && bfh == NULL) {
      fprintf(stderr, "Input file not found.\n");
      exit(1);
    }
    if(fh!=NULL){
      fclose(fh);
    }
    if(bfh != NULL){
      fclose(bfh);
    }
    fh =  NULL;
    bfh =  NULL;

    /* clear out_file if given */
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "w");
      if(ofile == NULL){
        fprintf(stderr, "Cannot open output file\n");
        exit(1);
      }
      fclose(ofile);
    }
    /**
       We get from files the requested data.
    **/
    //  int32_t mon_order   = 0;
    int32_t nr_vars     = 0;
    int32_t field_char  = 9001;
    int32_t nr_gens     = 0;
    data_gens_ff_t *gens = allocate_data_gens();

    get_data_from_file(files->in_file, &nr_vars, &field_char, &nr_gens, gens);
#ifdef IODEBUG
    display_gens(stdout, gens);
#endif

    gens->rand_linear           = 0;
    gens->random_linear_form = malloc(sizeof(int32_t)*(nr_vars));
    gens->elim = elim_block_len;

    if(0 < field_char && field_char < pow(2, 15) && la_option > 2){
        if(info_level){
            fprintf(stderr, "Warning: characteristic is too low for choosing \nprobabilistic linear algebra\n");
            fprintf(stderr, "\t linear algebra option set to 2\n");
        }
        la_option = 2;
    }

    /* data structures for parametrization */
    param_t *param  = NULL; //malloc(sizeof(param_t));
    mpz_param_t *mpz_paramp = malloc(sizeof(mpz_param_t));
    mpz_param_init(*mpz_paramp);

    long nb_real_roots      = 0;
    interval *real_roots    = NULL;
    real_point_t *real_pts  = NULL;

    /* main msolve functionality */
    int ret = core_msolve(la_option, use_signatures, nr_threads, info_level,
                          initial_hts, max_pairs, elim_block_len, update_ht,
                          generate_pbm, reduce_gb, print_gb, truncate_lifting, get_param,
                          genericity_handling, unstable_staircase, saturate, colon, normal_form,
                          normal_form_matrix, is_gb, lift_matrix, precision,
                          files, gens,
            &param, mpz_paramp, &nb_real_roots, &real_roots, &real_pts);

    /* free parametrization */
    if(param != NULL && gens->field_char){
        free_fglm_param(param);
    }
    mpz_param_clear(*mpz_paramp);
    free(mpz_paramp);

    if (nb_real_roots > 0) {
        for(long i = 0; i < nb_real_roots; i++){
          real_point_clear(real_pts[i]);
          mpz_clear(real_roots[i].numer);
        }
        free(real_pts);
    }
    free(real_roots);

    /* timings */
    if (info_level > 0) {
        double st1 = cputime();
        double rt1 = realtime();
        fprintf(stdout, "\n\n-------------------------------------------------\
-----------------------------------\n");
        fprintf(stdout, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(stdout, "-------------------------------------------------\
-----------------------------------\n");
    }
    free_data_gens(gens);

    free(files);
    return ret;
}
