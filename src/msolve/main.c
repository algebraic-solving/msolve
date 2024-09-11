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

static inline void display_help(char *str){
  fprintf(stdout, "\nmsolve library for polynomial system solving, version %s\n", VERSION);
  fprintf(stdout, "implemented by J. Berthomieu, C. Eder, M. Safey El Din\n");
  fprintf(stdout, "\n");
  /* fprintf(stdout, "commit hash: %s\n\n", GIT_COMMIT_HASH); */

  fprintf(stdout, "Basic call:\n");
  fprintf(stdout, "\t ./msolve -f [FILE1] -o [FILE2]\n\n");
  fprintf(stdout, "FILE1 and FILE2 are respectively the input and output files\n\n");

  fprintf(stdout, "Standard options\n\n");
  fprintf(stdout, "-f FILE  File name (mandatory).\n\n");
  fprintf(stdout, "-h       Prints this help.\n");
  fprintf(stdout, "-o FILE  Name of output file.\n");
  fprintf(stdout, "-t THR   Number of threads to be used.\n");
  fprintf(stdout, "         Default: 1.\n");
  fprintf(stdout, "-v n     Level of verbosity, 0 - 2\n");
  fprintf(stdout, "         0 - no output (default).\n");
  fprintf(stdout, "         1 - global information at the start and\n");
  fprintf(stdout, "             end of the computation.\n");
  fprintf(stdout, "         2 - detailed output for each step of the\n");
  fprintf(stdout, "             algorithm, e.g. matrix sizes, #pairs, ...\n");

  fprintf(stdout, "Input file format:\n");
  fprintf(stdout, "\t - first line: variables separated by a comma\n");
  fprintf(stdout, "\t - second line: characteristic of the field\n");
  fprintf(stdout, "\t - next lines provide the polynomials (one per line),\n");
  fprintf(stdout, "\t   separated by a comma\n");
  fprintf(stdout, "\t   (no comma after the final polynomial)\n\n");
  fprintf(stdout, "Output file format: \n");
  fprintf(stdout, "When there is no solution in an algebraic closure of the base field\n[-1]:\n");
  fprintf(stdout, "Where there are infinitely many solutions in \nan algebraic closure of the base field: \n[1, nvars, -1,[]]:\n");
  fprintf(stdout,"Else:\n");
  fprintf(stdout,"Over prime fields: a rational parametrization of the solutions\n");
  fprintf(stdout,"When input coefficients are rational numbers: \nreal solutions to the input system (see the -P flag to recover a parametrization of the solutions)\n");
  fprintf(stdout, "See the msolve tutorial for more details (https://msolve.lip6.fr)\n");


  fprintf(stdout, "\nAdvanced options:\n\n");
  fprintf(stdout, "-F FILE  File name encoding parametrizations in binary format.\n\n");
  fprintf(stdout, "-g GB    Prints reduced Groebner bases of input system for\n");
  fprintf(stdout, "         first prime characteristic w.r.t. grevlex ordering.\n");
  fprintf(stdout, "         One element per line is printed, commata separated.\n");
  fprintf(stdout, "         0 - Nothing is printed. (default)\n");
  fprintf(stdout, "         1 - Leading ideal is printed.\n");
  fprintf(stdout, "         2 - Full reduced Groebner basis is printed.\n");
  fprintf(stdout, "-c GEN   Handling genericity: If the staircase is not generic\n");
  fprintf(stdout, "         enough, msolve can automatically try to fix this\n");
  fprintf(stdout, "         situation via first trying a change of the order of\n");
  fprintf(stdout, "         variables and finally adding a random linear form\n");
  fprintf(stdout, "         with a new variable (smallest w.r.t. DRL)\n");
  fprintf(stdout, "         0 - Nothing is done, msolve quits.\n");
  fprintf(stdout, "         1 - Change order of variables.\n");
  fprintf(stdout, "         2 - Change order of variables, then try adding a\n");
  fprintf(stdout, "             random linear form. (default)\n");
  fprintf(stdout, "-d GEN   Handling genericity further: If the staircase is not generic\n");
  fprintf(stdout, "         enough, msolve can still try to perform the full computation\n");
  fprintf(stdout, "         by computing some normal forms and build the multiplication matrix,\n");
  fprintf(stdout, "         before fixing the situation via option -c\n");
  fprintf(stdout, "         0 - No normal forms are computed.\n");
  fprintf(stdout, "         1 - Few normal forms are computed.\n");
  fprintf(stdout, "         2 - Some normal forms are computed. (default)\n");
  fprintf(stdout, "         3 - Lots of normal forms are computed.\n");
  fprintf(stdout, "         4 - All the normal forms are computed.\n");
  fprintf(stdout, "-C       Use sparse-FGLM-col algorithm:\n");
  fprintf(stdout, "         Given an input file with k polynomials\n");
  fprintf(stdout, "         compute the quotient of the ideal\n");
  fprintf(stdout, "         generated by the first k-1 polynomials\n");
  fprintf(stdout, "         with respect to the kth polynomial.\n");
  fprintf(stdout, "-e ELIM  Define an elimination order: msolve supports two\n");
  fprintf(stdout, "         blocks, each block using degree reverse lexicographical\n");
  fprintf(stdout, "         monomial order. ELIM has to be a number between\n");
  fprintf(stdout, "         1 and #variables-1. The basis the first block eliminated\n");
  fprintf(stdout, "         is then computed.\n");
  fprintf(stdout, "-I       Isolates the real roots (provided some univariate data)\n");
  fprintf(stdout, "         without re-computing a Gröbner basis\n");
  fprintf(stdout, "         Default: 0 (no).\n");
  fprintf(stdout, "-l LIN   Linear algebra variant to be applied:\n");
  fprintf(stdout, "          1 - exact sparse / dense\n");
  fprintf(stdout, "          2 - exact sparse (default)\n");
  fprintf(stdout, "         42 - sparse / dense linearization (probabilistic)\n");
  fprintf(stdout, "         44 - sparse linearization (probabilistic)\n");
  fprintf(stdout, "-m MPR   Maximal number of pairs used per matrix.\n");
  fprintf(stdout, "         Default: 0 (unlimited).\n");
  fprintf(stdout, "-n NF    Given n input generators compute normal form of the last NF\n");
  fprintf(stdout, "         elements of the input w.r.t. a degree reverse lexicographical\n");
  fprintf(stdout, "         Gröbner basis of the first (n - NF) input elements.\n");
  fprintf(stdout, "         At the moment this only works for prime field computations.\n");
  fprintf(stdout, "         Combining this option with the \"-i\" option assumes that the\n");
  fprintf(stdout, "         first (n - NF) elements generate already a degree reverse\n");
  fprintf(stdout, "         lexicographical Gröbner basis.\n");
  fprintf(stdout, "-p PRE   Precision of the real root isolation.\n");
  fprintf(stdout, "         Default is 32.\n");
  fprintf(stdout, "-P PAR   Get also rational parametrization of solution set.\n");
  fprintf(stdout, "         Default is 0. For a detailed description of the output\n");
  fprintf(stdout, "         format please see the general output data format section\n");
  fprintf(stdout, "         above.\n");
  fprintf(stdout, "-L LIF   Controls lifting of multplication matrices over the rationals.\n");
  fprintf(stdout, "         Default is 0 (no lifting). \n");
  fprintf(stdout, "         Matrices are lifted when LIF is 1.\n");
  fprintf(stdout, "         Warning: when activated, this option may cause higher memory consumption.\n");
  fprintf(stdout, "-q Q     Uses signature-based algorithms.\n");
  fprintf(stdout, "         Default: 0 (no).\n");
  fprintf(stdout, "-r RED   Reduce Groebner basis.\n");
  fprintf(stdout, "         Default: 1 (yes).\n");
  /* fprintf(stdout, "-R       Refinement fo real roots.\n"); */
  /* fprintf(stdout, "         (not implemented yet).\n"); */
  fprintf(stdout, "-s HTS   Initial hash table size given\n");
  fprintf(stdout, "         as power of two. Default: 17.\n");
  fprintf(stdout, "-S       Use f4sat saturation algorithm:\n");
  fprintf(stdout, "         Given an input file with k polynomials\n");
  fprintf(stdout, "         compute the saturation of the ideal\n");
  fprintf(stdout, "         generated by the first k-1 polynomials\n");
  fprintf(stdout, "         with respect to the kth polynomial.\n");
  fprintf(stdout, "         Note: At the moment restricted to 32 bit\n");
  fprintf(stdout, "         prime fields.\n");
  fprintf(stdout, "-u UHT   Number of steps after which the\n");
  fprintf(stdout, "         hash table is newly generated.\n");
  fprintf(stdout, "         Default: 0, i.e. no update.\n");
  fprintf(stdout, "-V       Prints msolve's version\n");
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
        int32_t *info_level,
        files_gb *files){
  int opt, errflag = 0, fflag = 1;
  char *filename = NULL;
  char *bin_filename = NULL;
  char *out_fname = NULL;
  char *bin_out_fname = NULL;
  opterr = 1;
  char options[] = "hf:N:F:v:l:t:e:o:O:u:iI:p:P:L:q:g:c:s:SCr:R:m:M:n:d:Vf:";
  while((opt = getopt(argc, argv, options)) != -1) {
    switch(opt) {
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
          *precision = 128;
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
    int32_t precision             = 128;
    int32_t refine                = 0; /* not used at the moment */
    int32_t isolate               = 0; /* not used at the moment */

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
               &precision, &refine, &isolate, &generate_pbm, &info_level, files);

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
    param_t *param  = NULL;
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
    free(param);
    mpz_param_clear(*mpz_paramp);

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
        fprintf(stderr, "\n\n-------------------------------------------------\
-----------------------------------\n");
        fprintf(stderr, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(stderr, "-------------------------------------------------\
-----------------------------------\n");
    }
    free_data_gens(gens);

    free(files);
    return ret;
}
