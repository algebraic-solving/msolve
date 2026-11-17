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

#define LONG_OPT_LENGTH 18
#define ARG_OPT_LENGTH 4


static inline void display_option_help(char short_opt, char *long_opt,
				       char *arg_opt, char* str) {
  int long_opt_non_empty= strcmp (long_opt, "");

  if (short_opt == 0) {
    fprintf (OUTSTREAM, "    ");
  } else {
    fprintf (OUTSTREAM, "-%c", short_opt);
    if (long_opt_non_empty) {
      fprintf (OUTSTREAM, ", ");
    } else {
      fprintf (OUTSTREAM, "  ");
    }
  }

  if (long_opt_non_empty) {
    fprintf (OUTSTREAM, "--");
  } else {
    fprintf (OUTSTREAM, "  ");
  }
  fprintf (OUTSTREAM, "%-*s ",
	   LONG_OPT_LENGTH, long_opt);
  fprintf (OUTSTREAM, "%-*s  ",
	   ARG_OPT_LENGTH, arg_opt);
  fprintf (OUTSTREAM, "%s", str);
}

static inline void display_option_help_noopt(char* str) {
  display_option_help (0, "", "", str);
}

static inline void display_help(char *str){
  fprintf(OUTSTREAM, "\nmsolve library for polynomial system solving, version %s\n", VERSION);
  fprintf(OUTSTREAM, "implemented by J. Berthomieu, C. Eder, M. Safey El Din\n");
  /* fprintf(OUTSTREAM, "commit hash: %s\n\n", GIT_COMMIT_HASH); */
  fprintf(OUTSTREAM, "\nBasic call:\n\n");
  fprintf(OUTSTREAM, "  %s -f [FILE1] -o [FILE2]\n\n",
          str != NULL ? str : "msolve");
  fprintf(OUTSTREAM, "  FILE1 and FILE2 are respectively the input and output files.\n");
  fprintf(OUTSTREAM, "\nInput file format:\n\n");
  fprintf(OUTSTREAM, "  - first line: variables separated by a comma (no comma at end of line)\n");
  fprintf(OUTSTREAM, "  - second line: characteristic of the field, or 2^4 resp. 2^8 for\n");
  fprintf(OUTSTREAM, "    GF(16) = GF(2)[t]/(t^4+t+1) resp. GF(256) = GF(2)[t]/(t^8+t^4+t^3+t+1)\n");
  fprintf(OUTSTREAM, "    (only option -g), their elements are written as bytes in hexadecimal\n");
  fprintf(OUTSTREAM, "    notation, e.g. 0x1b for t^4+t^3+t+1, they have to be reduced, i.e. less\n");
  fprintf(OUTSTREAM, "    than 0x10 resp. 0x100, integers are taken modulo 2, or boolean for the\n");
  fprintf(OUTSTREAM, "    boolean ring over GF(2), i.e. x^2 = x for all variables x (only option -g)\n");
  fprintf(OUTSTREAM, "  - next lines provide the polynomials (one per line), separated by a comma (no\n");
  fprintf(OUTSTREAM, "    comma after the last polynomial)\n");
  fprintf(OUTSTREAM, "\nDefault output file format (changed by some command line options, see below):\n\n");
  fprintf(OUTSTREAM, "  - When there is no solution in an algebraic closure of the base field:\n");
  fprintf(OUTSTREAM, "      [-1]:\n");
  fprintf(OUTSTREAM, "  - When there are infinitely many solutions in an algebraic closure of the\n");
  fprintf(OUTSTREAM, "    base field:\n");
  fprintf(OUTSTREAM, "      [1, nvars, -1, []]:\n");
  fprintf(OUTSTREAM, "  - Else, over prime fields: a rational parametrization of the solutions. When\n");
  fprintf(OUTSTREAM, "    input coefficients are rational numbers: real solutions to the input system\n");
  fprintf(OUTSTREAM, "    (see the -P flag to recover a parametrization of the solutions).\n");
  fprintf(OUTSTREAM, "\n");
  fprintf(OUTSTREAM, "  See the msolve tutorial for more details (https://msolve.lip6.fr).\n");
  fprintf(OUTSTREAM, "\nStandard options:\n\n");
  display_option_help('f', "file", "FILE", "File name (mandatory).\n");
  display_option_help('o', "output-file", "FILE", "Name of output file.\n");
  display_option_help('t', "threads", "THR", "Number of threads to be used.\n");
  display_option_help_noopt("1       - one thread (default).\n");
  display_option_help_noopt("THR > 1 - THR threads.\n");
  display_option_help('v', "verbose", "VERB", "Level of verbosity, 0 - 2.\n");
  display_option_help_noopt("0 - no output (default).\n");
  display_option_help_noopt("1 - global information at the start and end of\n");
  display_option_help_noopt("    the computation.\n");
  display_option_help_noopt("2 - detailed output for each step of the\n");
  display_option_help_noopt("    algorithm, e.g. matrix sizes, #pairs, ...\n");
  fprintf(OUTSTREAM, "\nAdvanced options affecting the output:\n\n");
  display_option_help('P', "parametrization", "PAR", "Also outputs the rational parametrization of the\n");
  display_option_help_noopt("solution set; only affects systems with rational\n");
  display_option_help_noopt("coefficients.\n");
  display_option_help_noopt("0 - output is as above (default).\n");
  display_option_help_noopt("1 - outputs the rational parametrization in\n");
  display_option_help_noopt("    addition to the real solutions.\n");
  display_option_help_noopt("2 - outputs only the rational parametrization.\n");
  display_option_help('p', "precision", "PRE", "Precision (in bits) on the output of real root\n");
  display_option_help_noopt("isolation.\n");
  display_option_help_noopt("128 (default).\n");
  display_option_help('g', "groebner-basis", "GB", "Outputs information about the reduced Groebner\n");
  display_option_help_noopt("basis of the input system w.r.t. the grevlex\n");
  display_option_help_noopt("ordering. One element per line is printed,\n");
  display_option_help_noopt("comma-separated.\n");
  display_option_help_noopt("0 - output is as above (default).\n");
  display_option_help_noopt("1 - output is the leading ideal of the input\n");
  display_option_help_noopt("    system.\n");
  display_option_help_noopt("2 - output is the full reduced Groebner basis of\n");
  display_option_help_noopt("    the input system.\n");
  display_option_help('e', "elimination", "ELIM", "Define an elimination order: msolve supports two\n");
  display_option_help_noopt("blocks of variables, each block using the degree\n");
  display_option_help_noopt("reverse lexicographical monomial order. ELIM has\n");
  display_option_help_noopt("to be a number between 1 and #variables-1, and\n");
  display_option_help_noopt("gives the number of eliminated variables. The\n");
  display_option_help_noopt("output is then the DRL Groebner basis of the\n");
  display_option_help_noopt("input system with the first block of ELIM\n");
  display_option_help_noopt("variables eliminated.\n");
  display_option_help('r', "reduce-gb", "RED", "Reduce the Groebner basis computed with the \"-g\"\n");
  display_option_help_noopt("or \"-e\" option.\n");
  display_option_help_noopt("0 - no.\n");
  display_option_help_noopt("1 - yes (default).\n");
  display_option_help('n', "normal-form", "NF", "Given n input generators, outputs the normal\n");
  display_option_help_noopt("forms of the last NF elements of the input\n");
  display_option_help_noopt("w.r.t. a degree reverse lexicographical Groebner\n");
  display_option_help_noopt("basis of the first (n - NF) input elements. At\n");
  display_option_help_noopt("the moment this only works for prime field\n");
  display_option_help_noopt("computations. Combining this option with the\n");
  display_option_help_noopt("\"-i\" option assumes that the first (n - NF)\n");
  display_option_help_noopt("elements already form a degree reverse\n");
  display_option_help_noopt("lexicographical Groebner basis.\n");
  display_option_help('C', "", "", "Use the sparse-FGLM-col algorithm: given an\n");
  display_option_help_noopt("input file with k polynomials, outputs the\n");
  display_option_help_noopt("quotient of the ideal generated by the first k-1\n");
  display_option_help_noopt("polynomials with respect to the kth polynomial.\n");
  display_option_help('S', "", "", "Use the f4sat saturation algorithm: given an\n");
  display_option_help_noopt("input file with k polynomials, outputs the\n");
  display_option_help_noopt("saturation of the ideal generated by the first\n");
  display_option_help_noopt("k-1 polynomials with respect to the kth\n");
  display_option_help_noopt("polynomial.\n");
  display_option_help_noopt("Note: at the moment restricted to 32-bit prime\n");
  display_option_help_noopt("fields.\n");
  display_option_help(0, "boolean", "", "Computes in the boolean ring over GF(2), i.e.\n");
  display_option_help_noopt("with x^2 = x for all variables x. The field\n");
  display_option_help_noopt("characteristic of the input file has to be 2.\n");
  display_option_help_noopt("Only works with option -g set to 1 or 2.\n");
  display_option_help(0, "boolean-zeros", "", "Outputs the Groebner basis in the boolean ring\n");
  display_option_help_noopt("via the common zeros of the input polynomials\n");
  display_option_help_noopt("using the Buchberger-Moeller algorithm. All 2^n\n");
  display_option_help_noopt("points are evaluated, so this is meant for few\n");
  display_option_help_noopt("variables (at most 32). Implies --boolean. Only\n");
  display_option_help_noopt("works with option -g set to 1 or 2.\n");
  display_option_help(0, "lex", "", "Outputs the Groebner basis w.r.t. the\n");
  display_option_help_noopt("lexicographical monomial order. Only works with\n");
  display_option_help_noopt("option --boolean-zeros.\n");
  fprintf(OUTSTREAM, "\nAdvanced options affecting the computation:\n\n");
  display_option_help('c', "", "GEN", "Handling genericity: if the staircase is not\n");
  display_option_help_noopt("generic enough, msolve can automatically try to\n");
  display_option_help_noopt("fix this situation via first trying a change of\n");
  display_option_help_noopt("the order of variables and finally adding a\n");
  display_option_help_noopt("random linear form with a new variable (smallest\n");
  display_option_help_noopt("w.r.t. DRL).\n");
  display_option_help_noopt("0 - nothing is done, msolve quits.\n");
  display_option_help_noopt("1 - change order of variables.\n");
  display_option_help_noopt("2 - change order of variables, then try adding a\n");
  display_option_help_noopt("    random linear form (default).\n");
  display_option_help('d', "", "GEN", "Handling genericity further: if the staircase is\n");
  display_option_help_noopt("not generic enough, msolve can still try to\n");
  display_option_help_noopt("perform the full computation by computing some\n");
  display_option_help_noopt("normal forms and building the multiplication\n");
  display_option_help_noopt("matrix, before fixing the situation via\n");
  display_option_help_noopt("option -c.\n");
  display_option_help_noopt("0 - no normal forms are computed.\n");
  display_option_help_noopt("1 - few normal forms are computed.\n");
  display_option_help_noopt("2 - some normal forms are computed (default).\n");
  display_option_help_noopt("3 - lots of normal forms are computed.\n");
  display_option_help_noopt("4 - all the normal forms are computed.\n");
  display_option_help(0, "dense-route", "DEN", "With a separating linear form and an extra\n");
  display_option_help_noopt("variable, over the rationals: route of the\n");
  display_option_help_noopt("multi-modular steps.\n");
  display_option_help_noopt("0 - chosen from the sizes of the F4 matrices\n");
  display_option_help_noopt("    (default).\n");
  display_option_help_noopt("1 - on the system with the extra variable.\n");
  display_option_help_noopt("2 - on the original system, parametrizations\n");
  display_option_help_noopt("    from its multiplication matrix.\n");
  display_option_help(0, "elim-drl", "EDRL", "With option -e: compute a Groebner basis\n");
  display_option_help_noopt("w.r.t. the degree reverse lexicographical order\n");
  display_option_help_noopt("first and the one w.r.t. the elimination order\n");
  display_option_help_noopt("from it.\n");
  display_option_help_noopt("0 - never.\n");
  display_option_help_noopt("1 - always.\n");
  display_option_help_noopt("2 - if there are at least as many equations as\n");
  display_option_help_noopt("    variables and more than one variable is\n");
  display_option_help_noopt("    eliminated (default).\n");
  display_option_help('i', "", "", "Assume that the first input generators already\n");
  display_option_help_noopt("form a degree reverse lexicographical Groebner\n");
  display_option_help_noopt("basis and skip its computation. Only with -n\n");
  display_option_help_noopt("or -C over prime fields; the result is wrong if\n");
  display_option_help_noopt("the assumption does not hold.\n");
  display_option_help('l', "linear-algebra", "LIN", "Linear algebra variant to be applied:\n");
  display_option_help_noopt(" 1 - exact sparse / dense.\n");
  display_option_help_noopt(" 2 - exact sparse (default).\n");
  display_option_help_noopt("42 - sparse / dense linearization\n");
  display_option_help_noopt("     (probabilistic).\n");
  display_option_help_noopt("44 - sparse linearization (probabilistic).\n");
  display_option_help('L', "lifting-mulmat", "LIF", "Controls lifting of multiplication matrices over\n");
  display_option_help_noopt("the rationals.\n");
  display_option_help_noopt("0 - no lifting (default).\n");
  display_option_help_noopt("1 - matrices are lifted.\n");
  display_option_help_noopt("Warning: when activated, this option may cause\n");
  display_option_help_noopt("higher memory consumption.\n");
  display_option_help(0, "lifting-param", "LIFT", "Lifting strategy for the parametrizations over\n");
  display_option_help_noopt("the rationals.\n");
  display_option_help_noopt("0 - p-adic or multi-modular, chosen by a cost\n");
  display_option_help_noopt("    model (default); multi-modular with the\n");
  display_option_help_noopt("    probabilistic linear algebra (-l 42, 44).\n");
  display_option_help_noopt("1 - multi-modular.\n");
  display_option_help_noopt("2 - p-adic when applicable.\n");
  display_option_help_noopt("By default, msolve sets the parameters of the\n");
  display_option_help_noopt("p-adic lifting automatically; the options below\n");
  display_option_help_noopt("override them and are intended for experiments.\n");
  display_option_help(0, "hensel-precision", "N", "Precision p^N per prime.\n");
  display_option_help_noopt("0 - chosen by the cost model (default).\n");
  display_option_help(0, "hensel-cost-scale", "K", "Factor on the cost of a modular run in the cost\n");
  display_option_help_noopt("model; above 1, the p-adic lifting is chosen\n");
  display_option_help_noopt("more often.\n");
  display_option_help_noopt("1 (default).\n");
  display_option_help(0, "hensel-memory", "MB", "Memory limit of the p-adic lifting, in\n");
  display_option_help_noopt("megabytes, which bounds the precision.\n");
  display_option_help_noopt("0 - computed from the size of the modular\n");
  display_option_help_noopt("    computation (default).\n");
  display_option_help('m', "", "MPR", "Maximal number of pairs used per matrix.\n");
  display_option_help_noopt("0 - unlimited (default).\n");
  display_option_help('q', "", "Q", "Uses signature-based algorithms.\n");
  display_option_help_noopt("0 - no (default).\n");
  display_option_help_noopt("1 - yes.\n");
  display_option_help(0, "random-seed", "SEED", "Random seed to initialize the pseudo-random\n");
  display_option_help_noopt("generator.\n");
  display_option_help_noopt("-1       - time(0) will be used (default).\n");
  display_option_help_noopt("SEED ≥ 0 - use at your own risk; this is\n");
  display_option_help_noopt("           intended for developers and debug\n");
  display_option_help_noopt("           purposes only.\n");
  display_option_help('s', "", "HTS", "Initial hash table size given as a power of two.\n");
  display_option_help_noopt("17 (default).\n");
  display_option_help(0, "separating-element", "SEP", "Over the rationals, when the last variable does\n");
  display_option_help_noopt("not separate the solutions: search for a\n");
  display_option_help_noopt("separating variable or linear form on the\n");
  display_option_help_noopt("Groebner basis already computed.\n");
  display_option_help_noopt("0 - no search, see option -c.\n");
  display_option_help_noopt("1 - search (default).\n");
  display_option_help(0, "split", "SPL", "Over the rationals: splitting of the solutions\n");
  display_option_help_noopt("along equalities of coordinates, each component\n");
  display_option_help_noopt("being lifted separately.\n");
  display_option_help_noopt("0 - when the largest component has at most a\n");
  display_option_help_noopt("    quarter of the solutions (default).\n");
  display_option_help_noopt("1 - never.\n");
  display_option_help_noopt("2 - whenever there are several components.\n");
  display_option_help('u', "", "UHT", "Number of steps after which the hash table is\n");
  display_option_help_noopt("newly generated.\n");
  display_option_help_noopt("0 - no update (default).\n");
  fprintf(OUTSTREAM, "\nOther options:\n\n");
  display_option_help('h', "help", "", "Prints this help.\n");
  display_option_help('V', "version", "", "Prints msolve's version.\n");
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
        int32_t *sep_search,
        int32_t *sep_dense,
        int32_t *split_mode,
        int32_t *lifting,
        int32_t *hensel_prec,
        double *hensel_cost_scale,
        double *hensel_memory,
        int32_t *precision,
        int32_t *generate_pbm_files,
	int64_t *seed,
        int32_t *boolean,
        int32_t *boolean_zeros,
        int32_t *lex,
        int32_t *elim_drl,
        int32_t *info_level,
        files_gb *files){
  int opt, errflag = 0, fflag = 1;
  char *filename = NULL;
  char *out_fname = NULL;
  opterr = 1;
  char short_options[] = "c:Cd:e:f:g:hil:L:m:M:n:N:o:p:P:q:r:s:St:u:v:V";

  /* For long options that have no equivalent short option, use a
     non-character as a pseudo short option, starting with CHAR_MAX + 1.
     see https://cgit.git.savannah.gnu.org/cgit/coreutils.git/tree/src/ls.c */
  enum {
    RANDOM_SEED_OPTION = CHAR_MAX + 1,
    SEPARATING_ELEMENT_OPTION,
    DENSE_ROUTE_OPTION,
    SPLIT_OPTION,
    LIFTING_PARAM_OPTION,
    HENSEL_PRECISION_OPTION,
    HENSEL_COST_SCALE_OPTION,
    HENSEL_MEMORY_OPTION,
    BOOLEAN_OPTION,
    BOOLEAN_ZEROS_OPTION,
    LEX_OPTION,
    ELIM_DRL_OPTION,
    /* Below is the template for the next long option with
     * no equivalent short option */
    /* NEXT_OPTION */
  };
  struct option long_options[] = {
    {"boolean", no_argument, NULL, BOOLEAN_OPTION},
    {"boolean-zeros", no_argument, NULL, BOOLEAN_ZEROS_OPTION},
    {"elimination", required_argument, NULL, 'e'},
    {"elim-drl", required_argument, NULL, ELIM_DRL_OPTION},
    {"file", required_argument, NULL, 'f'},
    {"groebner-basis", required_argument, NULL, 'g'},
    {"help", no_argument, NULL, 'h'},
    {"lex", no_argument, NULL, LEX_OPTION},
    {"linear-algebra", required_argument, NULL, 'l'},
    {"lifting-param", required_argument, NULL, LIFTING_PARAM_OPTION},
    {"hensel-precision", required_argument, NULL, HENSEL_PRECISION_OPTION},
    {"hensel-cost-scale", required_argument, NULL, HENSEL_COST_SCALE_OPTION},
    {"hensel-memory", required_argument, NULL, HENSEL_MEMORY_OPTION},
    {"lifting-mulmat", required_argument, NULL, 'L'},
    {"normal-form", required_argument, NULL, 'n'},
    {"output-file", required_argument, NULL, 'o'},
    {"precision", required_argument, NULL, 'p'},
    {"parametrization", required_argument, NULL, 'P'},
    {"random-seed", required_argument, NULL, RANDOM_SEED_OPTION},
    {"separating-element", required_argument, NULL, SEPARATING_ELEMENT_OPTION},
    {"dense-route", required_argument, NULL, DENSE_ROUTE_OPTION},
    {"split", required_argument, NULL, SPLIT_OPTION},
    {"reduce-gb", required_argument, NULL, 'r'},
    {"threads", required_argument, NULL, 't'},
    {"verbose", required_argument, NULL, 'v'},
    {"version", no_argument, NULL, 'V'},
    /* Below is the template for the next long option with
     * no equivalent short option,
     * adapt the 2nd field accordingly */
    /* {"next-option", required_argument, NULL, NEXT_OPION}, */
    {NULL,0,NULL,0}
  };

  while(1) {
    int option_index = 0;
    opt = getopt_long(argc, argv, short_options, long_options, &option_index);
    if (opt == -1) {
      /* processed all command-line options */
      break;
    }

    switch(opt) {
    case BOOLEAN_OPTION:
      *boolean = 1;
      break;
    case BOOLEAN_ZEROS_OPTION:
      *boolean_zeros = 1;
      break;
    case LEX_OPTION:
      *lex = 1;
      break;
    case ELIM_DRL_OPTION:
      *elim_drl = strtol(optarg, NULL, 10);
      if (*elim_drl < 0 || *elim_drl > 2) {
          *elim_drl = 2;
      }
      break;
    case 'N':
      *truncate_lifting = strtol(optarg, NULL, 10);
      break;
    case 'h':
      display_help(argv[0]);
      exit(0);
    case 'V':
      fprintf(OUTSTREAM, "%s\n", VERSION);
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
    case 'i':
      *is_gb = 1;
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
    case SEPARATING_ELEMENT_OPTION:
      *sep_search = strtol(optarg, NULL, 10) != 0;
      break;
    case DENSE_ROUTE_OPTION:
      *sep_dense = strtol(optarg, NULL, 10);
      if (*sep_dense < 0 || *sep_dense > 2) {
          *sep_dense = 0;
      }
      break;
    case SPLIT_OPTION:
      *split_mode = strtol(optarg, NULL, 10);
      if (*split_mode < 0 || *split_mode > 2) {
          *split_mode = 0;
      }
      break;
    case LIFTING_PARAM_OPTION:
      *lifting = strtol(optarg, NULL, 10);
      if (*lifting < 0 || *lifting > 2) {
          *lifting = 0;
      }
      break;
    case HENSEL_PRECISION_OPTION:
      *hensel_prec = strtol(optarg, NULL, 10);
      if (*hensel_prec < 0) {
          *hensel_prec = 0;
      }
      if (*hensel_prec > 1024) {
          *hensel_prec = 1024;
      }
      break;
    case HENSEL_COST_SCALE_OPTION:
      *hensel_cost_scale = strtod(optarg, NULL);
      if (*hensel_cost_scale <= 0) {
          *hensel_cost_scale = 1.0;
      }
      break;
    case HENSEL_MEMORY_OPTION:
      *hensel_memory = strtod(optarg, NULL);
      if (*hensel_memory < 0) {
          *hensel_memory = 0;
      }
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
    case 'o':
      out_fname = optarg;
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
    case RANDOM_SEED_OPTION:
      *seed = strtoll(optarg, NULL, 10);
      break;
    /* Below is the template for the next long option with
     * no equivalent short option */
    /* case NEXT_OPTION: */
    /*   break; */
    default:
      errflag++;
      break;
    }
  }
  if(fflag){
    fprintf(ERRSTREAM,"No given file\n");
    display_help(argv[0]);
    exit(1);
  }
  if(errflag){
    fprintf(ERRSTREAM, "Invalid usage\n");
    display_help(argv[0]);
    exit(1);
  }
  files->in_file = filename;
  files->out_file = out_fname;
}

/**
 * Checks whether the first line of file fn ends with "\r\n" (DOS line ending)
 *
 * \return 1 if it does, 0 otherwise (including if fn cannot be read)
 */
static int has_crlf_line_ending(const char *fn) {
  FILE *fh = fopen(fn, "rb");
  if (fh == NULL) {
    return 0;
  }
  int c, prev = EOF;
  while ((c = fgetc(fh)) != EOF && c != '\n') {
    prev = c;
  }
  fclose(fh);
  return c == '\n' && prev == '\r';
}

/**
 * Rewrites file fn in place, replacing each "\n" which is not already
 * preceded by "\r" with "\r\n"
 *
 * \return 0 on success, 1 on failure
 */
static int convert_to_crlf_line_endings(const char *fn) {
  FILE *fh = fopen(fn, "rb");
  if (fh == NULL) {
    return 1;
  }
  FILE *tmp = tmpfile();
  if (tmp == NULL) {
    fclose(fh);
    return 1;
  }
  int c;
  while ((c = fgetc(fh)) != EOF) {
    fputc(c, tmp);
  }
  fclose(fh);
  rewind(tmp);

  fh = fopen(fn, "wb");
  if (fh == NULL) {
    fclose(tmp);
    return 1;
  }
  int prev = EOF;
  while ((c = fgetc(tmp)) != EOF) {
    if (c == '\n' && prev != '\r') {
      fputc('\r', fh);
    }
    fputc(c, fh);
    prev = c;
  }
  fclose(tmp);
  return fclose(fh) != 0;
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
    int32_t sep_search            = 1;
    int32_t sep_dense             = 0;
    int32_t split_mode            = 0;
    int32_t lifting               = 0;
    int32_t hensel_prec           = 0;
    double hensel_cost_scale      = 1.0;
    double hensel_memory          = 0;
    int32_t precision             = 64;
    int64_t seed                  = -1;
    int32_t boolean               = 0;
    int32_t boolean_zeros         = 0;
    int32_t lex                   = 0;
    int32_t elim_drl              = 2;

    files_gb *files = malloc(sizeof(files_gb));
    if(files == NULL) exit(1);
    files->in_file = NULL;
    files->out_file = NULL;
    getoptions(argc, argv, &initial_hts, &nr_threads, &max_pairs,
               &elim_block_len, &la_option, &use_signatures, &update_ht,
               &reduce_gb, &print_gb, &truncate_lifting, &genericity_handling,
               &unstable_staircase, &saturate, &colon,
               &normal_form, &normal_form_matrix, &is_gb, &lift_matrix, &get_param,
               &sep_search, &sep_dense, &split_mode,
               &lifting, &hensel_prec, &hensel_cost_scale,
               &hensel_memory, &precision, &generate_pbm,
	       &seed, &boolean, &boolean_zeros, &lex, &elim_drl, &info_level,
               files);

    /* srand initialization */
    uint32_t true_seed;
    if (seed < 0) {
      true_seed = time(0);
    } else {
      true_seed = (uint32_t) seed;
    }
    srand(true_seed);
    if (info_level) {
      fprintf (VERBSTREAM,"Initial seed for pseudo-random number generator ");
      fprintf (VERBSTREAM,"is %u\n",true_seed);
    }

    FILE *fh = fopen(files->in_file, "r");
    if (fh == NULL) {
      fprintf(ERRSTREAM, "Input file not found.\n");
      exit(1);
    }
    fclose(fh);

    /* clear out_file if given */
    if(files->out_file != NULL){
      FILE *ofile = fopen(files->out_file, "wb");
      if(ofile == NULL){
        fprintf(ERRSTREAM, "Cannot open output file\n");
        exit(1);
      }
      fclose(ofile);
    }
    /**
       We get from files the requested data.
    **/
    //  int32_t mon_order   = 0;
    int32_t nr_vars     = 0;
    int64_t field_char  = 9001;
    int32_t nr_gens     = 0;
    data_gens_ff_t *gens = allocate_data_gens();

    get_data_from_file(files->in_file, &nr_vars, &field_char, &nr_gens, gens);
    if (boolean) {
      if (field_char != 2) {
        fprintf(ERRSTREAM, "The boolean ring (option --boolean) needs ");
        fprintf(ERRSTREAM, "the field characteristic 2.\n");
        exit(1);
      }
      gens->boolean = 1;
    }
    if (boolean_zeros) {
      if (field_char != 2) {
        fprintf(ERRSTREAM, "The computation via the common zeros (option ");
        fprintf(ERRSTREAM, "--boolean-zeros) needs the field characteristic 2.\n");
        exit(1);
      }
      gens->boolean       = 1;
      gens->boolean_zeros = 1;
      gens->lex           = lex;
    } else if (lex) {
      fprintf(ERRSTREAM, "The lexicographical order (option --lex) is only ");
      fprintf(ERRSTREAM, "available with option --boolean-zeros.\n");
      exit(1);
    }
#ifdef IODEBUG
    display_gens(VERBSTREAM, gens);
#endif

    gens->rand_linear           = 0;
    gens->random_linear_form = malloc(sizeof(int32_t)*(nr_vars));
    gens->elim = elim_block_len;
    gens->elim_drl = elim_drl;
    /* the parametrization of the whole solution set is not built from
     * its components when it is not printed */
    gens->split_global = get_param > 0;
    gens->sep_search = sep_search;
    gens->sep_dense = sep_dense;
    gens->split_mode = split_mode;
    gens->lifting = lifting;
    gens->hensel_prec = hensel_prec;
    gens->hensel_cost_scale = hensel_cost_scale;
    gens->hensel_memory = hensel_memory;

    if(0 < field_char && field_char < pow(2, 15) && la_option > 2){
        if(info_level){
            fprintf(ERRSTREAM, "Warning: characteristic is too low for choosing \nprobabilistic linear algebra\n");
            fprintf(ERRSTREAM, "\t linear algebra option set to 2\n");
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

    /* output files are written with "\n" line endings; if the input file
     * uses "\r\n", the output file is converted afterwards to match it */
    int crlf_output = files->in_file != NULL && files->out_file != NULL
                      && has_crlf_line_ending(files->in_file);

    /* main msolve functionality */
    int ret = core_msolve(la_option, use_signatures, nr_threads, info_level,
                          initial_hts, max_pairs, elim_block_len, update_ht,
                          generate_pbm, reduce_gb, print_gb, truncate_lifting,
                          get_param,
                          genericity_handling, unstable_staircase, saturate,
                          colon, normal_form,
                          normal_form_matrix, is_gb, lift_matrix, precision,
                          files, gens,
                          &param, mpz_paramp, &nb_real_roots, &real_roots,
                          &real_pts);

    if (crlf_output && convert_to_crlf_line_endings(files->out_file)) {
        fprintf(ERRSTREAM, "Cannot convert line endings of output file\n");
    }

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
        fprintf(VERBSTREAM, "\n\n-------------------------------------------------\
-----------------------------------\n");
        fprintf(VERBSTREAM, "msolve overall time  %13.2f sec (elapsed) / %5.2f sec (cpu)\n",
                rt1-rt0, st1-st0);
        fprintf(VERBSTREAM, "-------------------------------------------------\
-----------------------------------\n");
    }
    free_data_gens(gens);

    free(files);
    return ret;
}
