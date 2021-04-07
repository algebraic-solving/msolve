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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <gmp.h>
#include <math.h>
#include <time.h>
#include<omp.h>

#define USEFLINT 1
#define CHECKSHIFT 0
/*
#include "bisection.h"
#include "print_usolve.h"
//#include "univ_io.h"
#include "refine.h"
#include "debug.h"
*/

#include "usolve.c"

//reads poly from file ; when rev = 1 coeffs are given in increasing degree order
static mpz_t *get_poly_from_file(FILE *file, unsigned long int *deg_ptr, int rev){
  if(fscanf(file, "%ldu\n", deg_ptr)){
  }
  else{
    fprintf(stderr, "Error when reading file containing the input polynomial\n");
    exit(1);
  }
  mpz_t *upoly = (mpz_t *)(malloc(sizeof(mpz_t) * ( (*deg_ptr) + 1)));
  if(rev == 1){
    for(unsigned long int i = 0; i <= *deg_ptr; i++){
      mpz_init(upoly[i]);
      mpz_inp_str(upoly[i], file, 10);
    }
  }
  else{
    for(long int i = *deg_ptr; i >= 0; i--){
      mpz_init(upoly[i]);
      mpz_inp_str(upoly[i], file, 10);
    }
  }
  return upoly;
}

static interval *get_roots_from_file(FILE *file, unsigned long int *nbroots_pos_ptr,
                              unsigned long int *nbroots_neg_ptr){
  unsigned long int nbroots_ptr = 0;
  if(fscanf(file, "%ldu\n", &nbroots_ptr)){
  }
  else{
    fprintf(stderr, "Error when reading file containing roots\n");
    exit(1);
  }
  interval *roots = (interval *)(malloc((nbroots_ptr) * sizeof(interval)));

  for(unsigned long int i = 0; i < (nbroots_ptr); i++){
    mpz_init((roots+i)->numer);
    mpz_inp_str((roots+i)->numer, file, 10);
    if(mpz_sgn((roots+i)->numer) <0){
      *nbroots_neg_ptr = *nbroots_neg_ptr +1;
    }
    else{
      *nbroots_pos_ptr = *nbroots_pos_ptr +1;
    }
    if(fscanf(file, "%ld", &((roots+i)->k))){
    }
    else{
      fprintf(stderr, "Error when reading file containing roots\n");
      exit(1);
    }
    if(fscanf(file, "%du", &((roots+i)->isexact))){
    }
    else{
      fprintf(stderr, "Error when reading file containing roots\n");
      exit(1);
    }
    (roots+i)->sign_left = 0;
  }
  return roots;
}




//TODO complete the documentation
static void print_help(char *progname){
  fprintf(stderr, "\nUsolve version %s\n\n", GIT_COMMIT_HASH);
  fprintf(stderr, "C library, implemented by M. Safey El DIn (Sorbonne Univ.)\n\n");
  fprintf(stderr,"Usage: %s [...] -f file\n", progname);
  fprintf(stderr,"where file contains first the degree of the input polynomial followed by \nits coefficients given by increasing degree\n\n");
  fprintf(stderr, "Warning: the program assumes that the input file defines a ** square-free ** polynomial\n");
  fprintf(stderr, "For instance, when the input file is \n3\n1\n-7\n1\n3\n");
  fprintf(stderr, "./usolve -f file \nreturns on the standard output the number and isolating intervals of the real roots of \nthe polynomial 1-7*x+x^2+3*x^3\n");
  fprintf(stderr, "3\n[[-115578/2^16, -115577/2^16], [9655/2^16, 9656/2^16], [336308/2^18, 336309/2^18]]\n");
  fprintf(stderr, "\n");
  fprintf(stderr,"-h\t displays this help\n");
  fprintf(stderr,"-V\t version number\n");
  fprintf(stderr,"-P\t isolates positive roots\n");
  fprintf(stderr,"-N\t isolates negative roots\n");
  fprintf(stderr,"-d\t coefficients are read in input file by decreasing degree\n");

  fprintf(stderr,"\nMore advanced options:\n");
  fprintf(stderr, "-v num     \tcontrols the verbosity (0<=num<=5)\n");
  fprintf(stderr,"-m          \tdisplays the list of roots ended with a semi-colon so that it can be read with many computer algebra systems \n");
  fprintf(stderr,"-s          \tdisplays statistics about the computation\n");
  fprintf(stderr,"-p num      \tincreases the requested (bit) precision on the isolating intervals (num must be an integer)\n");
  fprintf(stderr,"-e num      \tif one only wants to decide the existence of a real root, stops the computation as soon as one real root is found\n");
  fprintf(stderr,"-t num      \tmulti-threading\n");
  fprintf(stderr,"-o [ofile]  \tstores the real roots in the file 'ofile' using the internal format used by usolve\n");
  fprintf(stderr,"-r [rfile]  \treads a file containing real roots of the input polynomial in the input file and will refine those roots up to the requested precision (isolation is then skipped)\n");
  fprintf(stderr,"-b          \tinput file has binary format\n");

  fprintf(stderr,"\nExamples of use (to be completed):\n\n");
  fprintf(stderr,"./usolve -f file -o /tmp/sols16.out > sols16 \n--> will isolate the real roots up to 16 bits of precision ; those are stored in human readable format in sols16 (list of intervals) and in a less readable format in /tmp/sols16.out\n\n");
  fprintf(stderr,"./usolve -p 128 -f file /tmp/sols64.out -r /tmp/sols16.out > sols64 \n--> will refine the real roots up to 128 bits of precision, without repeating the isolation\n\n");
}


static char *getoptions(int argc, char **argv, usolve_flags *flags,
                 char **output_file, char **refine_file, int *system_format){
  int optch, errflag = 0, fflag = 1; 
	extern int opterr;

	char options[] = "NPDVmdhcsebv:t:f:p:o:r:";
  char *filename = NULL;
	opterr = 1;

	while ((optch = getopt(argc, argv, options)) != -1)
    switch (optch) {
    case 'h':
      print_help(argv[0]);
      exit(1);
    case 'm':
      *system_format=1;
      break;
    case 'e':
      flags->hasrealroots=1;
      break;
    case 'd':
      flags->revert = 0;
			break;
		case 'D':
      flags->debug = 1;
			break;
		case 'c':
      flags->classical_algo = 1;
			break;
		case 's':
      flags->print_stats = 1;
			break;
		case 'b':
      flags->bfile = 1;
			break;
    case 'V':
      fprintf(stderr, "Usolve version %s\n", GIT_COMMIT_HASH);
      exit(1);
      break;
    case 'N':
      flags->search = -1;
      break;
    case 'P':
      flags->search = 1;
      break;
		case 'f':
      fflag = 0;
      filename = optarg;
      optind--;
			break;
    case 'o':
      *output_file = optarg;
      optind--;
      break;
		case 'r':
      *refine_file = optarg;
      optind--;
			break;
		case 't':
      flags->nthreads = atoi(optarg);
      if(flags->nthreads<=0){
        flags->nthreads = 1;
      }
      optind--;
			break;
		case 'p':
      flags->prec_isole = atoi(optarg);
      if(flags->prec_isole<=0){
        flags->prec_isole = 16;
      }
      optind--;
			break;
		case 'v':
      flags->verbose = -1;
      flags->verbose = strtol(optarg, NULL, 10);
      optind--;
			break;
    case '?':
      errflag++;
      break;
    case ':':
      errflag++;
      break;
    default:
      errflag++;
      break;
    }

  if(fflag){
    fprintf(stderr,"No given file\n");
    print_help(argv[0]);
    exit(1);
  }
  if(errflag){
    fprintf(stderr,"Invalid usage\n");
    print_help(argv[0]);
    exit(1);
  }
	/* for (; optind < argc; ++optind) */
	/* 	printf ("argv[%d] : %s\n", optind, argv[optind]); */
  return filename;
}


int main(int argc, char **argv){

  usolve_flags *flags = (usolve_flags*)(malloc(sizeof(usolve_flags)));
  initialize_flags(flags);

  char *output_filename = NULL;
  char *refine_filename = NULL;
  int system_format = 0;
  char *input_filename =getoptions(argc, argv, flags, &output_filename, &refine_filename, &system_format);

  FILE *file = fopen(input_filename,"r");
  if(file==NULL){
    fprintf(stderr,"Problem: can not open file %s\n", input_filename);
    exit(1);
  }
  mpz_t *upoly = NULL;
  unsigned long int deg;

  if(flags->bfile==1){
    fprintf(stderr, "Not available yet\n");
  }
  else{
    upoly = get_poly_from_file(file, &deg, flags->revert);
  }
  while(mpz_sgn(upoly[deg])==0 && deg > 0){
    mpz_clear(upoly[deg]);
    deg--;
  }


  flags->cur_deg = deg;

  fprintf(stderr, "Univariate solving (Usolve version %s) \n", GIT_COMMIT_HASH);
  if(flags->verbose>=1 || flags->print_stats == 1){
    fprintf(stderr, "Degree = %ld \t Max bit size = %lu Min bit size = %lu \n",
            flags->cur_deg,
            USOLVEmpz_poly_max_bsize_coeffs(upoly, deg),
            USOLVEmpz_poly_min_bsize_coeffs(upoly, deg));
    fprintf(stderr, "nthreads = %d\n", flags->nthreads);
  }
  double e_time = omp_get_wtime ( );
  if(flags->debug==1){
    fprintf(stderr, "\n\tDEBUG MODE\n");
  }
  if(flags->classical_algo==1){
    fprintf(stderr, "WARNING: Using old algorithm\n");
  }

  #ifdef USEFLINT
  fprintf(stderr, "FLINT multiplication is used\n");
  #else
  //  fprintf(stderr, "WARNING: does not use FLINT multiplication\n");
  #endif


  unsigned long int nb_pos_roots = 0;
  unsigned long int nb_neg_roots = 0;
  unsigned long int nbroots = 0;
  interval *roots;

  if(refine_filename==NULL){
#ifdef NEWMODE
    roots = Bisection_Uspensky_from_root(upoly, deg, &nb_pos_roots, &nb_neg_roots, flags);
#else
    roots = Bisection_Uspensky(upoly, deg, &nb_pos_roots, &nb_neg_roots, flags);
#endif
    nbroots = nb_pos_roots + nb_neg_roots;
  }
  else{
    FILE *refinefile = fopen(refine_filename,"r");
    roots = get_roots_from_file(refinefile, &nb_pos_roots, &nb_neg_roots);
    nbroots = nb_neg_roots + nb_pos_roots;
    fclose(refinefile);
  }

  e_time = omp_get_wtime ( ) - e_time;

  if(flags->debug==1){
    fprintf(stderr, "\nRacines reelles apres isolation\n");
    USOLVEdisplay_roots(stderr, roots, nbroots);
    fprintf(stderr, "\n");
  }

  if(flags->debug==1){
    fprintf(stderr,"Checking intervals just after isolation\n");
    check_all_sign_changes(upoly, deg, roots, nbroots);
    fprintf(stderr, "Done\n");
  }

  if(flags->verbose>=1){
    fprintf(stderr, "\n");
  }

  if((flags->verbose>=1) || (flags->print_stats>=1)){
    fprintf(stderr,"Elapsed time (isolation): %f\n", e_time);
  }


  //Refinement starts
  if(flags->verbose==1){
    fprintf(stderr, "Refinement starts\n");
  }
  double refine_time = omp_get_wtime();
  double step = (e_time+1) / (deg) * 1000 * LOG2(MAX(16, flags->prec_isole)) * 2;

  if(nbroots>0 && flags->prec_isole>=0){
    if(flags->classical_algo>0){
      refine_all_roots_naive(upoly,deg, roots, nbroots,
                             flags->prec_isole, flags->classical_algo, flags->debug);
    }
    else{
      refine_QIR_roots(upoly, &deg, roots, nb_neg_roots, nb_pos_roots,
                       flags->prec_isole, flags->verbose, step, flags->nthreads);
    }
  }
  refine_time = omp_get_wtime() - refine_time;

  if(flags->print_stats>=1){
    display_stats(flags);
  }

  if(system_format==0){
    fprintf(stdout,"%lu\n", nbroots);
    USOLVEdisplay_roots(stdout, roots, nbroots);
  }
  else{
    display_roots_system(stdout, roots, nbroots);
  }
  if(flags->debug==1){
    fprintf(stderr,"Checking intervals just after refinement\n");
    check_all_sign_changes(upoly, deg, roots, nbroots);
    fprintf(stderr, "Done\n");
  }

  if((flags->verbose>=1) || (flags->print_stats>=1)){
    fprintf(stderr,"Elapsed time (isolation): %f\n", e_time);
    fprintf(stderr,"Elapsed time (refinement): %f\n", refine_time);
  }

  if(output_filename!=NULL){
    FILE *outputfile = fopen(output_filename,"w");
    USOLVEdisplay_roots_in_file(outputfile, roots, nb_neg_roots + nb_pos_roots);
    fclose(outputfile);
  }
  for(unsigned long int i=0; i<=deg;i++){
    mpz_clear(upoly[i]);
  }
  free(upoly);

  fclose(file);
  free(flags);
  return 0;
}



