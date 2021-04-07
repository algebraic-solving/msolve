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

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<gmp.h>
#include<math.h>

//#include "univ_io.h"

void mpz_poly_print_mathematica_file(FILE *file, mpz_t *pol, unsigned long int deg){
  fprintf(file, "pol =");
  for(long i=deg;i>=1;i--){
    gmp_fprintf(file, "(%Zd)", pol[i]);
    fprintf(file,"*x^%ld+",i);
  }
  gmp_fprintf(file,"(%Zd)", pol[0]);
  fprintf(file, ";\n");

  fprintf(file, "Timing[RootIntervals[pol]]");
}

mpz_t *myget_poly_from_file(FILE *file, unsigned long int *deg_ptr, int rev){
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

int main(int argc, char **argv){
  if(argc!=3 && argc!=4){
    fprintf(stderr, "Usage: %s file_in file_out\n%s file_in file_out b (if input file is in binary form)\n", argv[0], argv[0]);
    return 1;
  }
  FILE * file_in = fopen(argv[1],"r");
  FILE * file_out = fopen(argv[2],"w");

  unsigned long int deg;
  long int size = 0;
  mpz_t p;
  mpz_t *upoly = NULL;
  if(argc==3){
    upoly = myget_poly_from_file(file_in, &deg, 1);
  }
  else{
    //    binary_read_INT(file_in,p);
    //    upoly=binary_read(file_in,&size);
    deg = (unsigned long int)(size - 1);
    mpz_clear(p);
  }

  mpz_poly_print_mathematica_file(file_out, upoly, deg);

  for(unsigned long int i = 0; i <= deg; i++){
    mpz_clear(upoly[i]);
  }
  free(upoly);

  fclose(file_in);
  fclose(file_out);
}
