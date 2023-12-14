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

static void display_left(FILE * stream, interval *root){

  if (root->k <= 0){
    mpz_out_str(stream, 10, root->numer);
  }
  else{
    mpz_out_str(stream, 10, root->numer);
    fprintf(stream, "/2^%ld", root->k);
  }

}


static void display_right(FILE * stream, interval *root, mpz_t tmp){

  if (root->k <= 0){
    mpz_set_ui(tmp, 1);
    mpz_mul_2exp(tmp, tmp, -root->k);
    mpz_add(tmp, root->numer, tmp);
    mpz_out_str(stream, 10, tmp);
  }
  else{
    mpz_add_ui(tmp, root->numer, 1);
    mpz_out_str(stream, 10, tmp);
    fprintf(stream, "/2^%ld", root->k);
  }

  fprintf(stream, "]");
}

static void display_left_in_file(FILE * stream, interval *root){

  mpz_out_str(stream, 10, root->numer);
  fprintf(stream, " ");
  fprintf(stream, "%ld ", root->k);
  fprintf(stream, "%d\n", root->isexact);
}


void display_root(FILE *stream, interval* root){
  mpz_t tmp;

  mpz_init(tmp);

  fprintf(stream, "[");
  display_left(stream, root);

  fprintf(stream, ", ");

  if (root->isexact == 1) {
    display_left(stream, root);
    fprintf(stream, "]");
    return;
  }

  display_right(stream, root, tmp);

  mpz_clear(tmp);
}



void USOLVEdisplay_roots(FILE *stream, interval *roots,
                                       unsigned long int nbroots){
  fprintf(stream,"[");
  for(unsigned long int i=0;i<nbroots;i++){
    display_root(stream,(roots+i));
    if(i<nbroots-1){
      fprintf(stream, ", ");
    }
  }
  fprintf(stream, "]\n");
}

void display_roots_system(FILE *stream, interval *roots, unsigned long int nbroots){
  fprintf(stream,"[");
  for(unsigned long int i=0;i<nbroots;i++){
    display_root(stream,(roots+i));
    if(i<nbroots-1){
      fprintf(stream, ", ");
    }
  }
  fprintf(stream, "];\n");
}


static inline void USOLVEdisplay_roots_in_file(FILE *stream, interval* roots,
                           unsigned long int nb_roots){
  fprintf(stream,"%lu\n", nb_roots);
  for(unsigned long int i = 0; i < nb_roots; i++){
    display_left_in_file(stream, (roots+i));
  }
  return;
}







