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

uint32_t primes_table[10] = {2, 3, 4, 5, 7, 11, 13, 17, 19, 23};
int is_prime(uint32_t n){
  for(int i = 0; i < 10; i++){
    if(n % primes_table[i] == 0){
      return 0;
    }
  }
  for(uint32_t i = 5; i * i <= n; i += 6){
    if((n % i)==0) return 0;
    if( (n % (i+2)) == 0) return 0;
  }
  return 1;
}

// assumes n >= 2
uint32_t next_prime(uint32_t n){
  uint32_t cand = n + 1;
  while(!is_prime(cand)) cand++;
  return cand;
}
