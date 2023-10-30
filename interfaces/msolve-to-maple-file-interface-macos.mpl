# /* This file is part of msolve.
#  *
#  * msolve is free software: you can redistribute it and/or modify
#  * it under the terms of the GNU General Public License as published by
#  * the Free Software Foundation, either version 2 of the License, or
#  * (at your option) any later version.
#  *
#  * msolve is distributed in the hope that it will be useful,
#  * but WITHOUT ANY WARRANTY; without even the implied warranty of
#  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  * GNU General Public License for more details.
#  *
#  * You should have received a copy of the GNU General Public License
#  * along with msolve.  If not, see <https://www.gnu.org/licenses/>
#  *
#  * Authors:
#  * Christian Eder
#  * Jorge Garcia Fontan
#  * Huu Phuoc Le
#  * Mohab Safey El Din */

FormatOutputMSolve:=proc(ll, _Z)
local dim, deg, degquot, params, nvars, elim, den, cfs, i, varstr, linearform:
  dim:=ll[1]:
  if dim > 0 then
  return -1, [], [];
  fi;
  nvars:=ll[2]:
  degquot:=ll[3]:
  varstr:=ll[4]:
  linearform:=ll[5]:
  deg:=ll[6][1]:
  if nops(ll[6]) > 0 then
    elim:=PolynomialTools[FromCoefficientList](ll[6][2],_Z);
  else #computation failed
    return -2, [], [];
  fi:
  if nops(ll[5]) > 0 then
    den:=PolynomialTools[FromCoefficientList](ll[7][2],_Z);
  fi:
  params:=[]:
  cfs:=[1]:
  if degquot > 0 then
  for i from 1 to nvars-1 do
      params:=[op(params),
              PolynomialTools[FromCoefficientList](ll[8][i][1][2],_Z)]:
    cfs:=[op(cfs), ll[8][i][2]]:
  od:
  fi:
  return varstr, linearform, elim, den, [diff(elim, _Z), op(params)], cfs;
end:



GetRootsFromMSolve:=proc(l)
local _Z, e, d, p, c, v, lf, sols, vals, i, sols_and_vals, realroots, OldDigits,
mysols;
local exist_gcd_p_e_nontrivial, length_p, gcd_p_e;
  v, lf, e, d, p, c := FormatOutputMSolve(l, _Z);
  if degree(e) = 0 then #positive dimension (-1) or computation failed (-2)
    return e, [];
  else
    return 0, [v,lf,e,d,p,c];
  fi:
end proc:


ToMSolve:=proc(F, char, vars, fname)
local i, fd, F2, str;
   fd:=fopen(fname, WRITE):
   for i from 1 to nops(vars)-1 do
     fprintf(fd, "%a, ", vars[i]):
   end;
   fprintf(fd, "%a ", vars[nops(vars)]):
   fprintf(fd,"\n");
   fprintf(fd,"%d\n", char);
   if char = 0 then
     F2:=map(f->sort(expand(numer(f)), order=tdeg(op(vars))), F):
     F2:=remove(member, F2, [0]):
     for i from 1 to nops(F2)-1 do
       fprintf(fd, "%a,\n", F2[i]):
     od:
     fprintf(fd, "%a\n", F2[nops(F2)]):
   else
     F2:=map(f->sort(expand(numer(f)), order=tdeg(op(vars))) mod char, F):
     F2:=remove(member, F2, [0]):
     for i from 1 to nops(F2)-1 do
        fprintf(fd, "%a,\n", F2[i]):
     od:
     fprintf(fd, "%a\n", F2[nops(F2)]):
   fi:
   fclose(fd):
#   str := cat("sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\n//g' ", fname):
   str := cat("sed -i '' -e ':a' -e 'N' -e '$!ba' -e 's/\\\\\\n//g' ", fname):
   lprint(str);
   system(str):
end proc:

GetOptions:=proc(opts)
local str, msolve_path, fname1, fname2, file_dir, verb, param, nthreads, output, gb;
  str:=subs(opts,"verb");
  if type(str, integer) then
     verb:=str;
  else
      verb:=0:
  end if;

  str:=subs(opts,"leadmons");
  if type(str, integer) and str > 0 then
     gb:=str;
  else
     gb:=2:
  end if;

  str:=subs(opts,"elim");
  if type(str, integer) and str > 0 then
     elim:=str;
  else
     elim:=0:
  end if;

  str:=subs(opts,"output");
  if type(str, integer) then
     output:=str;
  else
     output:=0:
  end if;

  str:=subs(opts,"nthreads");
  if type(str, integer) then
     nthreads:=str;
  else
     nthreads:=1:
  end if;

  if member("mspath", map(lhs, opts)) then
     str:=subs(opts, "mspath");
     if type(str, string) then
       msolve_path:=str;
     else
       printf("Error in format options");
     end if;
  else
     msolve_path := "../msolve";
  end if;

  if member("file_dir", map(lhs, opts)) then
     str:=subs(opts, "file_dir");
     if type(str, string) then
        file_dir:=str;
     else
        printf("Error in format options");
     end if;
  else
     file_dir := "/tmp/";
  end if;

  if member("file_in", map(lhs, opts)) then
     str:=subs(opts, "file_in");
     if type(str, string) then
        fname1:=cat(file_dir,str);
     else
        printf("Error in format options");
     end if;
  else
      fname1 := cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");;
  end if;

  if member("file_out", map(lhs, opts)) then
     str:=subs(opts, "file_out");
     if type(str, string) then
        fname2:=cat(file_dir,str);
     else
        printf("Error in format options");
     end if;
  else
     fname2 := cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");;
  end if;

  param:=0;
  msolve_path, fname1, fname2, file_dir, verb, param, nthreads, output, gb, elim;
end proc;


#Input.
#F: list of polynomials with rational coefficients
#fc: field characteristic
#vars: list of variables
#Optional argument:
#   {"mspath"=<path to msolve binary>,
#    "verb"=<positive integer for verbosity>,
#    "file_dir"=<repository where intermediate files will be written>,
#    "file_in"=<name of file for input systems read by msolve>,
#    "file_out"=<name of file where msolve writes computed results>,
#    "leadmons"= <1 to get leading monomials only>
#    "elim"=<number of variables to eliminate>}
#Output.
#[] -> an error occured during the computation
#else returns a Groebner basis of the ideal generated by F for the 
# grevlex ordering over monomials involving variables in vars
# with vars[1] > ... > vars[n]
# if "leadmons"=1 is part of the third optional argument, then only
# the leading monomials are returned
MSolveGroebner:=proc(F, fc, vars, opts:={})
local results, dim, fname1, fname2, verb, param, msolve_path, file_dir,
field_char, lsols, nl, i, gb, output, nthreads, str, elim;
   if type(F, list(polynom(rational))) = false then
     printf("First argument is not a list of polynomials with rational coefficients\n");
   end if;

   if type(fc, integer)=true then
      if fc < 0 then
         error "Field characteristic cannot be negative";
      end if;

      if fc > 0 and isprime(fc)=false then
         error "Field characteristic should be a prime number";
      end if;

      if fc > 2^31 then
         error "Field characteristic is too large to be supported";
      end if;

      if fc > 0 then
         field_char := fc;
      end if;

   else
      printf("Second argument should be a prime integer < 2^31\n");
   end if;
   if not(indets(F) subset indets(vars)) then
     printf("Given variables do not match the variables in the input polynomials\n");
   end if;

   msolve_path:="../msolve";
   file_dir:="/tmp/";
   fname1:=cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");
   fname2:=cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");
   while evalb(fname1=fname2) do
         fname2:=cat(RandomTools[Generate](string(8,alpha)), ".ms");
   od:
   verb:=0;
   param:=0:
   nthreads:=1;
   output:=0;
   gb:=2;
   if nops(opts) > 0 then
     msolve_path, fname1, fname2, file_dir, verb, param, nthreads, output, gb, elim := GetOptions(opts);
   fi;
   ToMSolve(F, field_char, vars, fname1);
   str := cat(msolve_path, " -v ", verb, " -g ", gb, " -e ", elim, " -P ", param, " -t ", nthreads, " -f ", fname1, " -o ", fname2):

   try
   system(str):
   read(fname2):
   catch:
   printf("There has been an issue in msolve computation\n");
   return [];
   end:

   results:=%:
   system(cat("rm ", fname2));
   system(cat("rm ", fname1));
   return results;
end proc:


#Input.
#F: list of polynomials with rational coefficients
#vars: list of variables
#Optional argument:
#   {"mspath"=<path to msolve binary>,
#    "verb"=<positive integer for verbosity>,
#    "file_dir"=<repository where intermediate files will be written>,
#    "file_in"=<name of file for input systems read by msolve>,
#    "file_out"=<name of file where msolve writes computed results>,
#    "output"= <1 to get mid points of isolating boxes>}
#Output.
#[] -> an error occured during the computation
#[1] -> input system has infinitely many complex solutions
#[-1, []] -> input polynomial system has no real solution
#[0, [sols]] -> input polynomial system has finitely many complex solutions
# sols is the list of real solutions given with isolating boxes with the following format
# [vars[1] = [a1, b1], ..., vars[n] = [an,bn]]
#if "output"=1 is part of the third (optional) argument output format is 
# [ vars[1] = (a1+b1)/2, ..., vars[n] = (an+bn)/2 ]
MSolveRealRoots:=proc(F, vars, opts:={})
local results, dim, fname1, fname2, verb, param, msolve_path, file_dir,
lsols, nl, i, j, gb, output, nthreads, str, sols, prec;
   if type(F, list(polynom(rational))) = false then
     printf("First argument is not a list of polynomials with rational coefficients\n");
   end if;
   if not(indets(F) subset indets(vars)) then
     printf("Given variables do not match the variables in the input polynomials\n");
   end if;

   msolve_path:="../msolve";
   file_dir:="/tmp/";
   fname1:=cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");
   fname2:=cat(file_dir, RandomTools[Generate](string(8,alpha)), ".ms");
   while evalb(fname1=fname2) do
         fname2:=cat(RandomTools[Generate](string(8,alpha)), ".ms");
   od:
   verb:=0;
   param:=0:
   nthreads:=1;
   output:=0;
   gb:=0;
   if nops(opts) > 0 then
     msolve_path, fname1, fname2, file_dir, verb, param, nthreads, output, gb := GetOptions(opts);
   fi;
   ToMSolve(F, 0, vars, fname1);
   if Digits <= 10 then
     prec:=128:
   else
     prec:= iquo(Digits/10)*128:
   fi:
   str := cat(msolve_path, " -v ", verb, " -P ", param, " -p ", prec, " -t ", nthreads, " -f ", fname1, " -o ", fname2):
   gb:=0; #Needed to avoid the user stops GB comp once a prime computation is done
   param:=0;
   try
   system(str):
   read(fname2):
   catch:
   printf("There has been an issue in msolve computation\n");
   return [];
   end:

   results:=%:
   system(cat("rm ", fname2));
   system(cat("rm ", fname1));
   if nops(results) = 0 then
      printf("There has been an issue in msolve computation\n");
      return [];
   end if;

   dim := results[1];
   if dim = -1 then
      if verb >= 1 then
         printf("Algebraic set is empty\n");
      end if;
      return [-1, []];
   end if;
   if dim > 0 then
      if verb >= 1 then
         printf("Algebraic set has positive dimension\n");
      end if;
      return [1];
   end if;
   if dim = 0 then
      lsols := results[2];
      nl := lsols[1]:
      sols:=[]:
      for i from 1 to nl do
          sols:=[op(sols), op(map(_p->[seq(vars[j] = _p[j], j=1..nops(vars))], lsols[i+1]))];
      od:
      if output=1 then
         sols := map(_p->map(_c->lhs(_c)=(rhs(_c)[1]+rhs(_c)[2])/2, _p), sols);
      end if;
      return [0, sols];
   end if;
   return results;
end proc:


# #Example: Katsura-6

# #Input data
# F:=[x1+2*x2+2*x3+2*x4+2*x5+2*x6-1,
# x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2-x1,
# 2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6-x2,
# x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6-x3,
# 2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6-x4,
# x3^2+2*x2*x4+2*x1*x5+2*x2*x6-x5];
# vars:=[x1,x2,x3,x4,x5,x6];

# #Usage

# #with rational parametrization
# param, sols:=MSolveRealRoots(F,vars,"../binary/msolve","/tmp/in.ms","/tmp/out.ms",1);
# #or with solutions only
# sols:=MSolveRealRoots(F,vars);
