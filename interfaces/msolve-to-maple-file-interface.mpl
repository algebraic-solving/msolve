macro(_MSPATH="msolve"):

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
#  * Mohab Safey El Din 
#  * Bruno Salvy */

##Repositories to adapt 
homedir:=kernelopts(homedir);
savelibname:=cat(homedir, "/libs/"):
mladirname:=cat(homedir,"/libs/MSolve.mla");

# Installation 
# Run maple once on this file (after adjusting the above repositories)
# Recommended practice: adjust libname in your .mapleinit file to
# indicate maple to also look at savelibname

## Basic usage after installation: 
#with(MSolve);
#
#F:=[x+2*y+2*z-1,x^2+2*y^2+2*z^2-x,2*x*y+2*y*z-y]:
#G:=[x, x+z-2/3, y, x-1/3, y+x-1, y+z-1/3]:
#
## real root isolation of solutions of F and associated values of the
## polynomials in G
#sols:=MSolve:-MSolveRealRoots(F, [x,y,z], G);
#gb:=MSolve:-MSolveGroebner(F, 0, [x,y,z]):
#
## See below for more documentation on these functions (input/output +
## options)

MSolve:=module()
option package;

export MSolveGroebner, MSolveGroebnerLM,
MSolveRealRoots,MSolveParam;

local GetSystem, ToMSolve, GetOptions, CheckCharacteristic, ReadPolynomial,
ExtractParametrization, RemoveFiles, 
SplitParamCoord, 
CallMSolve, 
SelectVanishingSols, 
NonNegativeIntervalEvaluate, 
NonZeroIntervalEvaluate,
BuildSolution, Eval_linform, MakeBinaryInterval, 
Parametrization,
RefineSolutions,
SplitAndRefinePerCoordinates,
SplitAndRefinePerConstraints;



GetSystem:=proc()
local sys;
  sys:=kernelopts(system);
  if StringTools[Has](sys, "APPLE") then
    return "macOS";
  elif StringTools[Has](sys, "LINUX") then
    return "Linux";
  else
    return "Windows";
  fi;
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
   if evalb(GetSystem() = "macOS") then
     str := cat("sed -i '' -e ':a' -e 'N' -e '$!ba' -e 's/\\\\\\n//g' ", fname):
   elif evalb(GetSystem() = "Linux") then
     str := cat("sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\\\\n//g' ", fname):
   fi:
   ssystem(str);
end proc:

GetOptions:=proc(opts)
local seed, str, truncate, msolve_path, fname1, fname2, file_dir,
verb, param, nthreads, output, gb, elim, linalg;
  seed := randomize():
  str:=subs(opts,"verb");
  if type(str, integer) then
     verb:=min(str, 2);
     if verb < 2 then verb:= 0: end if;
  else
      verb:=0:
  end if;

  str:=subs(opts,"gb");
  if type(str, integer) and str > 0 then
     gb:=str;
  else
     gb:=0:
  end if;

  str:=subs(opts,"elim");
  if type(str, integer) and str > 0 then
     elim:=str;
  else
     elim:=0:
  end if;

  str:=subs(opts, "trunc");
  if type(str,integer) and str>0 then
      truncate := str:
  else truncate := -1:
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

  str:=subs(opts, "linalg"):
  if type(str, integer) then 
      linalg:=str:
  else 
      linalg:=2:
  end if;

  if member("mspath", map(lhs, opts)) then
     str:=subs(opts, "mspath");
     if type(str, string) then
       msolve_path:=str;
     else
       printf("Error in format options");
     end if;
  else
     msolve_path := _MSPATH;
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
      fname1 := cat(file_dir, RandomTools['Generate'](string(8,'alpha')), ".ms");;
  end if;

  if member("file_out", map(lhs, opts)) then
     str:=subs(opts, "file_out");
     if type(str, string) then
        fname2:=cat(file_dir,str);
     else
        printf("Error in format options");
     end if;
  else
     fname2 := cat(file_dir, RandomTools['Generate'](string(8,'alpha')), ".ms");;
  end if;

  str:=subs(opts,"param");
  if type(str,integer) then 
      param:=str 
  else 
      param:=0
  fi;

  str := cat(msolve_path, " -v ", verb, 
      `if`(gb=0 or param<>0,"",cat(" -g ", gb)), 
      `if`(elim=0,"",cat(" -e ", elim)),
      `if`(param=0,"",cat(" -P ", param)),
      `if`(truncate=-1,"",cat(" -N ", truncate)),
      `if`(nthreads=1,"",cat(" -t ", nthreads)),
      `if`(linalg=2,"",cat(" -l ", linalg)),
      " -f ", fname1, " -o ", fname2);

  str, fname1, fname2, verb, output
end proc:

RemoveFiles:=proc(fname1, fname2);
   ssystem(cat("rm ", fname2));
   ssystem(cat("rm ", fname1));
end proc:

CallMSolve:=proc(F, fc, vars, opts)
local nvars, newvars, prec, results, i, str, fname1, fname2, verb,
output, fd;

   str, fname1, fname2, verb, output := GetOptions(opts);

   nvars:=nops(vars);
   newvars:=[seq(_xx[i],i=1..nvars)];

   ToMSolve(subs([seq(vars[i]=_xx[i],i=1..nvars)],F), fc, newvars, fname1, _xx);
   if Digits <= 10 then prec:=64:
   else prec:= iquo(Digits,10)*64:
   fi:
   if fc = 0 then 
   str := cat(str, " -p ", prec);
   end if;

   try
      if verb=0 then ssystem(str)
      else system(str)
      fi;
      read(fname2);
      results:=subs([seq(_xx[i]=vars[i],i=1..nvars)],%)
   catch:
      fd:=fopen("/tmp/bug-call-msolve.mpl", WRITE):
      fprintf(fd, "F:=%a:\nfc:=%d\nvars:=%a:opts:=%a:\n", F, fc, vars,
      opts):
      fclose(fd);
      lprint(str);
      error "There has been an issue in msolve computation (see
      /tmp/bug-call-msolve.mpl)";
   end;

   RemoveFiles(fname1, fname2);

   return results, output;
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
MSolveGroebner:=proc(F::depends(list(polynom(rational,vars))), fc::nonnegint, vars::list(name), opts:={})
local fname1, fname2, verb, field_char, output, str, results, nvars, newvars, i;
   field_char := CheckCharacteristic(fc);

   str, fname1, fname2, verb, output := GetOptions(opts union {"gb"=2});
   nvars:=nops(vars);
   newvars:=[seq(_xx[i],i=1..nvars)];

   ToMSolve(subs([seq(vars[i]=_xx[i],i=1..nvars)],F), field_char, newvars, fname1);
   try
      if verb=0 then ssystem(str)
      else system(str)
      fi;
      read(fname2);
      results:=subs([seq(_xx[i]=vars[i],i=1..nvars)],%);
      RemoveFiles(fname1, fname2);
      return results;
   catch:
      error "There has been an issue in msolve computation"
   end
end proc:

MSolveGroebnerLM:=proc(F::depends(list(polynom(rational,vars))), fc::nonnegint, vars::list(name), opts:={})
local fname1, fname2, verb, field_char, output, str, results, nvars, newvars, i;
   field_char := CheckCharacteristic(fc);

   str, fname1, fname2, verb, output := GetOptions(opts union {"gb"=1});
   nvars:=nops(vars);
   newvars:=[seq(_xx[i],i=1..nvars)];

   ToMSolve(subs([seq(vars[i]=_xx[i],i=1..nvars)],F), field_char, newvars, fname1);
   try
      if verb=0 then ssystem(str)
      else system(str)
      fi;
      read(fname2);
      results:=subs([seq(_xx[i]=vars[i],i=1..nvars)],%);
      RemoveFiles(fname1, fname2);
      return results;
   catch:
      error "There has been an issue in msolve computation"
   end
end proc:

CheckCharacteristic:=proc(fc)
   if fc <> 0 and not isprime(fc) then
      error "Field characteristic should be a prime number"
   end if;
   if fc > 2^31 then
      error "Field characteristic > 2^31 not supported"
   end if;
   fc
end:

#pt is a list of intervals of rational numbers given in the format
#[x1 = [a1,b1], x2 = [a2, b2], ..., xn = [an, bn]]
#Assumes that all intervals defining pt are non-negative
#returns an interval containing the values that pol can reach at pt
NonNegativeIntervalEvaluate:=proc(pt, pol)
local monomials, cfs, poscfs, negcfs, posmonomials, negmonomials, lowpt, uppt, 
lowvalneg, i, lowval, upval, upvalneg, lowvalpos, upvalpos, vars:
  vars:= map(lhs, pt):
  cfs := coeffs(expand(pol), vars, 'monomials');
  cfs := [cfs]:
  monomials:=[monomials];
  poscfs:=map(c-> if c > 0 then c fi, cfs);
  negcfs:=remove(member, cfs, poscfs):
  posmonomials := [seq(if cfs[i]>0 then monomials[i] else fi, i=1..nops(cfs))]:
  negmonomials := remove(member, monomials, posmonomials):
  lowpt := map(p->lhs(p)=rhs(p)[1], pt); 
  uppt := map(p->lhs(p)=rhs(p)[2], pt); 

  lowval := subs(lowpt, posmonomials):
  upval := subs(uppt, posmonomials):
  lowvalpos := add(poscfs[i]*lowval[i],i=1..nops(poscfs));
  upvalpos := add(poscfs[i]*upval[i],i=1..nops(poscfs));

  lowval := subs(uppt, negmonomials):
  upval := subs(lowpt, negmonomials):
  lowvalneg := add(negcfs[i]*lowval[i],i=1..nops(negcfs));
  upvalneg := add(negcfs[i]*upval[i],i=1..nops(negcfs));
  return [lowvalpos + lowvalneg, upvalpos + upvalneg]:
end:

#pt is a list of intervals of rational numbers given in the format
#[x1 = [a1,b1], x2 = [a2, b2], ..., xn = [an, bn]]
#Assumes that none of the intervals [ai, bi] in pt contains 0 in its interior
NonZeroIntervalEvaluate:=proc(pt, pol)
local newpol, newpt, i, vars;
  if degree(pol)<= 0 then return [pol,pol] fi;
  vars := map(lhs, pt):
  newpol := pol:
  newpt := []:
  for i from 1 to nops(pt) do 
    if rhs(pt[i])[2] < 0 then 
      newpol:=subs(vars[i] = -vars[i], newpol):
      newpt:=[op(newpt), lhs(pt[i]) = [abs(rhs(pt[i])[2]),abs(rhs(pt[i])[1])]]:
    else 
      newpt:=[op(newpt),pt[i]]:
    end if;
  od:
  return NonNegativeIntervalEvaluate(newpt, newpol):
end:

ExtractParametrization:=proc(res, v)
local param, den, cfs, newres;

   param:=ReadPolynomial(res[1],v);
   den:=ReadPolynomial(res[2],v);
   cfs:=map2(op,-1,res[3]);
   newres:=map2(op,1,res[3]);
   newres:=map(ReadPolynomial,newres,v);

   return param, den, cfs, newres;
end proc;

SplitParamCoord:=proc(param, pol)
local param1, param2, lparam, newelim, _p;
  param1 := [pol, param[2], param[3], [op(param[4][1..c-1]), 0,
  op(param[4][c+1..-1])]]:
  param1 := [pol, param[2], param[3], map(_p->if divide(_p, pol) then
  0 else _p fi, param[4])]:
  lparam:=[param1];
  newelim := normal(param[1]/pol);
  if degree(newelim)>0 then 
  param2 := [newelim, param[2], param[3], param[4]]:
  lparam:=[op(lparam), param2];
  end if;
  return lparam;
end proc;

#sols is encoded such that all coordinates are given by intervals
#sols is the list of solutions from the original polynomial system
Eval_linform:=proc(lin_form, sols)
local v, lv, pt, vars, nlin, lc, _e;
  if nops(sols)=0 then return []; fi;
  lv := indets(lin_form);
  if nops(lv) = 1 then 
      v:= lv[1];
      return map(sol->subs(sol, v), sols);
  else 
      vars:=indets(sols);
      v:=remove(member, lv, vars):
      if nops(v) <> 1 then lprint(lin_form); error "lin_form incorrect" end
          if;
      v:=v[1]:
      lc:=lcoeff(lin_form, v):
      nlin:=solve(lin_form, v);
      return map(_e->[_e[1], _e[2]], 
                 map(pt->NonZeroIntervalEvaluate(pt, nlin), sols));
  end if;
end proc;

#vars is the list of variables from the original polynomial system
#param_in may depend on an extra variable (if lin_form is not a single variable)
SplitAndRefinePerCoordinates := proc(param_in, sols, lb, vars)
local i, idx, j, lparam, newlparam, g, splitparam, v, lsols, osols, k, vv,
param, interv, ldeg, vals, psols, intvl, newlb, elim, den, cfs, nums, lin_form,
psols2;

  lin_form := param_in[1]:
  elim := param_in[2]:
  den := param_in[3]:
  cfs := param_in[4]:
  nums := param_in[5]:
  idx := {}:
  for i from 1 to nops(lb) do 
      for j from 1 to nops(lb[i]) do
          if lb[i][j]=true then 
              idx :={op(idx), j}; 
          end if;
      end do;
  end do;

#idx is the set of indices of variables that have an undetermined sign
  lparam:=[[elim, den, cfs, nums]]:
  for i from 1 to nops(idx) do 
      newlparam :=[]:
      for j from 1 to nops(lparam) do 
        if idx[i]=nops(vars) and indets(lin_form)={vars[-1]}
          then 
          g:=gcd(lparam[j][1], lin_form):
        else 
          g:=gcd(lparam[j][1], lparam[j][-1][idx[i]]);
        fi:
        if degree(g) > 0 then 
          splitparam := SplitParamCoord(lparam[j], g);
          newlparam := [op(newlparam), op(splitparam)]:
        else
          newlparam := [op(newlparam), lparam[j]]:
        end if;
      end do;
      lparam := newlparam;
  end do;
  v := indets(elim)[1];
  lsols := []:
  osols := sols:
  for i from 1 to nops (lparam) do 
      param := lparam[i]:
      ldeg := [seq(if degree(param[-1][i])<0 then i else fi, i=1..nops(param[-1]))];
      interv := Eval_linform(lin_form, osols):
      vals := map(s->[subs(v=s[1], param[1]), subs(v=s[2], param[1])], interv);
      psols := []:
      for j from 1 to nops(vals) do 
        vv := vals[j];
        if ((vv[1]=0 and vv[2]=0) or (sign(vv[1])<>sign(vv[2]))) then 
            psols:=[op(psols), osols[j]]:
        end if;
      end do;
      osols:=remove(member, osols, psols);

      for j from 1 to nops(ldeg) do 
       for k from 1 to nops(psols) do 
           psols[k][ldeg[j]]:= lhs(psols[k][ldeg[j]])=[0,0]; 
       end do;
      end do;

      intvl := map(pt->map(rhs, pt), psols):
      newlb := map(_p->map(l->evalb((l[1] <= 0 and 0 <= l[2]) and (l[1]<>0 or l[2]<>0)), 
                        _p), intvl);
      psols2:=[]:
      for j from 1 to nops(psols) do 
        if not(member(true, newlb[j])) then 
            psols2:=[op(psols2), psols[j]]:
        end if;
      end do;
      psols:=remove(member, psols, psols2):
      if nops(map(l->select(member, true, l),newlb)) > 0 then 
#psols are the solutions to param and one knows that the coordinates we could
#not determine the sign are non-zero (splittings have been done)
          psols:=RefineSolutions(param, lin_form, psols, vars);
      end if;
      lsols:=[op(lsols), psols, psols2]:
  end do;

  return map(op, lsols);

end proc;

#return value boo is true if some coordinates have no determined sign yet
BuildSolution:=proc(param, interv, vars)
local i, j, coords, boo, var, den, nums, cfs, lb;
  var := indets(param[1])[1];
  den:=NonZeroIntervalEvaluate([var = interv], param[2]);
  cfs := param[3]:
  nums := [seq(NonZeroIntervalEvaluate([var=interv], param[-1][i]),
            i=1..nops(param[-1]))];
  coords := [seq([nums[i][1]/(den[2]*cfs[i]), nums[i][2]/(den[1]*cfs[i])], 
          i=1..nops(nums)), interv];
  lb := map(_c->(_c[1]<=0 and _c[2]>=0) and (_c[1]<>0 or _c[2]<>0), 
        coords);
  boo := member(true, lb, 'pos');
  if member(var, vars) then 
      return [op(map(i->if coords[i][2]<0 then
        vars[i]=[-coords[i][1],-coords[i][2]] else
        vars[i]=[-coords[i][2],-coords[i][1]] fi,
        [seq(j,j=1..nops(vars)-1)])),var=interv], 
        boo;
  else 
    return map(i->if coords[i][2]<0 then
            vars[i]=[-coords[i][1],-coords[i][2]] else
            vars[i]=[-coords[i][2],-coords[i][1]] fi, [seq(j,j=1..nops(vars))]),
            boo;
  end if;
end proc;

MakeBinaryInterval:=proc(interv)
local p1, p2;
  p1 := ilog2(denom(interv[1]));
  p2 := ilog2(denom(interv[2]));
  return [floor(interv[1]*2^p1)/2^p1, ceil(interv[2]*2^p2)/2^p2];
end proc;

RefineSolutions:=proc(param, lin_form, psols, vars)
local elim, v, i, newsols, sol, interv, prec, newsol, boo;
  elim := param[1]:
  v := indets(elim)[1]:
  newsols:=[]:
  for i from 1 to nops(psols) do 
    sol := psols[i]:
    lprint("sol",evalf(sol));
#Assumes all coefficients of lin_form are >=0
    interv := MakeBinaryInterval(Eval_linform(lin_form, [sol])[1]):
    lprint("interv", evalf(interv));
    prec := 2*(abs(ilog2(abs(interv[2]-interv[1])))+256);
    lprint("prec = ", prec):
    printf("{%d}",degree(elim));
    if interv[1]<>interv[2] then 
    interv := RootFinding[RefineRoot](interv, elim);#, digits = iquo(prec,4));
    end if;
    printf("[r]");
    newsol, boo := BuildSolution(param, interv, vars);
    gc();
    while boo = true and interv[2]<> interv[1] do 
      prec := 2*(abs(ilog2(abs(interv[2]-interv[1])))+256);
      printf("[r]");
      interv := RootFinding[RefineRoot](interv, elim);
      newsol, boo := BuildSolution(param, interv, vars);
    end do;
    newsols:=[op(newsols), newsol]:
  end do;
  return newsols;
end proc:

#Factoriser avec MSolveParam
Parametrization:=proc(F, vars, opts)
local verb, results, output, dim, v, i, res, nvars, newvars, lin_form, str, elim, den, cfs, nums;
   str := subs(opts, "param");
   verb := subs(opts, "verb");
   if not(type(verb, integer)) then 
       verb := 0;
   end if;
   if not(type(str, integer)) then
       error "Bad call to Parametrization (opts does not specify the param option)";
   end if;
   results, output := CallMSolve(F, 0, vars, opts);
   dim := results[1];
   if dim = -1 then
      if verb >= 2 then
         printf("Algebraic set is empty\n")
      end if;
      if str=2 then return vars, [1]; 
      else 
      return vars, [1], [-1, []];
      end if;
   end if;
   if dim > 0 then
      if verb >= 2 then
         printf("Algebraic set has positive dimension\n")
      end if;
      return vars, [0], [1, []];
   end if;

   res := results[2]:

   nvars:=res[2];
   newvars:=res[4];
   v:=newvars[-1];

   if not member(v,vars) then v:=_T; newvars[-1]:= _T; fi; # use a "clean" global name
   lin_form := add(res[5][i]*newvars[i], i=1..nops(newvars));
   res:=res[6];
   if res[1]<>1 then error "several output parametrizations" else res:=res[2] fi;
   elim, den, cfs, nums := ExtractParametrization(res, v);
   if str = 2 then 
   return newvars, [lin_form, elim, den, cfs, nums];
   else 
   return newvars, [lin_form, elim, den, cfs, nums], results[3];
   end if;
end proc;


SelectVanishingSols:=proc(pol, lin_form, sols)
local i, idx, interv, v1, v2, sol, var;
  idx:=[]:
  if degree(pol) = 0 then return []; end if;
  if pol = 0 then error "Input polynomial should not be 0"; end if;
  if nops(indets(pol)) > 1 then error "Input polynomial should be univariate";
  end if;
  var:=indets(pol)[1];
  interv := Eval_linform(lin_form, sols):
  for i from 1 to nops(sols) do 
      sol:=interv[i]:
      v1:=subs(var = sol[1], pol):
      v2:=subs(var = sol[2], pol):
      if sign(v1)*sign(v2)=-1 or (v1=0 or v2=0) then 
          idx:=[op(idx), i]:
      end if;
  end do;
  return idx;
end proc;

#elim is a univariate polynomial (the elimination polynomial of the
#parametrization)
#den is the denominator of the parametrization
#cfs is a list of integers used to multiply den when parametrizing coordinates
#nums are the numerators parametrizing the coordinates
#lin_form is the linear form used to compute the parametrization
#F is the list of polynomials for which the parametrization has been computed
#sols is the list of real solutions
#
#vals is the list of intervals of the values taken by some extra constraints per
#points encoded by sols
#cstr is the list of constraints for which the sign could not be determined
#idx is the index of these constraints, corresponding to elements in vals
SplitAndRefinePerConstraints:=proc(param_in,
                              sols, vals, F, vars, cstr, idx, opts)
local results, output, i, lparam, param, newsols, newvals, tsols, tvals, g,
idsols, j, newtsols, newtvals, lin_form, elim, den, cfs, nums, newvars, r,
pvars, n, gb, boo, interv, prec, newsol, val, vvar, new_lin_form;

  lin_form := param_in[1];
  elim := param_in[2]:
  den := param_in[3]:
  cfs := param_in[4]:
  nums :=param_in[5]:
  newsols:=[]:
  newvals:=[]:
  tsols:=[]:
  tvals:=[]:
#Identifiy those solutions (and corresponding values) for which the signs in
# vals are ambiguous
  for i from 1 to nops (sols) do 
      if member(true, 
          map(l->if ((l[1]<=0 and l[2]>=0) and not(l[1]=0 and l[2]=0)) then true
          fi, vals[i])) then 
          tvals := [op(tvals), vals[i]]:
          tsols := [op(tsols), sols[i]]:
      else 
          newvals:=[op(newvals), vals[i]]:
          newsols:=[op(newsols), sols[i]]:
      end if;
  end do;
#Ici on commence par calculer les parametrisations de F union cstr[i],i=1..#cstr
#Penser a regarder si on a bien la meme forme separante
   lparam := []:
   for i from 1 to nops(cstr) do 
       newvars, param := Parametrization([op(F), cstr[i]], vars, opts union
       {"param" = 2});
       if param[1]<>1 then #Intersection is not empty 
         if param[1] <> lin_form then 
#if the linear form changed, some elimination computations are needed
             if not(indets(lin_form) subset indets(vars)) then 
               pvars:=vars;
               vvar:=[op(indets(lin_form) minus indets(vars))][1];
               pvars:=[op(pvars), vvar];
               new_lin_form:=lin_form;
             else
               if nops(indets(lin_form)) <> 1 then
                 lprint(args);
                 error "Bug detected";
               end if;
               pvars:=remove(member, vars, indets(lin_form));
               vvar:=[op(indets(lin_form))][1];
               pvars:=[op(pvars), vvar];
               new_lin_form:=NULL;
             end if;
             gb:=MSolveGroebner([op(F),
             cstr[i],new_lin_form],0,[op(pvars),vvar], 
                                 opts union {"elim"=nops(pvars)}):
             param[2] := gb[1]:
         end if;
         g:=gcd(elim, param[2]):
#Identify those solutions in tsols which cancel g (and hence cstr[i])
         idsols:=SelectVanishingSols(g, lin_form, tsols):
#Update tvals accordingly
         for j from 1 to nops(idsols) do 
             tvals[idsols[j]][idx[i]]:=[0, 0]:
         end do;
         lparam:=[op(lparam), param];
       end if;
   end do;
   newtvals:=[]:
   newtsols:=[]:
#Identify those solutions in tsols for which the corresponding values in tvals
#   have ambiguous sign
   for i from 1 to nops(tsols) do 
      if not(member(true, 
          map(l->if ((l[1]<=0 and l[2]>=0) and not(l[1]=0 and l[2]=0)) then true
          fi, tvals[i]))) then 
          newsols:=[op(newsols), tsols[i]]:
          newvals:=[op(newvals), tvals[i]]:
      else 
          newtvals:=[op(newtvals), tvals[i]]:
          newtsols:=[op(newtsols), tsols[i]]:
      end if;
   end do;
   if nops(newtvals)=0 then 
       return newsols, newvals;
   end if;
#newtsols contains solutions for which the corresponding newtvals has ambiguous
#sign but one knows the sign is not zero

   for i from 1 to nops(newtsols) do 
    boo:= false:
    lprint("newtsols -> ", evalf(newtsols[i..i]));
    newsol:=newtsols[i..i]:
    while boo = false do 
      lprint("newsol", evalf(newsol));
      interv := MakeBinaryInterval(Eval_linform(lin_form, newsol)[1]):
      prec := 2*(abs(ilog2(abs(interv[2]-interv[1])))+256);
      lprint("prec", prec);
      printf("{%d}",degree(elim));
      interv := RootFinding[RefineRoot](interv, elim):#, digits = iquo(prec,4));
      printf("[r]");
      newsol, boo := BuildSolution([elim, den, cfs, nums] , interv, vars);
      newsol:=[newsol]:
      lprint("newsol", newsol);
      for j from 1 to nops(newtvals[i]) do 
        boo:=true:
        if (newtvals[i][j][1]<=0 and newtvals[i][j][2]>=0) and
            not(newtvals[i][j][1]=0 and newtvals[i][j][2]=0) then 
          if member(j, idx, 'r') = false then 
              error "Bug in refinement";
          end if;
          val := NonZeroIntervalEvaluate(newsol[1], cstr[r]):
          printf("[e]");
          if val[1]<=0 and val[2]>=0 then
              boo := false:
              printf("[!]");
              lprint("args", args);
              break;
          else 
              newtvals[i][j]:=val:
          end if:
        end if;
      end do;
      if boo = true then 
          newsols:=[op(newsols), newtsols[i]]:
          newvals:=[op(newvals), newtvals[i]]:
      end if;
    end do;
   end do;
   gc();
   return newsols, newvals;

end proc;

#Input.
#F: list of polynomials with rational coefficients
#vars: list of variables
#G: list of polynomials with rational coefficients
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
# when G=[g1, ..., gs] is not empty
# [vars[1] = [a1, b1], ..., vars[n] = [an,bn]], [list of intervals containing
# the values taken by G at the solutions encoded by sols]
#if "output"=1 is part of the third (optional) argument output format is 
# [ vars[1] = (a1+b1)/2, ..., vars[n] = (an+bn)/2 ]
MSolveRealRoots:=proc(F::depends(list(polynom(rational,vars))),
vars::list(name), G, opts:={})
local results, dim, fname1, fname2, verb, lsols, nl, i, j, output, str, sols,
prec, eqs, nvars, lb, vals, newvars, intervals, idx, cstr, lparam,
OldDigits;

   verb:=0:
   if type(subs(opts, "verb"), integer) then 
       verb:=subs(opts, "verb"):
   end if;
   if verb = 1 then 
       verb:=1:
   end if;
   if member(0, map(degree, F)) then return [-1,[]]; fi;
   output := subs(opts, "output");
   if nops(G) > 0 and output = 1 then 
       error "mid points cannot be returned with extra constraints: change your options";
   end if;
   OldDigits:=Digits:
   #Digits:=max(10,5*iquo(max(map(degree, map(expand,[op(F), op(G)])))*max(map(ilog2, 
   #[seq(coeffs(i), i in map(expand,[op(F), op(G)]))])), 256));
   if type(subs(opts, "prec"),string) then 
       eqs:=sort(F, (a, b)->degree(a) <= degree(b))[1..min(nops(vars),
       nops(F))]:
   Digits:=max(10,
   max(map(degree,map(expand,[op(eqs),op(G)])))+iquo(max(map(ilog2,[seq(coeffs(i), i in
   map(expand, [op(eqs), op(G)]))])),2));
   end if;
   
   if verb >= 1 then
       printf("->{");
   end if;
   newvars, lparam, lsols := Parametrization(F, vars, opts union {"param" = 1}):
   if verb >= 1 then
       printf("}");
   end if;
   Digits:=OldDigits;
   if lsols[1] = -1 then 
       return lsols;
   end if;
   if lparam=[0] then 
       return [1];
   end if;

   #nl = number of lists of solutions
   nl := lsols[1]:
   sols:=[]:
   for i from 1 to nl do
       sols:=[op(sols), op(map(_p->[seq(newvars[j] = _p[j], j=1..nops(vars))], 
                               lsols[i+1]))];
   od:
   if output=1 then
      sols := map(_p->map(_c->lhs(_c)=(rhs(_c)[1]+rhs(_c)[2])/2, _p), sols);
      return [0, sols];
   end if;
   intervals := map(pt->map(rhs, pt), sols):
   lb := map(_p->map(l->evalb((l[1] <= 0 and 0 <= l[2]) and (l[1]<>0 or l[2]<>0)), _p), 
   intervals);

   if nops(map(l->select(member, true, l),lb)) > 0 then 
       sols := SplitAndRefinePerCoordinates(lparam, sols, lb,
       newvars[1..nops(vars)]);
   end if;
   if nops(G)=0 then
       return [0, sols];
   end if;
   vals := [seq([seq(NonZeroIntervalEvaluate(sols[i], G[j]), 
   j=1..nops(G))], i=1..nops(sols))];
   idx:={op(
        map(lv->seq(if (lv[j][1]<=0 and lv[j][2]>=0) and not(lv[j][1]=0 and
        lv[j][2]=0) then j else fi, j=1..nops(G)),
        vals))};
   idx:=sort(convert(idx, list));
   if nops(idx)=0 then 
      return [0, sols, vals]:
   end if;
   cstr := map(i->G[i], idx);
   if nops(cstr) > 0 then 
       sols, vals := SplitAndRefinePerConstraints(lparam, sols,
       vals, F, newvars[1..nops(vars)], cstr,  idx, opts);
       return [0, sols, vals];
   end if;
   return [0, sols, vals];
end proc:

MSolveParam:=proc(F, fc, vars, opts:={})
local res,cfs,i,dim,den,newvars,nvars,param,v;

   res:=MSolveGroebner(F,fc,vars, opts union {"param"=2});

   if nops(res) = 0 then
      error "There has been an issue in msolve computation"
   end if;

   dim := res[1];
   if dim = -1 then error "empty set" fi;
   if dim > 0 then error "not 0-dimensional" fi;

   res:=res[2];
   nvars:=res[2];
   newvars:=res[4];
   v:=newvars[-1];
   if not member(v,vars) then v:=_T fi; # use a "clean" global name
   res:=res[6];
   if res[1]<>1 then error "several output parametrizations" else res:=res[2] fi;
   param, den, cfs, res := ExtractParametrization(res, v);
   (**
   param:=ReadPolynomial(res[1],v);
   den:=ReadPolynomial(res[2],v);
   cfs:=map2(op,-1,res[3]);
   res:=map2(op,1,res[3]);
   res:=map(ReadPolynomial,res,v);
   [param,[seq(newvars[i]=normal(-res[i]/(cfs[i]*den)),i=1..nvars-1)]]
   **)
   [param,[seq(newvars[i]=normal(-res[i]/(cfs[i]*den)),i=1..nvars-1)]]
end:

ReadPolynomial:=proc(L,var)
local i;
   add(L[2][i]*var^(i-1),i=1..L[1]+1)
end:

end module:

libname:=savelibname,libname:
ssystem(cat("mkdir -p ", savelibname)):
ssystem(cat("rm ", mladirname)):
march(`create`, mladirname);
savelib(`MSolve`);


# #Input data
(**

 F:=[x1+2*x2+2*x3+2*x4+2*x5+2*x6-1,
 x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2-x1,
 2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5+2*x5*x6-x2,
 x2^2+2*x1*x3+2*x2*x4+2*x3*x5+2*x4*x6-x3,
 2*x2*x3+2*x1*x4+2*x2*x5+2*x3*x6-x4,
 x3^2+2*x2*x4+2*x1*x5+2*x2*x6-x5];
 vars:=[x1,x2,x3,x4,x5,x6];

#Usage

 sols:=MSolve:-MSolveRealRoots(F,vars,[]);

#Other example 

F:=[x+2*y+2*z-1,x^2+2*y^2+2*z^2-x,2*x*y+2*y*z-y]:
S1:=MSolve:-MSolveRealRoots(F, [x,y,z], [x, x+z-2/3, y, x-1/3, y+x-1, y+z-1/3]);

**)
