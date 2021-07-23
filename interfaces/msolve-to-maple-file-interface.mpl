
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

# GetRootsFromMSolve:=proc(l)
# local _Z, e, p, c, sols, vals, i, sols_and_vals, realroots, OldDigits, mysols;
#   e, p, c := FormatOutputMSolve(l, _Z);
#   if degree(e) = 0 then #positive dimension (-1) or computation failed (-2)
#   return e, [];
#   fi;
#   OldDigits:=Digits:
#   Digits:=20:
#   sols_and_vals:=[RootFinding[Isolate](e, constraints=p, output=midpoint)];
#   sols:=sols_and_vals[1]:

#   mysols:=map(rhs, sols):

#   vals:=sols_and_vals[2]:
#   realroots:=[]:
#   for i from 1 to nops(sols) do
#     realroots:=[op(realroots),
# #         [seq(-rhs(vals[i][j])/rhs(vals[i][1]), j = 2..nops(p)), rhs(sols[i])]
#          [seq(-rhs(vals[i][j])/(c[j]*rhs(vals[i][1])), j = 2..nops(p)), rhs(sols[i])]
#          ]:
#   od:
#   Digits:=OldDigits:
#   return 0, realroots;
# end proc:


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
   str := cat("sed -i -e ':a' -e 'N' -e '$!ba' -e 's/\\\\\\n//g' ", fname):
   system(str):
end proc:

MSolveRealRoots:=proc()
local i, results, get_param, param, b, sols, str, F, vars, fname1, fname2, mspath;
  if nargs=6 then
    F:=args[1]:
    vars:=args[2]:
    mspath:=args[3]:
    fname1:=args[4]:
    fname2:=args[5]:
    get_param:=args[6]:
  elif nargs=5 then
    F:=args[1]:
    vars:=args[2]:
    mspath:=args[3]:
    fname1:=args[4]:
    fname2:=args[5]:
    get_param:=0:
  elif nargs=3 then
    F:=args[1]:
    vars:=args[2]:
    mspath:=args[3]:
    fname1:="/tmp/in.ms":
    fname2:="/tmp/out.ms":
    get_param:=0:
  else
    F:=args[1]:
    vars:=args[2]:
    mspath:="../binary/msolve":
    fname1:="/tmp/in.ms":
    fname2:="/tmp/out.ms":
    get_param:=0:
  fi:
  ToMSolve(F, 0, vars, fname1):
  str := cat(mspath," -v 1 -P ", get_param, " -f ", fname1," -o ", fname2):
  system(str):
  read(fname2):
  #The output of msolve should be modified to identify issues
  #(ideal not 0-dim, last variable not in generic position, etc.).
  results:=%:
  if get_param = 0 then
    sols:=map(s->[seq(vars[i]=s[i],i=1..nops(vars))], results):
    return sols;
  fi;
  b, param :=GetRootsFromMSolve(results[1]);
  if b = -1 then
    printf("System has infinitely many complex solutions\n");
    return [], [];
  fi;
  if b = -2 then
    printf("System not in generic position. You may add to your system\n");
    printf("A random linear form of your variables and a new variable");
    return [], [];
  fi;
  sols:=map(s->[seq(vars[i]=s[i],i=1..nops(vars))], results[2]):
  return param, sols;
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
