"""
17/02/2021
Adapted from the maple interface written by Huu Phuoc Le and Jorge Garcia-Fontan.
"""

import os

def ToMSolve(F, finput="/tmp/in.ms"):
    """Convert a system of sage polynomials into a msolve input file.

    Inputs :
    F (list of polynomials): system of polynomial to solve
    finput (string): name of the msolve input file.

    """
    A = F[0].parent()
    assert all(A1 == A for A1 in map(parent,F)),\
            "The polynomials in the system must belong to the same polynomial ring."

    variables, char = A.variable_names(), A.characteristic()
    s = (", ".join(variables) + " \n"
            + str(char) + "\n")

    B = A.change_ring(order = 'degrevlex') 
    F2 = [ str(B(f)).replace(" ", "") for f in F ]
    if "0" in F2:
        F2.remove("0")
    s += ",\n".join(F2) + "\n"

    fd = open(finput, 'w')
    fd.write(s)
    fd.close()


def FormatOutputMSolveOnlySolutions(foutput):
    """Convert a msolve output file into solutions

    Inputs :
    foutput (string): name of the msolve output file

    Output :
        The set of solutions computed by msolve.

    """
    f = open(foutput,'r')
    s = f.read()
    s = s.replace("\n","").replace(":","")
    R = sage_eval(s)
    intervals = R[1][1]
    S   =   []
    if len(intervals) > 0:
        nvars   =   len(intervals[0])
        for sol in intervals:
            s = []
            for i in range(nvars):
                s.append((sol[i][0]+sol[i][1])/2)
            S.append(s)
    return S

def FormatOutputMSolve(foutput):
    """Convert a msolve output file into a rational parametrization 

    Inputs :
    foutput (string): name of the msolve output file

    Output :
        A rational parametrization of the zero-dimensional ideal describing
    the solutions. Note : p[i] and c[i] stand for the (i+1)-th coordinate.

    """
    f = open(foutput,'r')
    s = f.read()
    s = s.replace("\n","").replace(":","")
    R = sage_eval(s)
    A.<t> = QQ[]
    # dimension
    dim = R[0]
    if dim > 0:
        return None, None, A(-1), None, None, None, None

    # parametrization
    nvars       = R[1][1]
    qdim        = R[1][2]
    varstr      = R[1][3]
    linearform  = R[1][4]
    elim        = R[1][5][1][0]
    den         = R[1][5][1][1]
    polys       = R[1][5][1][2]
    # solutions
    intervals   = R[2][1]

    #  nvars, degquot, deg = L[1], L[2], L[5][0]
    #  varstr      =   L[3]
    #  linearform  =   L[4]

    if len(elim) > 0:
        pelim = A(elim[1])
    else:
        return None, None, A(-2), None, None, None, None

    pden, p, c = A(1), [], []
    if qdim > 0:
        pden = A(den[1])
        for l in polys:
            p.append(A(l[0][1]))
            c.append( l[1] )

    S   =   []
    if len(intervals) > 0:
        for sol in intervals:
            s = []
            for i in range(nvars):
                s.append((sol[i][0]+sol[i][1])/2)
            S.append(s)
    return [varstr, linearform, pelim, pden, p, c, S]

def GetRootsFromMSolve(foutput, param):
    """Compute rational approximation roots from an msolve output file
    The rational number is chosen in the isolating interval to be
    the smallest in bitsize.

    Inputs :
    foutput (string): name of the msolve output file

    Output :
        b (integer): error code
    Qroots : list of rationals approximations of the roots

    """

    if param == 1:
        varstr, linearform, elim, den, p, c, sols = FormatOutputMSolve(foutput)
        if elim.degree() == 0:
            return elim, [], []
        return 0, [varstr, linearform, elim, den, p, c], sols
    else:
        sols = FormatOutputMSolveOnlySolutions(foutput)
        return 0, [], sols


def MSolveRealRoots(F, fname1="/tmp/in.ms", fname2="/tmp/out.ms",
        mspath="../binary/msolve", v=0, p=1):
    """Computes the a rational approximation of the real roots
    of a system of sage polynomials using msolve. 

    Inputs :
    F (list of polynomials): system of polynomials to solve
    fname1 (string): complete name of the msolve input file used
    fname2 (string): complete name of the msolve output file used
    mspath (string): path to the msolve binary
    v (in [0,1,2]): level of msolve verbosity

    Output :
        sols (list of lists): list of rational approximation roots to the system
    represented by F.

    """

    ToMSolve(F, fname1)

    os.system(mspath +" -v " + str(v) +" -P " + str(p) +  " -f " + fname1 + " -o " + fname2)

    b, param, sols = GetRootsFromMSolve(fname2,p)

    if b == -1:
        print("System has infinitely many complex solutions")
        return []
    if b == -2:
        print("System not in generic position. You may add to your system")
        print("a random linear form of your variables and a new variable")
        return []
    #New(s) variable(s) may have been introduced at the end, for genericity purposes.	
    n = len(F[0].parent().gens())
    if p == 0:
        return [ s[:n] for s in sols ]
    else:
        return param, [ s[:n] for s in sols ]

