# Extended Euclidean for the symbolic ring (assumes r is char 2 and with constant coefficients, f is with variables coefficients)
def SR_mod(f, r, X):
    LTr = r.lt()
    cpt = 10
    while f.degree() >= r.degree():
        LCf = f.lc()
        # Do f = f - f.lt()/r.lt()*r without using division
        f = f - LCf*X^(f.degree() - r.degree())*r
        cpt -= 1
    return f

# Did not find how to do this properly in Sage - so doing it manually
def remove_powers(f, field_equations):
    for equi in field_equations:
        I = Ideal(f, equi)
        f = I.groebner_basis()[0]
    return f

def generate_SAT_instance():
    P=[None] * m 
    for i in range(0,m-1):
        # Can not use P[i] = E.random_point() because we need points where the X coordinate is of degree at most l
        ok=0
        while ok == 0:
            listCoefs = []
            for a in range(0,l):
                listCoefs.insert(0,randint(0,1))
            X=K(listCoefs)
            f=y^2+X*y-X^3-X^2-1
            roots=f.roots()
            if len(roots) > 0:
                ok = 1
                Y=roots[0][0]
                P[i] = E(X, Y)
    R=sum(P[i] for i in range(m-1))
    X1=P[0][0]
    X2=P[1][0]
    X3=P[2][0]
    Xr=R[0]
    e1=X1+X2+X3
    e2=X1*X2+X2*X3+X1*X3
    e3=X1*X2*X3
    f3=Xr^4+e1^4+e3^4+e2^4*Xr^4+e3^3*Xr+e3*e2^2*Xr^3+e3*e1^2*Xr+e3*Xr^3+e1^2*e3^2*Xr^2+e3^2*Xr^4+e3^2+e2^2*Xr^2
    assert(f3 == 0)
    return Xr.polynomial(), X1.polynomial(), X2.polynomial(), X3.polynomial()

def generate_UNSAT_instance():
    R = E.random_point()
    Xr = R[0]
    return Xr.polynomial()

def get_INFO_output(isSat, n, l, X1, X2, X3, Xr):
    output_info = ""
    if isSat:
        lst = X1.coefficients(sparse=False)
        while(len(lst) < l):
            lst.append(0);
        strX1= ''.join(str(e) for e in lst)
        
        lst = X2.coefficients(sparse=False)
        while(len(lst) < l):
            lst.append(0);
        strX2= ''.join(str(e) for e in lst)

        lst = X3.coefficients(sparse=False)
        while(len(lst) < l):
            lst.append(0);
        strX3= ''.join(str(e) for e in lst)

        lst = Xr.coefficients(sparse=False)
        while(len(lst) < n):
            lst.append(0);
        strXr= ''.join(str(e) for e in lst)
        output_info =  str(n)+" "+str(l)+"\n"+strR+"\n"+strXr+"\nS\n"+strX1+"-"+strX2+"-"+strX3+"\n"
    return output_info

def get_ANF_output_X(nbVars, nbEqs, eqs_X):
    output = "p cnf " + str(nbVars) + " " + str(nbEqs) + "\n"
    i_e = nbVarsX
    for p in eqs_X:
        i_e += 1
        output += "x T " + str(i_e)
        for m in p.monomials():
            d = m.degree()
            if d > 1:
                output += " ." + str(d)
            temp = m.dict().keys()
            for v in temp: #once
                v = list(v)
                for i in range(nbVarsX):
                    if v[i] > 0:
                        output += " " + str(i + 1)
        output += " 0\n"
    return output

def get_ANF_output_E(offset, nbVarsE, eqs_E):
    output = ""
    for p in eqs_E:
        if p.constant_coefficient() == 1:
            output += "x"
        else:
            output += "x T"
        for m in p.monomials():
            d = m.degree()
            if d > 1:
                output += " ." + str(d)
            temp = m.dict().keys()
            for v in temp: #once
                v = list(v)
                for i in range(nbVarsE):
                    if v[i] > 0:
                        output += " " + str(i + offset + 1)
        output += " 0\n"
    return output

import sys
m=4 #it is actually the value of m+1
n=int(sys.argv[1])
l=int(sys.argv[2])
nb_inst=int(sys.argv[3])
K.<a> = GF(2^n)
IRR=K.modulus()

nbVarsX = (m-1)*l
nbVarsE = 0
for i in range(1, m):
    nbVarsE = nbVarsE + i * l - (i - 1)
nbVars = nbVarsX + nbVarsE
RxBase = PolynomialRing(K, nbVarsX, "x_")
Rx.<x> = PolynomialRing(RxBase)
xc = RxBase.gens()
XVarsStr = str(xc).strip('(').strip(')')
ReBase = PolynomialRing(K, nbVarsE, "e_")
Re.<e> = PolynomialRing(ReBase)
ec = ReBase.gens()
eVarsStr = str(ec).strip('(').strip(')')
field_equations = [ec[i]^2 + ec[i] for i in range(nbVarsE)]

#get irreducible with coeffs over Rx
IRRList = IRR.coefficients(sparse=False)
_IRR = sum(int(IRRList[k])*x^k for k in range(n+1))


#modelisation generic part: e-X relation
_X1 = sum(xc[k]*x^k for k in range(l))
_X2 = sum(xc[k+l]*x^k for k in range(l))
_X3 = sum(xc[k+2*l]*x^k for k in range(l))

e1 = _X1 + _X2 + _X3
e2 = _X1*_X2 + _X1*_X3 + _X2*_X3
e3 = _X1*_X2*_X3
e1 = SR_mod(e1, _IRR, x) #no need
e2 = SR_mod(e2, _IRR, x) #probably no need
e3 = SR_mod(e3, _IRR, x) #probably no need


eqs_X = [(expand(e1.coefficients(sparse=False)[i])) for i in range(e1.degree()+1)] + [(expand(e2.coefficients(sparse=False)[i])) for i in range(e2.degree()+1)] + [(expand(e3.coefficients(sparse=False)[i])) for i in range(e3.degree()+1)] #without the e_i

eqs_X_str = [str(ec[i])+ " + " + str((expand(e1.coefficients(sparse=False)[i]))) for i in range(e1.degree()+1)] + [str(ec[l + i])+ " + " + str((expand(e2.coefficients(sparse=False)[i]))) for i in range(e2.degree()+1)] + [str(ec[3*l - 1 + i])+ " + " + str((expand(e3.coefficients(sparse=False)[i]))) for i in range(e3.degree()+1)]

# Get Magma generic output
output_magma_start = "R<" + XVarsStr + ", " + eVarsStr + "> := BooleanPolynomialRing(" + str(nbVars) + ", \"grevlex\"); \n B := [ " + ','.join(eqs_X_str) + "," 
output_magma_end = "]; \n I := Ideal(B); \n Variety(I);"

# Get SAT ANF generic output
output_ANF_start = get_ANF_output_X(nbVars, len(eqs_X) + n, eqs_X)

#modelisation Weil descent part
R.<y> = PolynomialRing(K)
_IRR = sum(int(IRRList[k])*e^k for k in range(n+1))
E = EllipticCurve(K,[1,1,0,0,1]) #Koblitz curve
strR= ''.join(str(e) for e in IRRList)
offset = 0
_e1 = sum(ec[k + offset]*e^k for k in range(l))
offset += l
_e2 = sum(ec[k + offset]*e^k for k in range(2*l - 1))
offset += 2*l - 1
_e3 = sum(ec[k + offset]*e^k for k in range(3*l - 2))

for j in range(0, nb_inst):
    Xr, X1, X2, X3 = generate_SAT_instance()
    XrList = Xr.coefficients(sparse=False)
    _Xr = sum(int(XrList[k])*e^k for k in range(Xr.degree()+1))
    _f3 = _Xr^4 + _e1^4 + _e3^4 + _e2^4*_Xr^4 + _e3^3*_Xr + _e3*_e2^2*_Xr^3 + _e3*_e1^2*_Xr + _e3*_Xr^3 + _e1^2*_e3^2*_Xr^2 + _e3^2*_Xr^4 + _e3^2 + _e2^2*_Xr^2
    _f3 = SR_mod(_f3, _IRR, e)
    eqs_E = [(remove_powers(expand(_f3.coefficients(sparse=False)[i]), field_equations)) for i in range(_f3.degree()+1)]
    eqs_E_str = [str(e) for e in eqs_E]

    file_name_root = "n" + str(n) + "l" + str(l) + "-" + str(j) + "-S"

    # Get INFO output
    output_info = get_INFO_output(True, n, l, X1, X2, X3, Xr)
    file_name = "benchmarks/INFO" + file_name_root + ".dimacs"
    f = open(file_name, "w")
    f.write(output_info)
    f.close()

    # Get Magma output
    eqs_str = eqs_X_str + eqs_E_str
    output_magma = output_magma_start + ','.join(eqs_E_str) + output_magma_end
    file_name = "benchmarks/" + file_name_root + ".in"
    f = open(file_name, "w")
    f.write(output_magma)
    f.close()

    # Get SAT ANF output
    output_ANF_end = get_ANF_output_E(nbVarsX, nbVarsE, eqs_E)
    output_ANF = output_ANF_start + output_ANF_end
    file_name = "benchmarks/X" + file_name_root + ".anf"
    f = open(file_name, "w")
    f.write(output_ANF)
    f.close()

        
for j in range(0, nb_inst):
    Xr = generate_UNSAT_instance()
    XrList = Xr.coefficients(sparse=False)
    _Xr = sum(int(XrList[k])*e^k for k in range(Xr.degree()+1))
    _f3 = _Xr^4 + _e1^4 + _e3^4 + _e2^4*_Xr^4 + _e3^3*_Xr + _e3*_e2^2*_Xr^3 + _e3*_e1^2*_Xr + _e3*_Xr^3 + _e1^2*_e3^2*_Xr^2 + _e3^2*_Xr^4 + _e3^2 + _e2^2*_Xr^2
    _f3 = SR_mod(_f3, _IRR, e)
    eqs_E = [(remove_powers(expand(_f3.coefficients(sparse=False)[i]), field_equations)) for i in range(_f3.degree()+1)]
    eqs_E_str = [str(e) for e in eqs_E]

    file_name_root = "n" + str(n) + "l" + str(l) + "-" + str(j) + "-U"

    # Get INFO output
    output_info = get_INFO_output(False, n, l, None, None, None, Xr)
    file_name = "benchmarks/INFO" + file_name_root + ".dimacs"
    f = open(file_name, "w")
    f.write(output_info)
    f.close()

    # Get Magma output
    eqs_str = eqs_X_str + eqs_E_str
    output_magma = output_magma_start + ','.join(eqs_E_str) + output_magma_end
    file_name = "benchmarks/" + file_name_root + ".in"
    f = open(file_name, "w")
    f.write(output_magma)
    f.close()

    # Get SAT ANF output
    output_ANF_end = get_ANF_output_E(nbVarsX, nbVarsE, eqs_E)
    output_ANF = output_ANF_start + output_ANF_end
    file_name = "benchmarks/X" + file_name_root + ".anf"
    f = open(file_name, "w")
    f.write(output_ANF)
    f.close()
