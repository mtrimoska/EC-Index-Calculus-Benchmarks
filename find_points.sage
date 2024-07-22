import sys
m=4 #it is actually the value of m+1
n=int(sys.argv[1])
l=int(sys.argv[2])
nb_inst=int(sys.argv[3])
K.<a> = GF(2^n)
R.<y> = PolynomialRing(K)
E = EllipticCurve(K,[1,1,0,0,1])

IRR=K.modulus()
lst = IRR.coefficients(sparse=False)
strR= ''.join(str(e) for e in lst)

for j in range(0, nb_inst):
    P=[None] * m
    Y=[None] * m
    for i in range(0,m-1):
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
    R=P[0]
    for i in range(1,m-1):
        R=R+P[i]
    X1=P[0][0]
    X2=P[1][0]
    X3=P[2][0]
    x__1=R[0]
    e1=X1+X2+X3
    e2=X1*X2+X2*X3+X1*X3
    e3=X1*X2*X3
    f3=x__1^4+e1^4+e3^4+e2^4*x__1^4+e3^3*x__1+e3*e2^2*x__1^3+e3*e1^2*x__1+e3*x__1^3+e1^2*e3^2*x__1^2+e3^2*x__1^4+e3^2+e2^2*x__1^2
    #print e1, " " , e2, " " , e3
    lst = X1.polynomial().coefficients(sparse=False)
    while(len(lst) < l):
        lst.append(0);
    strX1= ''.join(str(e) for e in lst)
    lst = X2.polynomial().coefficients(sparse=False)
    while(len(lst) < l):
        lst.append(0);
    strX2= ''.join(str(e) for e in lst)
    lst = X3.polynomial().coefficients(sparse=False)
    while(len(lst) < l):
        lst.append(0);
    strX3= ''.join(str(e) for e in lst)
    lst = x__1.polynomial().coefficients(sparse=False)
    while(len(lst) < n):
        lst.append(0);
    strXr= ''.join(str(e) for e in lst)

    cmd="./script_benchmarks.sh "+str(n)+" "+str(l)+" "+strR+" "+strXr+" "+strX1+"-"+strX2+"-"+strX3+" S;"
    print(cmd)
for j in range(0, nb_inst):
    listCoefs = []
    for a in range(0,n):
        listCoefs.insert(0,randint(0,1))
    strXr= ''.join(str(e) for e in listCoefs)
    cmd="./script_benchmarks.sh "+str(n)+" "+str(l)+" "+strR+" "+strXr+" "+"UNSAT"+" U;"
    print(cmd)
