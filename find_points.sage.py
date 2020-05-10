
# This file was *autogenerated* from the file find_points.sage
from sage.all_cmdline import *   # import sage library

_sage_const_3 = Integer(3); _sage_const_2 = Integer(2); _sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_4 = Integer(4)
import sys
m=_sage_const_4  #it is actually the value of m+1
n=int(sys.argv[_sage_const_1 ])
l=int(sys.argv[_sage_const_2 ])
nb_inst=int(sys.argv[_sage_const_3 ])
K = GF(_sage_const_2 **n, names=('a',)); (a,) = K._first_ngens(1)
R = PolynomialRing(K, names=('y',)); (y,) = R._first_ngens(1)
E = EllipticCurve(K,[_sage_const_1 ,_sage_const_1 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ])

IRR=K.modulus()
lst = IRR.coefficients(sparse=False)
strR= ''.join(str(e) for e in lst)

for j in range(_sage_const_0 , nb_inst):
    P=[None] * m
    Y=[None] * m
    for i in range(_sage_const_0 ,m-_sage_const_1 ):
        ok=_sage_const_0 
        while ok == _sage_const_0 :
            listCoefs = []
            for a in range(_sage_const_0 ,l):
                listCoefs.insert(_sage_const_0 ,randint(_sage_const_0 ,_sage_const_1 ))
            X=K(listCoefs)
            f=y**_sage_const_2 +X*y-X**_sage_const_3 -X**_sage_const_2 -_sage_const_1 
            roots=f.roots()
            if len(roots) > _sage_const_0 :
                ok = _sage_const_1 
                Y=roots[_sage_const_0 ][_sage_const_0 ]
                P[i] = E(X, Y)
    R=P[_sage_const_0 ]
    for i in range(_sage_const_1 ,m-_sage_const_1 ):
        R=R+P[i]
    X1=P[_sage_const_0 ][_sage_const_0 ]
    X2=P[_sage_const_1 ][_sage_const_0 ]
    X3=P[_sage_const_2 ][_sage_const_0 ]
    x__1=R[_sage_const_0 ]
    e1=X1+X2+X3
    e2=X1*X2+X2*X3+X1*X3
    e3=X1*X2*X3
    f3=x__1**_sage_const_4 +e1**_sage_const_4 +e3**_sage_const_4 +e2**_sage_const_4 *x__1**_sage_const_4 +e3**_sage_const_3 *x__1+e3*e2**_sage_const_2 *x__1**_sage_const_3 +e3*e1**_sage_const_2 *x__1+e3*x__1**_sage_const_3 +e1**_sage_const_2 *e3**_sage_const_2 *x__1**_sage_const_2 +e3**_sage_const_2 *x__1**_sage_const_4 +e3**_sage_const_2 +e2**_sage_const_2 *x__1**_sage_const_2 
    #print e1, " " , e2, " " , e3
    #print X1.polynomial()
    lst = X1.polynomial().coefficients(sparse=False)
    while(len(lst) < l):
        lst.append(_sage_const_0 );
    strX1= ''.join(str(e) for e in lst)
    
    lst = X2.polynomial().coeffs()
    while(len(lst) < l):
        lst.append(_sage_const_0 );
    strX2= ''.join(str(e) for e in lst)

    lst = X3.polynomial().coeffs()
    while(len(lst) < l):
        lst.append(_sage_const_0 );
    strX3= ''.join(str(e) for e in lst)

    lst = x__1.polynomial().coeffs()
    while(len(lst) < n):
        lst.append(_sage_const_0 );
    strXr= ''.join(str(e) for e in lst)
    
    cmd="./script_benchmarks.sh "+str(n)+" "+str(l)+" "+strR+" "+strXr+" "+strX1+"-"+strX2+"-"+strX3+" S;"
    print cmd
for j in range(_sage_const_0 , nb_inst):
    listCoefs = []
    for a in range(_sage_const_0 ,n):
        listCoefs.insert(_sage_const_0 ,randint(_sage_const_0 ,_sage_const_1 ))
    strXr= ''.join(str(e) for e in listCoefs)
    cmd="./script_benchmarks.sh "+str(n)+" "+str(l)+" "+strR+" "+strXr+" "+"UNSAT"+" U;"
    print cmd
