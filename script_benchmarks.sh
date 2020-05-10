#!/bin/bash
n=$1
l=$2
r=$3
x=$4
res=$5
sat=$6
idFile=id-n${n}l${l}

cd Weil_descent
gcc *.c -o weill > /dev/null 2>&1
./weill -a 2 -n $n -l $l -r $r -x $x -o dimacs
gcc *.c -o weill > /dev/null 2>&1
./weill -a 22 -n $n -l $l -r $r -x $x -o dimacs
gcc *.c -o weill > /dev/null 2>&1
./weill -a 3 -n $n -l $l -r $r -x $x -o dimacs
cd out
./assemble_dimacs.sh
cp weil.dimacs ../../
./assemble_anf.sh
cp weil.anf ../../
cd ../../
if [ ! -f $idFile ]
then
echo 1 > $idFile
fi
id=$(head -n 1 $idFile)
newId=$(($id+1))
echo $newId > $idFile
dimacsFile=Xn${n}l${l}-${id}-${sat}.dimacs
anfFile=Xn${n}l${l}-${id}-${sat}.anf
infoFile=INFOn${n}l${l}-${id}-${sat}.dimacs
echo "${n} ${l}" > benchmarks/$infoFile
echo "${r}" >> benchmarks/$infoFile
echo "${x}" >> benchmarks/$infoFile
echo "${sat}" >> benchmarks/$infoFile
echo "${res}" >> benchmarks/$infoFile
cat weil.dimacs > benchmarks/$dimacsFile
cat weil.anf > benchmarks/$anfFile

cd Weil_descent
gcc *.c -o weill > /dev/null 2>&1
./weill -a 2 -n $n -l $l -r $r -x $x
gcc *.c -o weill > /dev/null 2>&1
./weill -a 22 -n $n -l $l -r $r -x $x
gcc *.c -o weill > /dev/null 2>&1
./weill -a 3 -n $n -l $l -r $r -x $x
cd out
./assemble_grobner.sh
grobnerFile=n${n}l${l}-${id}-${sat}.in
cp magma.in ../../benchmarks/$grobnerFile

cd ../../benchmarks
nbVars=$(head -n 1 $dimacsFile | cut -d' ' -f 3)
newNbVars=$(../XORtoCNF.sh $nbVars $dimacsFile)
newNbVars=$(($newNbVars-1))
newNbLines=$(wc -l < F_3_norm.dimacs)
newNbLines=$(($newNbLines-1))
cnfFile=n${n}l${l}-${id}-${sat}.dimacs
echo "p cnf $newNbVars $newNbLines" > $cnfFile
tail -n+2 F_3_norm.dimacs >> $cnfFile

cd ..
