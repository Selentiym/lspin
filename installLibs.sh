#!/bin/bash
#скачиваем и распаковывем LAPACK
mkdir ../buildLapack
cd ../buildLapack
wget http://www.netlib.org/lapack/lapack-3.7.0.tgz
tar zxvf lapack-3.7.0.tgz
cd lapack-3.7.0
#Используем makefile по умолчанию
mv make.inc.example make.inc
#Собираем BLAS. Появится librefblas.a
make blaslib
#Собираем LAPCK. Появится liblapack.a
make lapacklib
mv librefblas.a ../../lspin/libblas.a
mv liblapack.a ../../lspin/liblapack.a
