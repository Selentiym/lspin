#!/bin/bash
#скачиваем и распаковываем LAPACK
mkdir ../buildLapack
cd ../buildLapack
wget http://www.netlib.org/lapack/lapack-3.7.0.tgz
tar zxvf lapack-3.7.0.tgz
cd lapack-3.7.0
#Используем makefile по умолчанию, но можно его менять, если нужно что-то особенное
mv make.inc.example make.inc
#Собираем BLAS. После этого в папке lapack-3.7.0 должен появиться файл librefblas.a
make blaslib
#Собираем LAPACK. Результат liblapack.a
make lapacklib
mv librefblas.a ../../lspin/libblas.a
mv liblapack.a ../../lspin/liblapack.a