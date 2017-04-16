################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../constants.f90 \
../integrate.f90 \
../laguerre.f90 \
../lapack_prb.f90 \
../lspinors.f90 \
../main.f90 \
../matrix.f90 \
../test.f90 

OBJS += \
./constants.o \
./integrate.o \
./laguerre.o \
./lapack_prb.o \
./lspinors.o \
./main.o \
./matrix.o \
./test.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

constants.o: ../constants.f90

integrate.o: ../integrate.f90 constants.o laguerre.o

laguerre.o: ../laguerre.f90

lapack_prb.o: ../lapack_prb.f90

lspinors.o: ../lspinors.f90 constants.o laguerre.o

main.o: ../main.f90

matrix.o: ../matrix.f90

test.o: ../test.f90 constants.o integrate.o laguerre.o lspinors.o matrix.o


