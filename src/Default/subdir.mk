################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../ADDAmain.c \
../CalculateE.c \
../GenerateB.c \
../Romberg.c \
../calculator.c \
../chebyshev.c \
../comm.c \
../crosssec.c \
../debug.c \
../fft.c \
../interaction.c \
../io.c \
../iterative.c \
../linalg.c \
../make_particle.c \
../matvec.c \
../memory.c \
../mt19937ar.c \
../oclcore.c \
../oclmatvec.c \
../param.c \
../prec_time.c \
../sinint.c \
../somnec.c \
../timing.c \
../vars.c 

OBJS += \
./ADDAmain.o \
./CalculateE.o \
./GenerateB.o \
./Romberg.o \
./calculator.o \
./chebyshev.o \
./comm.o \
./crosssec.o \
./debug.o \
./fft.o \
./interaction.o \
./io.o \
./iterative.o \
./linalg.o \
./make_particle.o \
./matvec.o \
./memory.o \
./mt19937ar.o \
./oclcore.o \
./oclmatvec.o \
./param.o \
./prec_time.o \
./sinint.o \
./somnec.o \
./timing.o \
./vars.o 

C_DEPS += \
./ADDAmain.d \
./CalculateE.d \
./GenerateB.d \
./Romberg.d \
./calculator.d \
./chebyshev.d \
./comm.d \
./crosssec.d \
./debug.d \
./fft.d \
./interaction.d \
./io.d \
./iterative.d \
./linalg.d \
./make_particle.d \
./matvec.d \
./memory.d \
./mt19937ar.d \
./oclcore.d \
./oclmatvec.d \
./param.d \
./prec_time.d \
./sinint.d \
./somnec.d \
./timing.d \
./vars.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


