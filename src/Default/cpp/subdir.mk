################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../cpp/fft_execute.cpp \
../cpp/fft_kernelstring.cpp \
../cpp/fft_setup.cpp 

OBJS += \
./cpp/fft_execute.o \
./cpp/fft_kernelstring.o \
./cpp/fft_setup.o 

CPP_DEPS += \
./cpp/fft_execute.d \
./cpp/fft_kernelstring.d \
./cpp/fft_setup.d 


# Each subdirectory must supply rules for building sources it contributes
cpp/%.o: ../cpp/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


