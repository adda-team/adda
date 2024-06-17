@echo off
:: Compiles ADDA in seq mode using Intel compiler for Windows (which uses some of the tools from Visual Studio).
:: Note that Intel compilers behave very differently (in terms of options and predefined macros) on Windows than
:: on Linux/macOS. Thus, this file is designed only for testing Intel compiler, but not for production compilation
:: in various modes. For the latter use the main Makefile with MinGW64 toolchain.
::
:: Copyright (C) ADDA contributors
:: GNU General Public License version 3
::
:: Tested with Intel Compilers 2021.2 (v.2021.1) on top of Visual Studio Build Tools 2019 (v.16.9.0) on Windows 10 x64
:: Both "Classic" and "oneAPI" (LLVM) compilers were tested. Note that Visual Studio Build Tools should include
:: Desktop development with C++. More recently, it was tested with VS 2022 (v.19.40) and Intel oneAPI compilers 
:: version 202401.2. They were obtained from:
:: https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2022
:: https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-undefined
:: https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#inpage-nav-6-1
::
:: To perform the compilation:
:: 1) Find and click "x64 Native Tools Command Prompt for VS ..." in Windows Run menu
:: 2) (only once) Additionally set up FFTW, assuming you already downloaded it - see
::    https://github.com/adda-team/adda/wiki/InstallingFFTW3
::    a) Set the path to it below (fftw-path)
::    b) cd %fftw_path%   (into its directory)
::       C:   (change drive letter if needed)
::       lib /def:libfftw3-3.def
:: 3) Set environment for Intel compilers, by running
::    <Intel>\compiler\latest\env\vars.bat (<Intel> is the install directory)
:: 4) cd into the script directory (.../adda/src)
::    (for repeated use, you may combine steps 3 and 4 into a batch file placed into a directory, where you end on step 1)
:: 5) Check other settings in this script below (compilation flags, use of Fortran, source files,...)
:: 6) Run this script

setlocal

:: comment out the following line to avoid Fortran sources (if any problems appear)
set with_fortran=true
:: uncomment the following to use Classic Intel compilers (but they are deprecated)
::set intel_type=classic
:: uncomment the following to compile in debug mode
::set debug=true
set fftw_path=C:\local\fftw-3.3.5

:: many diagnostics on Windows are due to failure to understand GCC attributes like noreturn and unused
if "%intel_type%"=="classic" (
  set cc=icl
  set cwarn=-O1 -W5 -Wcheck -D_CRT_SECURE_NO_WARNINGS^
  -Qdiag-disable:193,593,869,1011,1418,1419,1572,11074,11075,13000,13046
) else (
  set cc=icx
  set cwarn=-W4 -D_CRT_SECURE_NO_WARNINGS -Wno-return-type -Wno-unused-parameter -Wno-unused-function^
  -Wno-sometimes-uninitialized -Wno-implicit-const-int-float-conversion -Wno-misleading-indentation
)
:: It may be possible to use Visual Studio C++ compiler (cl) in the future (using, e.g., "/std:c11"), 
:: but even VS 2022 doesn't yet fully support complex numbers - see
:: https://learn.microsoft.com/en-us/cpp/c-runtime-library/complex-math-support?view=msvc-170
set cstd=-Qstd=c99
if "%debug%"=="true" (
  set opt=-O1
  set defs=-DDEBUG
) else (
  set opt=-Ofast
  set defs=
  set cwarn=-w
)

:: C sources are everything except prec_time.c and those with "ocl" or "_test" in its name
set cfiles=
for /f %%G in ('dir /b *.c ^| find /V "ocl" ^| find /V "_test" ^| find /V "prec_time.c"') do (
  call set cfiles=%%cfiles%% ..\%%G
)

if "%with_fortran%"=="true" (
  if "%intel_type%"=="classic" (
    set cf=ifort
  ) else (
    set cf=ifx
  )
  set ffiles=
  for /f %%G in ('dir /b fort\*.f fort\*.f90 ^| find /V "cfft99D.f"') do (
    call set ffiles=%%ffiles%% ..\fort\%%G
  )
rem The following is required, since ifx on Windows makes function names in object files differently than on other OS
rem Starting from Fortran 2003 it is recommended to use ISO_C_BINDING, but it is not a perfect solution for legacy code
rem For now ADDA requires only a few single cross-language function call, so we replace the naming in C sources
rem P.S. comments starting with :: do not work inside ()
  call set defs=%%defs%% -Dpropaespacelibreintadda_=PROPAESPACELIBREINTADDA -Dbjndd_=BJNDD
) else (
  call set defs=%%defs%% -DNO_FORTRAN
)

cd seq
:: to avoid influence of previous compilations (since we do not have a proper Makefile)
del *.obj
%cc% -c -I"%fftw_path%" %cstd% %opt% %cwarn% %defs% %cfiles%
if "%with_fortran%"=="true" %cf% -c %opt% %ffiles%
%cc% -o adda *.obj "%fftw_path%\libfftw3-3.lib"
