@echo off

set COM=program_files
set RUN=%~1
set BIN=%COM%\bin\%~2\badchimp.exe
set GEO=%COM%\geo\%~3
set OUT=%RUN%\out
set INP=%RUN%\input.dat
set OPT=%RUN%\options.dat
set BAS=%COM%\basis
set DAT=%COM%\data

:: check if out-folder exists, and increment name if true

set /a "c=1"
:while
if exist %OUT% (
  set OUT=%RUN%\out-%c%
  set /a "c=c + 1"
  goto :while
)

:: extract process vector from input.dat and calculate total number of processes

::for /f "tokens=2,3,4 delims= " %%A in ('findstr /r /c:" *num_proc.*" %~1\input.dat') do set /a NUM=%%A*%%B*%%C


for /f "tokens=2,3,4 delims= " %%A in ('findstr /r /c:" *num_proc.*" %INP%') do set /a NUM=%%A*%%B*%%C



:: increase cmd window

mode con:cols=130 lines=500

mpiexec -n %NUM% %BIN% -basis %BAS% -data %DAT% -out %OUT% -geo %GEO% -input %INP% -options %OPT%
