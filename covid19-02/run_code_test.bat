REM run_code_test.bat
REM -----------------
REM コードが動くことのテスト
@echo off

REM SET DAY=503-07-30

SET DAY=%date%
for  /f "delims=/, tokens=1,2,3" %%a in ("%date%") do SET DAY=%%a-%%b-%%c

set cs=1.0
SET A=%time%
for %%c in (%cs%) do mpiexec -np 6 metapop_da.exe 7 beta-ini.csv -inffile init_inf_tokyo_osaka.csv -outdir %DAY%-test-odset-conn%%c -t-upto 14 -conn %%c

SET B=%time%
ECHO %B% %A%
