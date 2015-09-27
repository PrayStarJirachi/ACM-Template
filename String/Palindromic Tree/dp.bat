@echo off
:again
call data
call std > output.txt
call prt > outputbc.txt
fc output.txt outputbc.txt
if errorlevel 1 echo WA! && pause > nul & exit
goto again