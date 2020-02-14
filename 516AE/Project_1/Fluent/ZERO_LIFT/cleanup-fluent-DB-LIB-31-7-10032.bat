echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v192\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v192\fluent\ntbin\win64\tell.exe" DB-LIB-31-7 50533 CLEANUP_EXITING
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 1896) 
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 11104) 
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 1028) 
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 12048) 
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 10032) 
if /i "%LOCALHOST%"=="DB-LIB-31-7" (%KILL_CMD% 8172)
del "P:\AE516\Project_1\Fluent\ZERO_LIFT\cleanup-fluent-DB-LIB-31-7-10032.bat"
