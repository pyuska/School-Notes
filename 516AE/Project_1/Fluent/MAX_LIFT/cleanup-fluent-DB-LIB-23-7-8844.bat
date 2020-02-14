echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\v192\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\v192\fluent\ntbin\win64\tell.exe" DB-LIB-23-7 62422 CLEANUP_EXITING
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 10872) 
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 9144) 
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 6388) 
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 7256) 
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 8844) 
if /i "%LOCALHOST%"=="DB-LIB-23-7" (%KILL_CMD% 8440)
del "P:\AE516\Project_1\Fluent\MAX_LIFT\cleanup-fluent-DB-LIB-23-7-8844.bat"
