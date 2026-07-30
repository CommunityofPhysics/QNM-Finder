@echo off
setlocal enabledelayedexpansion

echo ============================================
echo   Starting 2D sweep over m = 0..10 and k = 0.0..2.0
echo ============================================
echo.

REM Path to finder executable
set EXE="QNM Finder.exe"

REM Outer loop: m = 0, 1, ..., 10
for /L %%M in (0,1,10) do (

    echo ############################################
    echo   Sweeping m = %%M
    echo ############################################
    echo.

    REM Update m in Config.xml
    powershell -Command "(Get-Content 'Config.xml') -replace '<m>.*</m>', '<m>%%M</m>' | Set-Content 'Config.xml'"

    REM Inner loop: k = 0.0, 0.1, ..., 2.0
    for /L %%I in (0,1,20) do (
	REM Integer division and modulo
	set /a "intPart=%%I/10"
    	set /a "fracPart=%%I%%10"

	REM Build k_val = intPart.fracPart
	set "k_val=!intPart!.!fracPart!"


        echo --------------------------------------------
        echo Running finder with m = %%M, k = !k_val!
        echo --------------------------------------------

        REM Update k in Config.xml
        powershell -Command "(Get-Content 'Config.xml') -replace '<k>.*</k>', '<k>!k_val!</k>' | Set-Content 'Config.xml'"

        REM Run finder
        %EXE%

        echo Finished run for m = %%M, k = !k_val!
        echo.
    )
)

echo ============================================
echo   Sweep complete
echo ============================================
pause