@rem test.varmod.bat
@rem Stephen Milborrow Dec 2014 Shrewsbury

@"C:\PROGRA~1\R\R-3.5.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.varmod.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.varmod.Rout:
@echo.
@tail test.varmod.Rout
@echo test.varmod.R
@exit /B 1
:good1
mks.diff test.varmod.Rout test.varmod.Rout.save
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.varmod.save.ps
@exit /B 1
:good2
@rem test.varmod.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.varmod.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.varmod.Rout
@rm -f Rplots.ps
@exit /B  0
