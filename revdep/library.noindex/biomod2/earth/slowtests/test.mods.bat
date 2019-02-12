@Rem test.mods.R: test earth's ability to build various models
@rem Stephen Milborrow Jan 2014 Berea

@"C:\PROGRA~1\R\R-3.5.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.mods.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.mods.Rout:
@echo.
@tail test.mods.Rout
@echo test.mods.R
@exit /B 1
:good1
@echo mks.diff -w test.mods.Rout test.mods.Rout.save
@rem egreps to deal with times
@C:\Rtools\bin\echo -n "new "
@egrep "^\[total time" test.mods.Rout
@C:\Rtools\bin\echo -n "old "
@egrep "^\[total time" test.mods.Rout.save
@egrep -v "^\[total time" test.mods.Rout      >test.mods.Rout1
@egrep -v "^\[total time" test.mods.Rout.save >test.mods.Rout.save1
@rem -w to treat \n same as \r\n
@mks.diff -w test.mods.Rout1 test.mods.Rout.save1
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.mods.save.ps
@exit /B 1
:good2
@rem test.mods.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.mods.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.mods.Rout test.mods.Rout1 test.mods.Rout.save1
@rem @rm -f test.mods.pdf
@exit /B  0
