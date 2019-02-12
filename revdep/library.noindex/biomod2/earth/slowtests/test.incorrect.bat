@Rem test.incorrect.R: example incorrect model built by earth
@rem Stephen Milborrow May 2015 Berea

@"C:\PROGRA~1\R\R-3.5.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.incorrect.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.incorrect.Rout:
@echo.
@tail test.incorrect.Rout
@echo test.incorrect.R
@exit /B 1
:good1
@echo mks.diff -w test.incorrect.Rout test.incorrect.Rout.save
@rem egreps to deal with times
@C:\Rtools\bin\echo -n "new "
@egrep "^\[total time" test.incorrect.Rout
@C:\Rtools\bin\echo -n "old "
@egrep "^\[total time" test.incorrect.Rout.save
@egrep -v "^\[total time" test.incorrect.Rout      >test.incorrect.Rout1
@egrep -v "^\[total time" test.incorrect.Rout.save >test.incorrect.Rout.save1
@rem -w to treat \n same as \r\n
@mks.diff -w test.incorrect.Rout1 test.incorrect.Rout.save1
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.incorrect.save.ps
@exit /B 1
:good2
@rem test.incorrect.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.incorrect.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.incorrect.Rout test.incorrect.Rout1 test.incorrect.Rout.save1
@exit /B  0
