@rem test.big.bat: This tests earth on a biggish model
@rem This is the test mentioned in the earth man page "Big Models" section
@rem Stephen Milborrow Mar 2008 Durban

@"C:\PROGRA~1\R\R-3.5.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.big.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.big.Rout:
@echo.
@tail test.big.Rout
@echo test.big.R
@exit /B 1
:good1
@echo mks.diff -w test.big.Rout test.big.Rout.save
@rem egreps to deal with times
@C:\Rtools\bin\echo -n "new "
@egrep "^\[total time" test.big.Rout
@C:\Rtools\bin\echo -n "old "
@egrep "^\[total time" test.big.Rout.save
@egrep -v "^\[total time" test.big.Rout      >test.big.Rout1
@egrep -v "^\[total time" test.big.Rout.save >test.big.Rout.save1
@mks.diff test.big.Rout1 test.big.Rout.save1
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.big.save.ps
@exit /B 1
:good2
@rem test.big.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.big.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.big.Rout test.big.Rout1 test.big.Rout.save1
@rm -f Rplots.ps
@exit /B  0
