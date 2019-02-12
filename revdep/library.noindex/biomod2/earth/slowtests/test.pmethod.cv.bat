@Rem test.pmethod.cv.R: example pmethod.cv model built by earth
@rem Stephen Milborrow May 2015 Berea

@"C:\PROGRA~1\R\R-3.5.0\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.pmethod.cv.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.pmethod.cv.Rout:
@echo.
@tail test.pmethod.cv.Rout
@echo test.pmethod.cv.R
@exit /B 1
:good1
@echo mks.diff -w test.pmethod.cv.Rout test.pmethod.cv.Rout.save
@rem egreps to deal with times
@C:\Rtools\bin\echo -n "new "
@egrep "^\[total time" test.pmethod.cv.Rout
@C:\Rtools\bin\echo -n "old "
@egrep "^\[total time" test.pmethod.cv.Rout.save
@egrep -v "^\[total time" test.pmethod.cv.Rout      >test.pmethod.cv.Rout1
@egrep -v "^\[total time" test.pmethod.cv.Rout.save >test.pmethod.cv.Rout.save1
@rem -w to treat \n same as \r\n
@mks.diff -w test.pmethod.cv.Rout1 test.pmethod.cv.Rout.save1
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.pmethod.cv.save.ps
@exit /B 1
:good2
@rem test.pmethod.cv.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.pmethod.cv.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.pmethod.cv.Rout test.pmethod.cv.Rout1 test.pmethod.cv.Rout.save1
@exit /B  0
