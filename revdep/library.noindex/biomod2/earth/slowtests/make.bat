@rem earth/inst/slowtests/make.bat

time /T
@call test.earthmain.gcc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.earthmain.clang.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.earthmain.vc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.earthc.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.mods.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.incorrect.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.big.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.weights.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.glm.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.full.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.pmethod.cv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.varmod.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.varmod.mgcv.bat
                        @if %errorlevel% NEQ 0 goto error
@call test.plotd.bat
                        @if %errorlevel% NEQ 0 goto error
@goto done
:error
@echo ==== ERROR ====
@exit /B %errorlevel%
:done
@rm -f ../../src/earth_res.rc ../Makedeps
@rm -f test.*.pdf *.dll *.lib *.pdb
time /T
@exit /B  0
