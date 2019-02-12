@rem test.earthmain.gcc64.bat: test 64 bit standalone earth.c with main()
@rem
@rem TODO I haven't yet been able to get this to work:
@rem      Crashes in daxpy_ call in FindKnot, ok with USE_BLAS = 0.

cp "C:/Program Files/R/R-3.5.0/bin/x64/R.dll" .
                                @if %errorlevel% neq 0 goto error
cp "C:/Program Files/R/R-3.5.0/bin/x64/Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
cp "C:/Program Files/R/R-3.5.0/bin/x64/Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
cp "C:/Program Files/R/R-3.5.0/bin/x64/Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@rem cp "C:/Program Files/R/R-3.5.0/bin/x64/Rzlib.dll" .
@rem                                 @if %errorlevel% neq 0 goto error

@rem you may have to create Rdll_x64.lib and Rblas_x64.lib beforehand
@cp "../../.#/Rdll_x64.lib" Rdll.lib
                                @if %errorlevel% neq 0 goto error
@cp "../../.#/Rblas_x64.lib" Rblas.lib
                                @if %errorlevel% neq 0 goto error

@rem modify the path to include gcc, if needed
@rem only do it if needed
@set | egrep -i "PATH=[^;]*Rtools.mingw_64" >NUL && goto :donesetpath
@echo Modifying path for 64 bit Rtools and R
@set PATH=C:\Rtools\mingw_64\bin;^
C:\Rtools\bin;^
C:\Program Files\R\R-3.2.2\bin\x64;^
C:\Program Files\gs\gs9.19\bin;^
%PATH%
:donesetpath

gcc -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -O3 -std=gnu99^
 -m64^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 Rdll.lib Rblas.lib -o earthmain-gcc64.exe
                                @if %errorlevel% neq 0 goto error
@rem earthmain-gcc64.exe > test.earthmain-gcc64.out
@rem                                 @if %errorlevel% neq 0 goto error
earthmain-gcc64.exe
                                @if %errorlevel% neq 0 goto error

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w test.earthmain-gcc64.out test.earthmain.out64.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll Rdll.lib Rblas.lib earthmain-gcc.* test.earthmain-gcc64.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
