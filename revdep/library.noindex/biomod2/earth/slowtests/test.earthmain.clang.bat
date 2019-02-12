@rem test.earthmain.clang.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Dec 2014 Shrewsbury

@cp "d:/bin/R320dll/i386/R.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rblas.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Riconv.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rgraphapp.dll" .
                                @if %errorlevel% neq 0 goto error
@cp "d:/bin/R320dll/i386/Rzlib.dll" .
                                @if %errorlevel% neq 0 goto error
@rem you may have to create Rdll.lib and Rblas.lib beforehand
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error
@cp "../../.#/Rblas.lib" .
                                @if %errorlevel% neq 0 goto error
@rem get iconv.dll from /a/r/ra/src/gnuwin32/unicode
@cp "../../.#/Rdll.lib" .
                                @if %errorlevel% neq 0 goto error

@rem modify the path to include clang, if needed
@set | egrep -i "^PATH=.*LLVM" >NUL && goto donesetpath
@echo Modifying path for clang
@set path=C:\Program Files (x86)\LLVM\bin;%PATH%
:donesetpath

@rem Flags same as gcc where possible.
@rem -U_MSC_VER is needed because some clang executables define this inherently https://stackoverflow.com/questions/38499462/how-to-tell-clang-to-stop-pretending-to-be-other-compilers
@rem Some of these warning suppressions are necessary because we use R and BLAS routines.
clang -DSTANDALONE -DMAIN -Wall -pedantic -Wextra -Weverything -O3 -std=gnu99^
 -U_MSC_VER^
 -Wno-strict-prototypes -Wno-reserved-id-macro -Wno-cast-qual -Wno-unknown-pragmas^
 -Wno-float-equal -Wno-format-nonliteral -Wno-padded -Wno-sign-conversion -Wno-undef^
 -Wno-shadow -Wno-missing-prototypes -Wno-deprecated-declarations -Wno-implicit-function-declaration^
 -Wno-missing-noreturn^
 -I"/a/r/ra/include" -I../../inst/slowtests ../../src/earth.c^
 Rdll.lib Rblas.lib -o earthmain-clang.exe
                                @if %errorlevel% neq 0 goto error
@earthmain-clang.exe > test.earthmain-clang.out
                                @if %errorlevel% neq 0 goto error

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w test.earthmain-clang.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Riconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll Rdll.lib Rblas.lib earthmain-clang.* test.earthmain-clang.* *.o
@exit /B 0

:error
@exit /B %errorlevel%
