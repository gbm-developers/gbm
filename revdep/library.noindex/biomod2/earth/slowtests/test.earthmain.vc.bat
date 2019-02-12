@rem test.earthmain.bat: test the standalone earth.c with main()
@rem
@rem Stephen Milborrow Apr 2007 Petaluma

@set CYGWIN=nodosfilewarning
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

@md Debug

@rem Use -W4 for lint like warnings
cl -nologo -DSTANDALONE -DMAIN -TP -Zi -W3 -I"/a/r/ra/include" -I. -FpDebug\vc60.PCH -Fo"Debug/" -c ..\..\src\earth.c
                                @if %errorlevel% neq 0 goto error
link -nologo -debug -out:earthmain.exe Debug\earth.obj Rdll.lib Rblas.lib
                                @if %errorlevel% neq 0 goto error
earthmain.exe > Debug\test.earthmain.out
                                @if %errorlevel% neq 0 goto error

@rem we use -w on mks.diff so it treats \r\n the same as \n
mks.diff -w Debug\test.earthmain.out test.earthmain.out.save
                                @if %errorlevel% neq 0 goto error

@rm -f R.dll Rblas.dll Rdll.lib Rblas.lib iconv.dll Riconv.dll Rgraphapp.dll Rzlib.dll earthmain.exe *.map *.ilk *.pdb
@rm -rf Debug
@exit /B 0

:error
@exit /B %errorlevel%
