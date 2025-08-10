@echo off
setlocal enabledelayedexpansion

REM --- Configuration ---
set "PSF_FILE=G:\mitochondria_FITC\generated_psf.tif"
set "INPUT_DIR=G:\mitochondria_FITC\Y333 ATP6 ATP3\FITC"
set "OUTPUT_DIR=G:\mitochondria_FITC\Y333 ATP6 ATP3\dw_30"

REM --- Script Logic ---
echo Starting batch processing...
echo Input Directory: %INPUT_DIR%
echo Output Directory: %OUTPUT_DIR%
echo PSF File: %PSF_FILE%

REM Create the output directory if it doesn't exist
if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
    echo Created directory: %OUTPUT_DIR%
)

REM Iterate over all .TIF files in the input directory
for %%f in ("%INPUT_DIR%\*.TIF") do (
    echo -------------------------------------
    echo Processing file: %%~nxf
    echo Full input path: %%f
    
    set "output_file=!OUTPUT_DIR!\dw_%%~nxf"
    echo Full output path: !output_file!
    
    dw --iter 30 --out "!output_file!" "%%f" "%PSF_FILE%"
)

echo -------------------------------------
echo Batch processing complete.

endlocal