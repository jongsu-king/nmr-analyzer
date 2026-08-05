@echo off
REM Double-click this file to start NMR Analyzer on Windows.
cd /d "%~dp0"
python -m nmranalyzer %*
if errorlevel 1 pause
