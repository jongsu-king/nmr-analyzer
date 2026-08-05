#!/bin/bash
# Double-click this file in Finder to start NMR Analyzer on macOS.
cd "$(dirname "$0")" || exit 1
exec python3 -m nmranalyzer "$@"
