#!/bin/bash
# Start NMR Analyzer on Linux or macOS.
cd "$(dirname "$0")" || exit 1
exec python3 -m nmranalyzer "$@"
