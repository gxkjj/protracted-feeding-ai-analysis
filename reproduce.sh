#!/bin/sh
set -eu

cd "$(dirname "$0")"
python3 -c 'import sys; sys.exit("Python 3.9 or newer is required") if sys.version_info < (3, 9) else None'
python3 -m unittest discover -s tests -v
python3 reproduce.py
