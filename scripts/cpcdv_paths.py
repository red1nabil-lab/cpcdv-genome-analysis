r"""Path configuration for the CpCDV analysis pipeline.

The scripts in this repository were written against a fixed local directory. For
publication that single constant is resolved here instead, so the pipeline can be run
anywhere without editing the analysis code. Nothing else about the scripts was changed.

Set the data directory with the CPCDV_DATA environment variable:

    export CPCDV_DATA=/path/to/data          # Linux / macOS
    $env:CPCDV_DATA = "C:\path\to\data"    # PowerShell

It defaults to the `data/` directory of this repository, which is where the committed
inputs live, so in most cases nothing needs to be set at all.

MUSCLE (v3.8.31) is required by steps 05 and 08 only. Point MUSCLE_BIN at the binary,
or leave it and those two steps will tell you what is missing; every later step reads
the committed alignment instead.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)

DATA_DIR = os.environ.get("CPCDV_DATA") or os.path.join(_REPO, "data")
FIG_DIR = os.environ.get("CPCDV_FIGURES") or os.path.join(_REPO, "figures")
MUSCLE_BIN = os.environ.get("MUSCLE_BIN", "muscle")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
