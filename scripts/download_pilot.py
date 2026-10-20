"""Download a Modal run into an existing directory to preserve every remote file."""
import argparse
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("run_id")
    p.add_argument("destination",type=Path)
    args=p.parse_args()
    args.destination.mkdir(parents=True,exist_ok=True)
    subprocess.run([str(ROOT/".venv/bin/modal"),"volume","get",
                    "tfr1-pilot-results",args.run_id,str(args.destination)],check=True)

if __name__=="__main__":
    main()
