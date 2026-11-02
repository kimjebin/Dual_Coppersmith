"""Generate release tables without requiring SageMath or large data unpacking."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from core.reporting import regenerate

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--output",type=Path);args=parser.parse_args()
    status=regenerate(args.output)
    if args.output is None:(ROOT/"results/status.json").write_text(json.dumps(status,indent=2)+"\n")
    print(json.dumps(status,sort_keys=True))
