"""Download public reference inputs with URLs and SHA256 provenance."""
import hashlib
import json
import urllib.request
from pathlib import Path
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw"
SOURCES = {
    **{f"{p}.cif": f"https://files.rcsb.org/download/{p}.cif"
       for p in ("11QC", "6W3H", "6WRV", "3S9L")},
    **{f"{p}.fasta": f"https://rest.uniprot.org/uniprotkb/{p}.fasta"
       for p in ("P02786", "Q62351")},
}


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, url in SOURCES.items():
        path = RAW / name
        if not path.exists():
            req = urllib.request.Request(url, headers={"User-Agent": "TfR1-research-pilot/0.1"})
            with urllib.request.urlopen(req, timeout=120) as response:
                payload = response.read()
            path.write_bytes(payload)
        manifest.append({
            "file": str(path.relative_to(ROOT)), "url": url,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "recorded_utc": datetime.now(timezone.utc).isoformat(),
        })
        print(name, path.stat().st_size)
    (RAW / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
