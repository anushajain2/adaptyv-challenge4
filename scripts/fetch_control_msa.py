"""Fetch a public TfR1 target MSA using the official Boltz/ColabFold API protocol.

Only the public receptor sequence is sent; no designed binder sequences are sent.
"""
import hashlib
import io
import json
import tarfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SERVER = "https://api.colabfold.com"


def request(path, data=None):
    req = urllib.request.Request(
        SERVER+path, data=urllib.parse.urlencode(data).encode() if data else None,
        headers={"User-Agent": "boltz"})
    with urllib.request.urlopen(req, timeout=45) as response:
        return response.read()


def main():
    out = ROOT / "data/controls_msa"
    out.mkdir(exist_ok=True)
    manifest = json.loads((ROOT / "data/controls/manifest.json").read_text())
    seq = next(c["target_sequence"] for c in manifest["cases"] if c["id"] == "3ds18_current")
    record = out / "msa_job.json"
    if record.exists():
        job = json.loads(record.read_text())
        assert job["query_sha256"] == hashlib.sha256(seq.encode()).hexdigest()
    else:
        result = json.loads(request("/ticket/msa", {"q": f">101\n{seq}\n", "mode": "env"}))
        job = {"query_sha256": hashlib.sha256(seq.encode()).hexdigest(),
               "server": SERVER, "mode": "env", **result}
        record.write_text(json.dumps(job, indent=2))
    started = time.monotonic()
    while job.get("status") in {"UNKNOWN", "RUNNING", "PENDING"} and job.get("id"):
        if time.monotonic()-started > 900:
            raise TimeoutError("MSA server wait exceeded 15 minutes; ticket retained for resuming.")
        print("MSA status:", job["status"], flush=True)
        time.sleep(10)
        job.update(json.loads(request("/ticket/"+job["id"])))
        record.write_text(json.dumps(job, indent=2))
    if job.get("status") != "COMPLETE":
        raise RuntimeError(f"MSA server did not complete: {job.get('status')}")
    archive = request("/result/download/"+job["id"])
    (out / "server_result.tar.gz").write_bytes(archive)
    # Read named members without extracting paths from an external tar archive.
    chunks = []
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
        for name in ["uniref.a3m", "bfd.mgnify30.metaeuk30.smag30.a3m"]:
            try:
                member = tar.getmember(name)
            except KeyError:
                continue
            text = tar.extractfile(member).read().decode().replace("\x00", "")
            chunks.append(text.rstrip()+"\n")
    if not chunks:
        raise ValueError("No expected A3M files in server result")
    # Keep the query only once and bound the stored MSA to 512 unique sequences.
    records = []
    for text in chunks:
        header = None
        for line in text.splitlines():
            if line.startswith(">"):
                header = line
            elif line and not line.startswith("#"):
                records.append((header, line))
    assert records[0][1] == seq
    unique, seen = [], set()
    for header, alignment in records:
        aligned = "".join(c for c in alignment if not c.islower())
        if len(aligned) != len(seq):
            raise ValueError("MSA alignment width does not match target")
        if alignment not in seen:
            seen.add(alignment)
            unique.append((header, alignment))
        if len(unique) == 512:
            break
    path = out / "target.a3m"
    path.write_text("".join(f"{h}\n{s}\n" for h, s in unique))
    job.update({"stored_sequences": len(unique), "target_length": len(seq),
                "a3m_sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
    record.write_text(json.dumps(job, indent=2))
    print(json.dumps(job, indent=2))


if __name__ == "__main__":
    main()
