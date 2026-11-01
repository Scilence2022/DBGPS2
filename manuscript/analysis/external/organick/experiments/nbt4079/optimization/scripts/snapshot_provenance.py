import hashlib, json
from pathlib import Path
def snapshot_commit(repo):
    return json.loads((Path(repo)/"SOURCE_SNAPSHOT.json").read_text())["source_commit"]
def snapshot_modified(repo):
    repo=Path(repo)
    expected=json.loads((repo/"SOURCE_SNAPSHOT.json").read_text())["source_sha256"]
    return any(hashlib.sha256((repo/name).read_bytes()).hexdigest()!=sha for name,sha in expected.items())
