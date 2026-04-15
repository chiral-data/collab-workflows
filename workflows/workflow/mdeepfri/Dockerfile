FROM python:3.12-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl build-essential && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir mdeepfri

# Patch mDeepFRI pipeline.py to fix IndexError in score dict comprehension
RUN python3 - <<'PYEOF'
import re, pathlib
p = pathlib.Path('/usr/local/lib/python3.12/site-packages/mDeepFRI/pipeline.py')
t = p.read_text()
# Fix IndexError: len(scores) < len(terms) in dict comprehension.
# Handles multi-line form:  terms[i]: float(scores[i])\n    for i in range(len(terms)) ...
t = re.sub(
    r'terms\[i\]:\s*float\(scores\[i\]\)(\s*\n\s*)for i in range\(len\(terms\)\) if float\(scores\[i\]\) >= 0\.1',
    r'term: float(score)\1for term, score in zip(terms, scores) if float(score) >= 0.1',
    t,
)
# Also handle single-line form
t = t.replace(
    'terms[i]: float(scores[i]) for i in range(len(terms)) if float(scores[i]) >= 0.1',
    'term: float(score) for term, score in zip(terms, scores) if float(score) >= 0.1',
)
p.write_text(t)
changed = 'zip(terms, scores)' in t
print(f'pipeline.py patch {"applied" if changed else "NOT applied — pattern not found"}')
assert changed, 'patch failed'
PYEOF

# Download model weights at image build time (v1.1: MF, BP, CC modes)
RUN mDeepFRI get-models -o /opt/mdeepfri-weights -v 1.1

# Pre-build PDB100 MMseqs2 index at build time so runtime skips download
RUN python3 -c "from mDeepFRI.pdb import create_pdb_mmseqs; create_pdb_mmseqs(threads=4)"

WORKDIR /workspace
