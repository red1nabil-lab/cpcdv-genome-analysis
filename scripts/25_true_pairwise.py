#!/usr/bin/env python
"""Step 25: recompute the ENTIRE identity matrix from true pairwise alignments.

Reviewer criticism 1: the ICTV mastrevirus demarcation procedure specifies pairwise
alignment with pairwise deletion of gaps, whereas our matrix was derived from a single
MUSCLE multiple sequence alignment.  Identity estimated from an MSA can differ from
identity estimated pairwise, and the difference matters most near 78% and 94%.

So we do not argue about it: all 35,511 pairs are realigned independently by global
Needleman-Wunsch (affine gaps) and identity recomputed under both gap conventions.  The
MSA-derived and pairwise-derived matrices are then compared, and every threshold
violation is re-counted under the pairwise values.  If the violations are an artefact of
the MSA they will disappear here.

Robustness: the whole matrix is computed under two independent scoring schemes, so the
conclusions cannot rest on one arbitrary gap penalty.
"""
import os, sys, json, time, itertools
import numpy as np
from multiprocessing import Pool

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
AN = DATA_DIR
PARAMS = {
    "A": dict(match=1.0, mismatch=-2.0, open_gap=-10.0, extend_gap=-0.5),
    "B": dict(match=5.0, mismatch=-4.0, open_gap=-12.0, extend_gap=-1.0),  # EDNAFULL-like
}

_G = {}


def load():
    from Bio import SeqIO
    cfg = json.load(open(os.path.join(AN, "10_analysis_set.json")))
    keep = list(cfg["keep"])
    seqs = {r.id.split()[0]: str(r.seq).upper()
            for r in SeqIO.parse(os.path.join(AN, "04_full_oriented.fasta"), "fasta")}
    missing = [k for k in keep if k not in seqs]
    if missing:
        raise SystemExit("missing oriented genomes: " + ",".join(missing))
    return keep, {k: seqs[k] for k in keep}


def init(accs, seqs, pname):
    from Bio import Align
    p = PARAMS[pname]
    al = Align.PairwiseAligner()
    al.mode = "global"
    al.match_score = p["match"]; al.mismatch_score = p["mismatch"]
    al.open_gap_score = p["open_gap"]; al.extend_gap_score = p["extend_gap"]
    _G["al"] = al
    _G["accs"] = accs
    _G["arr"] = {k: np.frombuffer(v.encode(), dtype=np.uint8) for k, v in seqs.items()}


def one(ij):
    i, j = ij
    accs = _G["accs"]; arr = _G["arr"]; al = _G["al"]
    a, b = accs[i], accs[j]
    A, B = arr[a], arr[b]
    aln = al.align(A.tobytes().decode(), B.tobytes().decode())[0]
    ba, bb = aln.aligned
    matches = 0; ungapped = 0
    for (s1, e1), (s2, e2) in zip(ba, bb):
        matches += int(np.count_nonzero(A[s1:e1] == B[s2:e2]))
        ungapped += e1 - s1
    L = aln.length
    return i, j, (100.0 * matches / ungapped if ungapped else 0.0), \
           (100.0 * matches / L if L else 0.0)


def run(pname, accs, seqs):
    n = len(accs)
    pairs = list(itertools.combinations(range(n), 2))
    pw = np.full((n, n), 100.0); gp = np.full((n, n), 100.0)
    t0 = time.time()
    nproc = max(1, (os.cpu_count() or 4) - 1)
    with Pool(nproc, initializer=init, initargs=(accs, seqs, pname)) as pool:
        for c, (i, j, a, b) in enumerate(pool.imap_unordered(one, pairs, chunksize=200), 1):
            pw[i, j] = pw[j, i] = a
            gp[i, j] = gp[j, i] = b
            if c % 5000 == 0:
                el = time.time() - t0
                print(f"  [{pname}] {c}/{len(pairs)} pairs  {el:.0f}s "
                      f"(eta {el/c*(len(pairs)-c):.0f}s)", flush=True)
    print(f"  [{pname}] done in {time.time()-t0:.0f}s on {nproc} processes", flush=True)
    return pw, gp


if __name__ == "__main__":
    accs, seqs = load()
    print(f"genomes: {len(accs)}  pairs: {len(accs)*(len(accs)-1)//2}", flush=True)
    out = {}
    for pname in ("A", "B"):
        pw, gp = run(pname, accs, seqs)
        np.save(os.path.join(AN, f"25_truepw_pid_pw_{pname}.npy"), pw)
        np.save(os.path.join(AN, f"25_truepw_pid_gap_{pname}.npy"), gp)
        out[pname] = (pw, gp)
    json.dump({"accessions": accs, "params": PARAMS},
              open(os.path.join(AN, "25_truepw_meta.json"), "w"), indent=1)
    # quick agreement check between the two scoring schemes
    a = out["A"][0]; b = out["B"][0]
    iu = np.triu_indices(len(accs), 1)
    d = np.abs(a[iu] - b[iu])
    print(f"\nscheme A vs B (pairwise deletion): mean|diff| = {d.mean():.3f} pp, "
          f"max = {d.max():.3f} pp, r = {np.corrcoef(a[iu], b[iu])[0,1]:.5f}")
    print("wrote 25_truepw_* to", AN)
