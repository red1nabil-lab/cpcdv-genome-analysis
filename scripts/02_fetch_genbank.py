#!/usr/bin/env python
"""Step 2: fetch GenBank records for the CpCDV accessions and extract metadata.

Only accession numbers are sent to NCBI. A generic client contact is used in
the Entrez email field (NCBI requires the field to be populated); no personal
or project information is transmitted.
"""
import os, sys, time, json, csv, re

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
CACHE = os.path.join(OUT, "02_genbank_cache.gb")

from Bio import Entrez, SeqIO

Entrez.email = "plant-virology-client@example.org"   # generic, deliberately not the user's
Entrez.tool = "cpcdv-diversity-analysis"

accs = [l.strip() for l in open(os.path.join(OUT, "01_accessions.txt"), encoding="utf8") if l.strip()]
print(f"{len(accs)} accessions to fetch")

if os.path.exists(CACHE) and os.path.getsize(CACHE) > 10000:
    print("using cached GenBank download:", CACHE)
else:
    chunks = [accs[i:i+50] for i in range(0, len(accs), 50)]
    with open(CACHE, "w", encoding="utf8") as out:
        for n, ch in enumerate(chunks, 1):
            for attempt in range(4):
                try:
                    h = Entrez.efetch(db="nuccore", id=",".join(ch),
                                      rettype="gb", retmode="text")
                    data = h.read()
                    h.close()
                    if "LOCUS" not in data:
                        raise RuntimeError("no LOCUS in reply")
                    out.write(data)
                    print(f"  batch {n}/{len(chunks)}: {len(ch)} ids, {len(data)} bytes")
                    break
                except Exception as e:
                    print(f"  batch {n} attempt {attempt+1} failed: {e}")
                    time.sleep(5 * (attempt + 1))
            else:
                sys.exit(f"batch {n} failed permanently")
            time.sleep(0.5)

# ---------------------------------------------------------------- parse
def qual(feat, key):
    v = feat.qualifiers.get(key)
    return v[0].strip() if v else ""

rows = []
seen = set()
for rec in SeqIO.parse(CACHE, "genbank"):
    acc = rec.id.split(".")[0]
    if acc in seen:
        continue
    seen.add(acc)
    src = next((f for f in rec.features if f.type == "source"), None)
    d = dict(accession=acc,
             version=rec.id,
             definition=rec.description,
             organism=rec.annotations.get("organism", ""),
             length=len(rec.seq),
             gb_date=rec.annotations.get("date", ""),
             taxonomy=" ; ".join(rec.annotations.get("taxonomy", [])))
    if src:
        d.update(country_raw=qual(src, "country") or qual(src, "geo_loc_name"),
                 host=qual(src, "host") or qual(src, "isolation_source"),
                 strain=qual(src, "strain"),
                 isolate=qual(src, "isolate"),
                 clone=qual(src, "clone"),
                 collection_date=qual(src, "collection_date"),
                 collected_by=qual(src, "collected_by"),
                 note=qual(src, "note"),
                 mol_type=qual(src, "mol_type"),
                 lat_lon=qual(src, "lat_lon"))
    else:
        d.update({k: "" for k in ("country_raw host strain isolate clone "
                                  "collection_date collected_by note mol_type "
                                  "lat_lon").split()})
    # references give a publication year when collection_date is absent
    refs = rec.annotations.get("references", [])
    d["journal"] = refs[0].journal if refs else ""
    d["title"] = refs[0].title if refs else ""
    # CDS inventory
    cds = [qual(f, "gene") or qual(f, "product") for f in rec.features if f.type == "CDS"]
    d["cds"] = " | ".join(x for x in cds if x)
    d["n_cds"] = len([f for f in rec.features if f.type == "CDS"])
    rows.append(d)

print(f"\nparsed {len(rows)} GenBank records")
missing = sorted(set(accs) - seen)
if missing:
    print(f"MISSING from GenBank reply ({len(missing)}): {missing}")

flds = list(rows[0].keys())
p = os.path.join(OUT, "02_genbank_metadata.csv")
with open(p, "w", newline="", encoding="utf8") as fh:
    w = csv.DictWriter(fh, fieldnames=flds)
    w.writeheader()
    w.writerows(rows)
print("wrote", p)

from collections import Counter
print("\ncountry_raw values:")
for k, v in Counter(r["country_raw"] for r in rows).most_common():
    print(f"  {k or '(empty)'}: {v}")
print("\nstrain qualifier values:")
for k, v in Counter(r["strain"] for r in rows).most_common():
    print(f"  {k or '(empty)'}: {v}")
print("\nhost values:")
for k, v in Counter(r["host"] for r in rows).most_common():
    print(f"  {k or '(empty)'}: {v}")
print("\norganism values:")
for k, v in Counter(r["organism"] for r in rows).most_common():
    print(f"  {k}: {v}")
