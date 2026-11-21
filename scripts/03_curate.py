#!/usr/bin/env python
"""Step 3: merge FASTA-header and GenBank metadata into one curated table."""
import csv, re, os, json
from collections import Counter, defaultdict

import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))  # noqa: E402
from cpcdv_paths import DATA_DIR, FIG_DIR, MUSCLE_BIN  # noqa: E402
OUT = DATA_DIR
inv = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "01_inventory_raw.csv"), encoding="utf8"))}
gb  = {r["accession"]: r for r in csv.DictReader(open(os.path.join(OUT, "02_genbank_metadata.csv"), encoding="utf8"))}

# ---------------------------------------------------------------- geography
SUBREGION = {
    "Morocco": "North Africa", "Tunisia": "North Africa", "Egypt": "North Africa",
    "Sudan": "East Africa", "Eritrea": "East Africa", "Kenya": "East Africa",
    "Burkina Faso": "West Africa", "Nigeria": "West Africa",
    "South Africa": "Southern Africa",
    "Iran": "West Asia", "Turkey": "West Asia", "Syria": "West Asia", "Jordan": "West Asia",
    "Yemen": "Arabian Peninsula", "Oman": "Arabian Peninsula",
    "India": "South Asia", "Pakistan": "South Asia",
}
REGION = {  # coarse scheme, matching Table 2 of the review
    "North Africa": "Africa", "East Africa": "Africa", "West Africa": "Africa",
    "Southern Africa": "Africa",
    "West Asia": "West Asia / Arabian Peninsula",
    "Arabian Peninsula": "West Asia / Arabian Peninsula",
    "South Asia": "South Asia",
}
COUNTRY_FIX = {"Turkey": "Türkiye"}

def split_country(raw):
    """'India: A.P., Nandyal' -> ('India', 'A.P., Nandyal')"""
    raw = (raw or "").strip()
    if not raw:
        return "", ""
    if ":" in raw:
        c, loc = raw.split(":", 1)
    else:
        c, loc = raw, ""
    return c.strip(), loc.strip()

# ---------------------------------------------------------------- hosts
HOST = {
    "cicer arietinum": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "chickpea": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "chickpea leaf": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "chickpea cultivar jackey": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "chickpea cultivar pusa bheema (bg5028)": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "chickpea cv. pusa 1053": ("Cicer arietinum", "chickpea", "Fabaceae"),
    "lens culinaris": ("Lens culinaris", "lentil", "Fabaceae"),
    "lentil": ("Lens culinaris", "lentil", "Fabaceae"),
    "vicia faba": ("Vicia faba", "faba bean", "Fabaceae"),
    "faba bean": ("Vicia faba", "faba bean", "Fabaceae"),
    "pisum sativum": ("Pisum sativum", "pea", "Fabaceae"),
    "field pea": ("Pisum sativum", "pea", "Fabaceae"),
    "phaseolous vulgaris": ("Phaseolus vulgaris", "common bean", "Fabaceae"),
    "aeschynomene virginica": ("Aeschynomene virginica", "Virginia jointvetch (weed)", "Fabaceae"),
    "gossypium hirsutum": ("Gossypium hirsutum", "cotton", "Malvaceae"),
    "gossypium arboreum": ("Gossypium arboreum", "tree cotton", "Malvaceae"),
    "okra": ("Abelmoschus esculentus", "okra", "Malvaceae"),
    "okra [abelmoschus esculentus (l.) moench]": ("Abelmoschus esculentus", "okra", "Malvaceae"),
    "wild jute [corchorus olitorius l.]": ("Corchorus olitorius", "wild jute", "Malvaceae"),
    "citrullus lanatus": ("Citrullus lanatus", "watermelon", "Cucurbitaceae"),
    "citrullus lanatus (watermelon)": ("Citrullus lanatus", "watermelon", "Cucurbitaceae"),
    "watermelon": ("Citrullus lanatus", "watermelon", "Cucurbitaceae"),
    "squash": ("Cucurbita pepo", "squash", "Cucurbitaceae"),
    "cucumber": ("Cucumis sativus", "cucumber", "Cucurbitaceae"),
    "solanum lycopersicum": ("Solanum lycopersicum", "tomato", "Solanaceae"),
    "solanum lycopersicum (tomato)": ("Solanum lycopersicum", "tomato", "Solanaceae"),
    "tomato": ("Solanum lycopersicum", "tomato", "Solanaceae"),
    "capsicum annuum": ("Capsicum annuum", "pepper", "Solanaceae"),
    "pepper plant": ("Capsicum annuum", "pepper", "Solanaceae"),
    "pepper fruit": ("Capsicum annuum", "pepper", "Solanaceae"),
    "pepper leaves": ("Capsicum annuum", "pepper", "Solanaceae"),
    "morelle jaune leaves": ("Solanum elaeagnifolium", "silverleaf nightshade (weed)", "Solanaceae"),
    "spinacia oleracea": ("Spinacia oleracea", "spinach", "Amaranthaceae"),
    "brassica nigra": ("Brassica nigra", "black mustard", "Brassicaceae"),
    "carica papaya": ("Carica papaya", "papaya", "Caricaceae"),
    "xanthium strumarium": ("Xanthium strumarium", "rough cocklebur (weed)", "Asteraceae"),
    "heteropogan contartus": ("Heteropogon contortus", "black speargrass (weed)", "Poaceae"),
    "orosius albicinctus": ("Orosius albicinctus", "leafhopper vector", "Cicadellidae (insect)"),
    "amrasca biguttula": ("Amrasca biguttula", "leafhopper (candidate vector)", "Cicadellidae (insect)"),
}

# ---------------------------------------------------------------- strain
def strain_from_isolate(iso):
    """CpCDV_H_SD_SD263_2014 -> H ; 'K ErC347-05' -> K ; 'A IC 148-99' -> A"""
    if not iso:
        return ""
    m = re.match(r"^CpCDV[_-]([A-Z])[_-]", iso)
    if m:
        return m.group(1)
    m = re.match(r"^([A-Z])\s+\S", iso)
    if m:
        return m.group(1)
    return ""

rows = []
for a in sorted(inv):
    i, g = inv[a], gb.get(a, {})
    country, locality = split_country(g.get("country_raw", ""))
    note = []
    if not country and a == "HE610413":
        country, note = "Pakistan", note + ["country inferred from source publication (not in GenBank source qualifier)"]
    sub = SUBREGION.get(country, "")          # look up before renaming
    reg = REGION.get(sub, "")
    country = COUNTRY_FIX.get(country, country)

    hraw = (g.get("host") or "").strip()
    hsp, hcom, hfam = HOST.get(hraw.lower(), ("", "", ""))
    if hraw and not hsp:
        note.append(f"host '{hraw}' unmapped")

    # year
    cd = g.get("collection_date", "")
    m = re.search(r"(19|20)\d\d", cd)
    year = m.group(0) if m else ""
    yr_src = "collection_date" if year else ""
    if not year:
        m = re.search(r"-(\d\d)$", (g.get("isolate") or ""))   # ...ErC347-05
        if m:
            yy = int(m.group(1))
            year = str(1900 + yy) if yy > 50 else str(2000 + yy)
            yr_src = "isolate code (inferred)"
            note.append("year inferred from isolate code")
    if not year:
        m = re.search(r"\((\d{4})\)", g.get("journal", ""))
        if m:
            year = m.group(1)
            yr_src = "publication year (upper bound)"
            note.append("collection year not given; publication year shown")

    # strain: GenBank qualifier > FASTA header > isolate-name pattern
    s_gb, s_hdr = (g.get("strain") or "").strip(), (i.get("strain") or "").strip()
    s_iso = strain_from_isolate(g.get("isolate") or i.get("isolate") or "")
    srcs = []
    strain = ""
    for v, src in ((s_gb, "GenBank strain qualifier"), (s_hdr, "FASTA header"), (s_iso, "isolate designation")):
        if v and re.fullmatch(r"[A-Z]", v):
            if not strain:
                strain, = (v,)
            srcs.append(src)
    if not strain:
        for v, src in ((s_gb, "GenBank strain qualifier"), (s_hdr, "FASTA header")):
            if v:
                strain = v
                srcs.append(src)
                break
    rows.append(dict(
        accession=a,
        virus_name="chickpea chlorotic dwarf virus",
        abbreviation="CpCDV",
        species="Mastrevirus cicerparvi",
        gb_organism=g.get("organism", ""),
        strain=strain,
        strain_evidence="; ".join(srcs),
        isolate=(g.get("isolate") or i.get("isolate") or g.get("clone") or i.get("clone") or ""),
        country=country,
        locality=locality,
        subregion=sub,
        region=reg,
        host_species=hsp,
        host_common=hcom,
        host_family=hfam,
        host_raw=hraw,
        collection_year=year,
        year_source=yr_src,
        length=i["length"],
        molecule_class=i["class_"],
        genome_completeness=i["completeness"],
        genes_covered=i["genes"],
        ambiguous_bases=i["ambiguous"],
        reference_title=g.get("title", ""),
        reference_journal=g.get("journal", ""),
        notes="; ".join(note),
    ))

p = os.path.join(OUT, "03_curated_metadata.csv")
with open(p, "w", newline="", encoding="utf-8-sig") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("wrote", p, f"({len(rows)} rows)")

def tab(key, label):
    print(f"\n== {label} ==")
    for k, v in Counter(r[key] for r in rows).most_common():
        print(f"  {k or '(unassigned)':<34s} {v}")

tab("country", "country")
tab("subregion", "subregion")
tab("region", "region (review scheme)")
tab("host_family", "host family")
tab("strain", "strain (from annotation only)")
tab("molecule_class", "molecule class")
print("\n== completeness of the curated fields ==")
for f in ("country", "subregion", "host_species", "collection_year", "strain", "isolate"):
    n = sum(1 for r in rows if r[f])
    print(f"  {f:<18s} {n}/{len(rows)}  ({100*n/len(rows):.1f}%)")
