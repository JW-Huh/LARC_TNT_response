"""Provenance: how input/metabolome/fecal_metabolites.csv was rebuilt in November 2026.

The laboratory re-issued its two result workbooks after finding a calculation
error in the first report (v1): cells that should have been reported as not
detected (N.D.) carried a repeated constant instead (all 32 histamine values,
most hyodeoxycholic / ursodeoxycholic / tauroursodeoxycholic+taurohyodeoxycholic
and many taurochenodeoxycholic acid values), a few further cells were N.D., and
taurodeoxycholic acid was re-quantified.  Short-chain fatty acids did not change.

This script copies the corrected values (v2) into the existing table layout:
same 32 rows in the same order, same id columns (AnalysisOrder,
MetaboliteSampleID, SubjectID, Timepoint, SampleID), same 40 metabolite columns.
N.D. is written as 0; nothing else is transformed.  The raw workbooks are not
part of the repository (MetaboLights MTBLS16367 and MTBLS16378 hold the raw data).

Usage:
    python dev/update_fecal_metabolites_v2.py OLD.csv SCFA_RAW_v2.xlsx Tryptophan_Bile_v2.xlsx NEW.csv
Requires openpyxl.
"""
import sys, csv, openpyxl

old_csv, scfa_xlsx, tb_xlsx, out_csv = sys.argv[1:5]

with open(old_csv, newline="", encoding="utf-8") as f:
    rd = csv.reader(f); hdr = next(rd); old_rows = list(rd)
id_cols, met_cols = hdr[:5], hdr[5:]


def read_sheet(path, sheet, first_data_row, name_row):
    """{laboratory sample code: (laboratory row number, {metabolite: value})}"""
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    rows = list(ws.iter_rows(values_only=True))
    names = [str(n).replace(" + ", "_").replace("-", "_").replace(" ", "_") for n in rows[name_row][3:]]
    out = {}
    for r in rows[first_data_row:]:
        if r[1] is None:
            continue
        sid = str(r[1]).strip()
        sid = sid[:-2] if sid.endswith("-S") else sid          # the SCFA sheet uses FD-599-S
        vals = {}
        for n, v in zip(names, r[3:]):
            if isinstance(v, str):
                if v.strip() != "N.D.":
                    raise ValueError(f"unexpected text {v!r} in {sheet} {sid} {n}")
                v = 0.0                                          # not detected
            if v is None:
                raise ValueError(f"empty cell in {sheet} {sid} {n}")
            vals[n] = float(v)
        out[sid] = (int(r[0]), vals)
    return names, out


n1, scfa = read_sheet(scfa_xlsx, "SCFA", 2, 1)
n2, trp = read_sheet(tb_xlsx, "Tryptophan", 3, 2)
n3, ba = read_sheet(tb_xlsx, "Biel acid", 3, 2)                 # sheet name as spelled by the laboratory
assert n1 + n2 + n3 == met_cols, "metabolite columns / order differ from the existing input"

with open(out_csv, "w", newline="", encoding="utf-8") as f:
    f.write(",".join('"%s"' % h for h in hdr) + "\n")
    for r in old_rows:
        ao, sid, subj, tp, samp = r[:5]
        for d in (scfa, trp, ba):
            assert d[sid][0] == int(ao), (sid, ao, d[sid][0])    # laboratory row number = AnalysisOrder
        vals = {**scfa[sid][1], **trp[sid][1], **ba[sid][1]}
        cells = ["0" if vals[c] == 0 else "%.15g" % vals[c] for c in met_cols]
        f.write(",".join([ao, '"%s"' % sid, '"%s"' % subj, '"%s"' % tp, '"%s"' % samp] + cells) + "\n")
print("written", out_csv, "rows", len(old_rows), "cols", len(hdr))
