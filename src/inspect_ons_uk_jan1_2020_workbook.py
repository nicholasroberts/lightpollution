#!/usr/bin/env python3
"""Read-only inventory of the ONS *1 January 2020* NUTS XLSX workbook.

DO NOT extract/impute a UK population series from an uninspected workbook.
Inspect source sheets, NUTS3 identifiers and example rows first. The ONS
workbook is NOT committed; it is a 2021 official publication for the
1 January 2020 reference date (not a 30 June mid-year estimate).

Dependency-free XLSX reader: Python stdlib zipfile + ElementTree; no need
to install openpyxl. Does not modify the frozen 1,273-region analyses.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import posixpath
import re
import zipfile
from collections import defaultdict
from pathlib import Path
import xml.etree.ElementTree as ET

from historical_prc_pixel_nuts3 import load_nuts
from project_config import PROJECT_ROOT as ROOT, load_config

FILE = ROOT / "data/raw/eurostat/ons_jan1_2020_nuts_unidemo.xlsx"
OUT = ROOT / "data/processed/publication_figures/uk_population_nuts2016_audit"
ONS_URL = (
    "https://www.ons.gov.uk/file?uri=/peoplepopulationandcommunity/"
    "populationandmigration/populationestimates/adhocs/"
    "13079populationestimatesfor1january2020onnutsboundaries"
    "unifieddemographiccollectionforeurostat/unidemopublishing.xlsx"
)
NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
NS_PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
UK3 = re.compile(r"^UK[A-Z0-9]{3}$")


def rootxml(z: zipfile.ZipFile, path: str):
    return ET.fromstring(z.read(path))


def col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref)
    if not letters:
        return 0
    n=0
    for letter in letters.group():
        n = 26*n+ord(letter)-64
    return n


def strings(z: zipfile.ZipFile):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    obj=rootxml(z, "xl/sharedStrings.xml")
    result=[]
    for si in obj.findall(f"{{{NS_MAIN}}}si"):
        result.append("".join(t.text or "" for t in si.iter(f"{{{NS_MAIN}}}t")))
    return result


def sheets(z: zipfile.ZipFile):
    base=rootxml(z, "xl/workbook.xml")
    relations=rootxml(z, "xl/_rels/workbook.xml.rels")
    targets={
        item.attrib["Id"]: item.attrib["Target"]
        for item in relations.findall(f"{{{NS_PKG}}}Relationship")
    }
    names=[]
    for element in base.findall(f".//{{{NS_MAIN}}}sheet"):
        relid=element.attrib.get(f"{{{NS_REL}}}id")
        if relid not in targets:
            continue
        target=targets[relid]
        if target.startswith("/"):
            member=target.lstrip("/")
        elif target.startswith("xl/"):
            member=target
        else:
            member=posixpath.normpath(posixpath.join("xl", target))
        if member not in z.namelist():
            raise ValueError(f"Missing worksheet XML for {element.attrib['name']}: {member}")
        names.append((element.attrib["name"], member))
    return names


def sheet_rows(z: zipfile.ZipFile, member: str, ss: list[str]):
    """Yield row number and sparse column-position->text values."""
    root=rootxml(z, member)
    for row in root.findall(f".//{{{NS_MAIN}}}sheetData/{{{NS_MAIN}}}row"):
        values={}
        for cell in row.findall(f"{{{NS_MAIN}}}c"):
            ref=cell.attrib.get("r", "")
            col=col_index(ref)
            if col==0:
                continue
            kind=cell.attrib.get("t")
            if kind=="inlineStr":
                value="".join(t.text or "" for t in cell.iter(f"{{{NS_MAIN}}}t"))
            else:
                v=cell.find(f"{{{NS_MAIN}}}v")
                value="" if v is None else (v.text or "")
                if kind=="s" and value:
                    idx=int(value)
                    value=ss[idx] if 0<=idx<len(ss) else f"ERR_SHARED_STRING_{idx}"
            if value:
                values[col]=value.strip()
        if values:
            yield int(row.attrib.get("r", "0")), values


def analyse(file: Path, known_ids: set[str]):
    with zipfile.ZipFile(file) as z:
        ss=strings(z)
        tab_meta=[]
        raw_rows=[]
        matched_codes=set()
        examples=[]
        for name,member in sheets(z):
            counts=set()
            numrows=0
            headers=[]
            for rowno, cols in sheet_rows(z,member,ss):
                numrows+=1
                if len(headers)<4:
                    headers.append({
                        "row":rowno,
                        "cells":{str(k):v for k,v in sorted(cols.items())[:16]},
                    })
                # A row may repeat a code in separate columns, count each code once.
                codes={str(v).strip().upper() for v in cols.values()
                       if UK3.fullmatch(str(v).strip().upper())}
                for code in codes:
                    counts.add(code)
                    matched_codes.add(code)
                    fields={
                        "sheet":name,
                        "sheet_row":rowno,
                        "NUTS_ID":code,
                        "on_locked_NUTS2016":code in known_ids,
                        "row_values_json":json.dumps(
                            {str(k):v for k,v in sorted(cols.items())},
                            ensure_ascii=False
                        ),
                    }
                    raw_rows.append(fields)
                    if code=="UKK11" and len(examples)<10:
                        examples.append(fields)
            tab_meta.append({
                "sheet":name, "xml_path":member,
                "nonempty_rows":numrows,
                "UK_nuts3_codes":len(counts),
                "UK_codes_matching_locked_NUTS2016":len(counts & known_ids),
                "UK_codes_not_in_locked_NUTS2016":sorted(counts-known_ids)[:40],
                "sample_start_rows":headers,
            })
    return tab_meta,raw_rows,matched_codes,examples


def self_test():
    assert col_index("A1")==1
    assert col_index("Z17")==26
    assert col_index("AA123")==27
    assert col_index("BC9")==55
    assert UK3.fullmatch("UKK11") and not UK3.fullmatch("UKK1")
    print("ONS XLSX sheet/cell code inspection tests: PASS")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workbook", type=Path, default=FILE)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--self-test", action="store_true")
    args=p.parse_args()
    self_test()
    if args.self_test:return
    file=args.workbook if args.workbook.is_absolute() else ROOT/args.workbook
    if not file.is_file():
        raise FileNotFoundError(
            f"ONS 1 January 2020 workbook missing: {file}\n"
            f"Official 311-kB XLSX: {ONS_URL}\n"
            "Download that workbook to the path above; do not use ONS mid-year XLS."
        )
    out=args.out if args.out.is_absolute() else ROOT/args.out
    geo=load_nuts(load_config())
    ids=set(geo.loc[geo.CNTR_CODE.eq("UK"),"NUTS_ID"].astype(str))
    if len(ids)!=179:
        raise ValueError(f"Previously audited UK NUTS2016 region count changed: {len(ids)}")
    tab,rows,codes,examples=analyse(file,ids)
    out.mkdir(parents=True,exist_ok=True)
    for name,records,columns in [
        ("ons_2020_uk_nuts3_sheet_code_rows.csv", rows,
         ["sheet","sheet_row","NUTS_ID","on_locked_NUTS2016","row_values_json"]),
    ]:
        with (out/name).open("w",newline="",encoding="utf8") as fh:
            writer=csv.DictWriter(fh,fieldnames=columns)
            writer.writeheader()
            writer.writerows(records)
    summary={
        "official_reference_date":"2020-01-01",
        "official_ONS_release_date":"2021-03-31",
        "original_url":ONS_URL,
        "local_file":str(file),
        "file_sha256":hashlib.sha256(file.read_bytes()).hexdigest(),
        "fixed_uk_nuts2016_nuts3_ids":len(ids),
        "all_workbook_uk_nuts3_codes":len(codes),
        "codes_in_fixed_NUTS2016":len(codes & ids),
        "fixed_codes_not_found":sorted(ids-codes),
        "workbook_codes_not_fixed":sorted(codes-ids),
        "per_sheet":tab,
        "UKK11_example_rows":examples,
        "IMPORTANT":(
            "Code presence alone cannot establish exact 2020 1-Jan population. "
            "Identify correct sex=T age=TOTAL, NUTS3-level and year in original "
            "workbook, check against 2019 Eurostat and aggregate totals. "
            "No NUTS2016 demographic cohort has been created or modified."
        ),
    }
    (out/"ons_2020_uk_nuts3_workbook_inventory.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8"
    )
    print("\nONS 1 JANUARY 2020 WORKBOOK INVENTORY",flush=True)
    print("Sheets:",len(tab)," | archived NUTS2016 UK regions:",len(ids),flush=True)
    for item in tab:
        print(f" {item['sheet']!r}: rows={item['nonempty_rows']}, "
              f"UK NUTS3 codes={item['UK_nuts3_codes']}, "
              f"matching NUTS2016={item['UK_codes_matching_locked_NUTS2016']}",
              flush=True)
        if item["UK_nuts3_codes"]:
            print("  First header rows:",item["sample_start_rows"][:2],flush=True)
    print("\nDistinct region codes matching NUTS2016:",len(codes&ids),"/",len(ids))
    print("Missing codes:",summary["fixed_codes_not_found"][:30],flush=True)
    print("Outside 2016 geography:",summary["workbook_codes_not_fixed"][:30],flush=True)
    print("\nUKK11 (Bristol) cells (up to 3 rows):",flush=True)
    for row in examples[:3]:
        print(row["sheet"],row["sheet_row"],row["row_values_json"][:650],flush=True)
    print("\nAUDIT OUTPUTS",flush=True)
    for f in sorted(out.glob("ons_2020_uk_*")):
        print(f"  {f} ({f.stat().st_size/1024:.1f} KiB)",flush=True)


if __name__=="__main__":
    main()
