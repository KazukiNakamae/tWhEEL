#!/usr/bin/env python3
"""
Map leaf-level GO IDs to higher-level “parent” terms and
summarise which biological features are fully knock-out-able.

Outputs
  1) counts_table.csv          ← parent_GO, description, n
  2) counts_plot.png           ← 横棒グラフ
       *行数が多ければ treemap でも出力 (--treemap オプション)

Required files
  --go_table   : CSV  (列に GO ID がある)
  --parent_go  : テキスト (1 行 1 GO ID)  …★任意 (無ければ generic GO-Slim)
  --obo        : GO DAG (go-basic.obo, gzip OK) …初回自動 DL
"""
import argparse, csv, gzip, os, sys, textwrap
import pandas as pd
from pathlib import Path
from collections import Counter
from goatools.obo_parser import GODag
import requests, pathlib, matplotlib.pyplot as plt

from go_slim_mapping import map_to_parent, sort_count_rows

GENERIC_SLIM_URL = ("http://current.geneontology.org/ontology/"
                    "subsets/goslim_generic.obo")

def download(url, out):
    print(f"[INFO] Download {url}")
    r = requests.get(url, stream=True, timeout=60)
    r.raise_for_status()
    with open(out, "wb") as w:
        for chunk in r.iter_content(8192):
            w.write(chunk)

def load_go_dag(obo_path: str):
    if not os.path.exists(obo_path):
        download("http://purl.obolibrary.org/obo/go/go-basic.obo", obo_path)
    return GODag(obo_path, optional_attrs={'relationship'})

def load_parent_set(go_dag, parent_file: str | None):
    if parent_file:
        parents = {l.strip() for l in open(parent_file) if l.strip()}
    else:
        slim_path = pathlib.Path("goslim_generic.obo")
        if not slim_path.exists():
            download(GENERIC_SLIM_URL, slim_path)
        slim_dag = GODag(str(slim_path))
        parents = {go for go in slim_dag if go in go_dag}
    return parents

def main():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
        Example:
          python count_go_parent_terms.py \\
            --go_table go_kfi_table_onlyF.csv \\
            --go_col 0 \\
            --parent_go custom_slimgos.txt
        """))
    p.add_argument("--go_table", required=True, help="CSV produced earlier")
    p.add_argument("--go_col", type=int, default=0,
                   help="0-based column index of GO ID (default 0)")
    p.add_argument("--parent_go", help="txt: list of higher-level GO IDs")
    p.add_argument("--obo", default="go-basic.obo",
                   help="GO DAG file (downloaded if missing)")
    p.add_argument("--out_csv", default="parent_counts.csv")
    p.add_argument("--out_png", default="parent_counts.png")
    p.add_argument("--treemap", action="store_true",
                   help="also output treemap_parent_counts.png")
    args = p.parse_args()

    dag = load_go_dag(args.obo)
    parent_set = load_parent_set(dag, args.parent_go)
    print(f"[INFO] Parent term count: {len(parent_set)}")

    df = pd.read_csv(args.go_table)
    go_list = df.iloc[:, args.go_col].tolist()
    mapped = [map_to_parent(go, parent_set, dag) for go in go_list]
    counts = Counter(m for m in mapped if m)

    # --- 出力 CSV ---------------------------------------------------
    rows = sort_count_rows(
        (go, dag[go].name, n) for go, n in counts.items()
    )
    pd.DataFrame(rows,
                 columns=["Parent GO", "Description", "count"]
                ).to_csv(args.out_csv, index=False)
    print(f"[DONE] wrote {args.out_csv}")

    # --- 可視化 (bar) ----------------------------------------------
    labels  = [f"{go} {dag[go].name}" for go, _, _ in rows]
    values  = [n for _, _, n in rows]
    plt.figure(figsize=(10, max(2, 0.3*len(rows))))
    plt.barh(labels, values)
    plt.gca().invert_yaxis()
    plt.xlabel("Number of fully targetable leaf GO terms")
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=300)
    print(f"[DONE] wrote {args.out_png}")

    # --- treemap オプション ----------------------------------------
    if args.treemap and rows:
        try:
            import squarify
            plt.figure(figsize=(10, 6))
            squarify.plot(sizes=values, label=labels, alpha=.8)
            plt.axis('off')
            out_png_path = Path(args.out_png)
            treemap_png_filename = "treemap_" + out_png_path.name
            treemap_png = out_png_path.with_name(treemap_png_filename)
            plt.savefig(treemap_png, dpi=300)
            print(f"[DONE] wrote {treemap_png}")
        except ImportError:
            print("[WARN] squarify missing – install to draw treemap")

if __name__ == "__main__":
    main()
