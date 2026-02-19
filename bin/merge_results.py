#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
from io import StringIO
import sys

# ---------- Tool parsers ----------

def read_signalp(path: Path) -> pd.DataFrame:
    with path.open("r", encoding="utf-8") as fh:
        lines = fh.readlines()
    if len(lines) < 2:
        return pd.DataFrame()
    header = lines[1].lstrip("#").strip()
    cols = header.split("\t")
    data_lines = [l.rstrip("\n") for l in lines[2:] if l.strip()]
    if not data_lines:
        return pd.DataFrame(columns=cols + ['ID_merge', 'signalp_ID_full'])
    df = pd.read_csv(StringIO("\n".join(data_lines)),
                     sep="\t",
                     names=cols,
                     dtype=str,
                     header=None)
    df['ID_merge'] = df['ID'].astype(str).str.strip().str.split().str[0]  # first token
    df['signalp_ID_full'] = df['ID']  # full sequence header
    return df

def read_deeplocpro(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    if df.columns[0].startswith("Unnamed") or df.columns[0] == '':
        df = df.drop(columns=df.columns[0])
    df = df.rename(columns={"ACC": "ID"})
    df['ID_merge'] = df['ID'].astype(str).str.strip().str.split().str[0]
    return df

def read_phobius(path: Path) -> pd.DataFrame:
    rows = []

    with path.open() as fh:
        header = fh.readline()  # skip header
        for line in fh:
            if not line.strip():
                continue

            # split from the right: ID | TM | SP | PREDICTION
            parts = line.rstrip().rsplit(maxsplit=3)
            if len(parts) != 4:
                continue

            seq_id, tm, sp, pred = parts
            rows.append({
                "ID": seq_id.strip(),
                "TM": tm,
                "SP": sp,
                "PREDICTION": pred
            })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows, dtype=str)
    df["ID_merge"] = df["ID"].str.strip().str.split().str[0]
    return df

def read_deeploc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    if "Protein_ID" not in df.columns:
        return pd.DataFrame()
    df = df.rename(columns={"Protein_ID": "ID_merge"})
    return df

def read_deeptmhmm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str)
    if "Sequence_ID" not in df.columns:
        return pd.DataFrame()
    df["ID_merge"] = df["Sequence_ID"].astype(str).str.strip().str.split().str[0]
    return df

# ---------- Utility functions ----------

def make_prefix(path: Path) -> str:
    stem = path.stem
    for suffix in ['_results', '_summary', '_out']:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return stem

def prefix_columns(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    df = df.copy()
    new_cols = {}
    for c in df.columns:
        if c in ['ID_merge', 'signalp_ID_full']:
            new_cols[c] = c
        else:
            cleaned = str(c).strip().replace(' ', '_')
            new_cols[c] = f"{prefix}_{cleaned}"
    df.rename(columns=new_cols, inplace=True)
    # Reorder: ID_merge first, signalp_ID_full next if present
    cols = ['ID_merge']
    if 'signalp_ID_full' in df.columns:
        cols.append('signalp_ID_full')
    cols += [c for c in df.columns if c not in ('ID_merge', 'signalp_ID_full')]
    df = df[cols]
    return df

def load_file(path: Path) -> pd.DataFrame:
    prefix = make_prefix(path)
    try:
        if prefix.startswith('signalp'):
            df = read_signalp(path)
        elif prefix.startswith('deeplocpro') or prefix.startswith('deeploc_pro'):
            df = read_deeplocpro(path)
        elif prefix.startswith('deeploc'):  # DeepLoc (not DeepLocPro)
            df = read_deeploc(path)
        elif prefix.startswith('deeptmhmm'):
            df = read_deeptmhmm(path)
        elif prefix.startswith('phobius'):
            df = read_phobius(path)
        else:
            df = pd.read_csv(path, dtype=str)
            if 'Sequence_ID' in df.columns:
                df = df.rename(columns={'Sequence_ID': 'ID'})
            df['ID_merge'] = df['ID'].astype(str).str.strip().str.split().str[0]
    except Exception as e:
        print(f"Error reading {path}: {e}", file=sys.stderr)
        df = pd.DataFrame(columns=['ID_merge'])

    if df.empty or 'ID_merge' not in df.columns:
        return pd.DataFrame()

    df = prefix_columns(df, prefix)
    return df

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description="Merge protein localization tool results.")
    parser.add_argument('files', nargs='+', help='Result files to merge')
    parser.add_argument('-o', '--out', default='final_merged.tsv', help='Output TSV file')
    args = parser.parse_args()

    dfs = []
    for f in args.files:
        p = Path(f)
        if not p.exists():
            print(f"Warning: file not found, skipping {p}", file=sys.stderr)
            continue
        df = load_file(p)
        if df.shape[0] == 0:
            print(f"Note: {p} had no data rows, skipping.", file=sys.stderr)
            continue
        dfs.append(df)
        print(f"Loaded {p}: {df.shape[0]} rows, {df.shape[1]} cols", file=sys.stderr)

    if not dfs:
        print("No valid files to merge. Writing empty outputs.", file=sys.stderr)

        # Write empty but valid outputs
        pd.DataFrame(columns=["ID_merge"]).to_csv(args.out, sep="\t", index=False)
        pd.DataFrame(columns=["ID_merge"]).to_csv("final_concise.tsv", sep="\t", index=False)

        sys.exit(0)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on='ID_merge', how='outer')

    # After merging all dfs
    merged = merged.sort_values('ID_merge', na_position='last')

    # Keep only ID_merge and signalp_ID_full; drop all other redundant ID columns
    redundant_ids = [c for c in merged.columns if c.endswith('_ID') and c not in ['signalp_ID_full']]
    merged = merged.drop(columns=redundant_ids)

    merged = merged.drop(columns=['signalp_ID_full'], errors='ignore')

    merged.to_csv(args.out, sep='\t', index=False)
    print(f"Wrote merged results to {args.out}")

    # ----- concise summary -----
    concise_cols = [
        "ID_merge",
        "deeplocpro_Localization",
        "signalp_Prediction",
        "signalp_CS_Position",
        "deeptmhmm_Num_TM_helices",
        "deeptmhmm_TM_helices(start-end)",
        "deeptmhmm_Prediction",
        "deeptmhmm_Topology",
        "deeploc_Localizations",
        "deeploc_Signals",
        "deeploc_Membrane_types",
        "phobius_TM",
        "phobius_SP",
        "phobius_PREDICTION",
    ]

    keep = [c for c in concise_cols if c in merged.columns]
    concise = merged[keep]
    concise.to_csv("final_concise.tsv", sep="\t", index=False)
    print("Wrote concise summary to final_concise.tsv")

if __name__ == "__main__":
    main()
