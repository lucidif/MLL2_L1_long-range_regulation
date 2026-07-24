#!/usr/bin/env python3
import argparse
import os
import sys

import numpy as np
import cooler


def json_sanitize(x):
    """
    Converte ricorsivamente oggetti numpy (es. np.int64) in tipi Python nativi
    per renderli JSON-serializzabili (cooler salva metadata come JSON).
    """
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, (np.ndarray,)):
        return [json_sanitize(i) for i in x.tolist()]
    if isinstance(x, dict):
        return {str(k): json_sanitize(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [json_sanitize(i) for i in x]
    return x


def lift_one(in_uri: str, out_path: str, target_chrom: str, offset: int, assembly: str | None):
    c = cooler.Cooler(in_uri)

    # --- bins: sposta coordinate e rinomina chrom (preserva colonne extra tipo weight)
    bins = c.bins()[:].copy()
    bins["chrom"] = target_chrom
    bins["start"] = bins["start"].astype("int64") + offset
    bins["end"]   = bins["end"].astype("int64") + offset

    # --- pixels: NON cambiano
    pixels = c.pixels()[:]

    # --- metadata: sanitize per evitare np.int64 non serializzabili
    meta = json_sanitize(dict(getattr(c, "info", {}) or {}))
    if assembly:
        meta["assembly"] = assembly

    # Crea il nuovo cooler
    cooler.create_cooler(
        out_path,
        bins=bins,
        pixels=pixels,
        metadata=meta,
        ordered=True,
        ensure_sorted=True
    )


def main():
    ap = argparse.ArgumentParser(
        description="Lift-back di un cooler: aggiunge offset a bins.start/end e imposta bins.chrom=--chrom. Pixels invariati."
    )
    ap.add_argument(
        "-i", "--input", required=True,
        help="Input cooler URI. Esempi: sample.cool oppure sample.mcool::/resolutions/10000"
    )
    ap.add_argument(
        "-o", "--output", required=True,
        help="Output file .cool (verrà creato)."
    )
    ap.add_argument("--chrom", required=True, help="Cromosoma target (es. chr10)")
    ap.add_argument("--offset", required=True, type=int, help="Offset genomico (es. 51627809)")
    ap.add_argument("--assembly", default=None, help="Metadata assembly (es. mm10). NON è un path.")

    args = ap.parse_args()

    # evita errori per output in directory non esistente
    out_dir = os.path.dirname(os.path.abspath(args.output))
    if out_dir and not os.path.isdir(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # evita overwrite accidentale
    if os.path.exists(args.output):
        sys.stderr.write(f"ERROR: output esiste già: {args.output}\n")
        sys.exit(2)

    lift_one(args.input, args.output, args.chrom, args.offset, args.assembly)


if __name__ == "__main__":
    main()
