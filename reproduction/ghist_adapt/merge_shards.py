"""Merge outputs/<slide>/shard_*/ results into outputs/<slide>/ (runs.csv, per_gene.csv.gz, cores.json)."""
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "outputs"
for slide in sys.argv[1:]:
    main = OUT / slide
    shards = sorted(main.glob("shard_*"))
    runs = [pd.read_csv(main / "runs.csv")] + [pd.read_csv(s / "runs.csv") for s in shards if (s / "runs.csv").exists()]
    genes = [pd.read_csv(main / "per_gene.csv.gz")] + [pd.read_csv(s / "per_gene.csv.gz") for s in shards if (s / "per_gene.csv.gz").exists()]
    cores = json.loads((main / "cores.json").read_text())
    for s in shards:
        if (s / "cores.json").exists():
            c = json.loads((s / "cores.json").read_text())
            cores["cores"].update(c["cores"])
            cores["skipped"] += c["skipped"]
    cores["skipped"] = [dict(t) for t in {tuple(sorted(d.items())) for d in cores["skipped"]}]
    r = pd.concat(runs).drop_duplicates(subset=["core_area_mm2", "n_cores", "placement", "seed", "method"])
    g = pd.concat(genes).drop_duplicates(subset=["config", "method", "gene_idx"])
    r.to_csv(main / "runs.csv", index=False, float_format="%.6g")
    g.to_csv(main / "per_gene.csv.gz", index=False, float_format="%.5g")
    (main / "cores.json").write_text(json.dumps(cores, indent=1))
    for s in shards:
        shutil.rmtree(s)
    print(f"{slide}: {r.groupby(['core_area_mm2','n_cores','placement','seed']).ngroups} configs merged from {len(shards)} shards")
