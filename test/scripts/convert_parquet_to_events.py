import argparse
from pathlib import Path

import polars as pl
import pod5 as p5
import numpy as np

from tqdm import tqdm

def main(args):

    df = pl.read_parquet(args.parquet)

    borders_series = (
        df
        .group_by("read_id", maintain_order=True)
        .agg(pl.col("event_start"))
        .rename({"event_start": "borders"})
    )

    borders_map = {
        rid: borders
        for rid, borders in zip(borders_series["read_id"], borders_series["borders"])
    }

    pod5_dir = Path(args.pod5)
    pod5_files = sorted(pod5_dir.glob("*.pod5"))
    if not pod5_files:
        raise FileNotFoundError(f"No .pod5 files found in directory: {pod5_dir}")

    with open(args.target, 'w') as out_f:
        for pod5_file in tqdm(pod5_files, desc="pod5 files"):
            with p5.Reader(pod5_file) as reader:
                for r in reader.reads(preload="samples"):
                    rid = str(r.read_id)
                    borders = borders_map.get(rid)
                    if borders is None or borders.is_empty():
                        continue

                    sig_len = len(r.signal)
                    starts = np.asarray(borders)

                    # Drop any starts >= sig_len (out-of-bounds for reduceat)
                    starts = starts[starts < sig_len]
                    if len(starts) == 0:
                        continue

                    ends = np.append(starts[1:], sig_len)
                    lengths = ends - starts

                    sums = np.add.reduceat(r.signal.astype(np.float32), starts)
                    means = sums / lengths

                    # Normalize to mean 0, std 1
                    means = (means - means.mean()) / means.std()        

                    out_f.write(rid + '\t' + '\t'.join(f'{m:.4f}' for m in means) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--parquet', type=str, required=True, help='Path to parquet file with predicted borders')
    parser.add_argument('--pod5', type=str, required=True, help='Path to directory of .pod5 files with the corresponding signals')
    parser.add_argument('--target', type=str, required=True, help='Path to target file')

    main(parser.parse_args())