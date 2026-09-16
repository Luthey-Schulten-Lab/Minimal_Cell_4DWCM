#!/usr/bin/env python3
"""
Wall-clock cost of each 4DWCM process, per second of biological time.

Reads the run log (the job's stdout) and sums the per-hook timing lines the
Hook prints, then divides by the biological time the run reached. Give it one
log for a breakdown, or two to compare them.

    analyze_timing.py run.out
    analyze_timing.py baseline.out optimized.out
    analyze_timing.py run.out --until 3600 --plot breakdown.png

Categories match the published breakdown in Thornburg et al., so the numbers
are directly comparable to it.
"""

import argparse
import os
import re
import sys

# Each process, and the line it is timed by. The value is the first number
# after the colon; trailing text is ignored, so lines that carry extra context
# still parse:
#
#   DNA time:  0.418 (next DNA @ 8.0s, interval=4.0s)
#     CME time:      1.2309s
#
# CME is printed indented under a "CME time breakdown:" header, so the pattern
# allows leading whitespace but still requires the colon straight after the
# label -- otherwise the header itself would match and parse as empty.
NUM = r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
PATTERNS = [
    ("Non-hook RDME",  re.compile(r"^\s*Time between hook steps:\s*" + NUM)),
    ("DNA",            re.compile(r"^\s*DNA time:\s*" + NUM)),
    ("Ribosomes",      re.compile(r"^\s*Ribo time(?: \(skipped\))?:\s*" + NUM)),
    ("CME",            re.compile(r"^\s*(?:Total )?CME time:\s*" + NUM)),
    ("ODE",            re.compile(r"^\s*ODE time:\s*" + NUM)),
    ("Comm. + write",  re.compile(r"^\s*Communication and file write time:\s*" + NUM)),
]
BIO_TIME = re.compile(r"^\s*Current biological time:\s*" + NUM)


def parse(path, until=None):
    """Sum each process's times. Stops at `until` biological seconds if given.

    Returns (totals, counts, bio_time).
    """
    totals = {name: 0.0 for name, _ in PATTERNS}
    counts = {name: 0 for name, _ in PATTERNS}
    bio_time = 0.0

    with open(path, errors="replace") as fh:
        for line in fh:
            m = BIO_TIME.match(line)
            if m:
                t = float(m.group(1))
                if until is not None and t > until:
                    break
                bio_time = max(bio_time, t)
                continue
            for name, pat in PATTERNS:
                m = pat.match(line)
                if m:
                    totals[name] += float(m.group(1))
                    counts[name] += 1
                    break

    return totals, counts, bio_time


def report(path, totals, counts, bio_time):
    print("{}  ({:.0f} s biological)".format(os.path.basename(path), bio_time))
    accounted = sum(totals.values())
    print("  {:<16} {:>12} {:>10} {:>9} {:>8}".format(
        "process", "total (s)", "per bio s", "calls", "share"))
    print("  " + "-" * 58)
    for name, _ in PATTERNS:
        per_s = totals[name] / bio_time if bio_time else 0.0
        share = 100 * totals[name] / accounted if accounted else 0.0
        flag = "   <- never matched" if counts[name] == 0 else ""
        print("  {:<16} {:>12.1f} {:>10.3f} {:>9d} {:>7.1f}%{}".format(
            name, totals[name], per_s, counts[name], share, flag))
    print("  " + "-" * 58)
    print("  {:<16} {:>12.1f} {:>10.3f}".format(
        "accounted", accounted, accounted / bio_time if bio_time else 0.0))
    print("  (accounted time excludes startup, replication and teardown, so it"
          " is less than the job's wall clock)")
    print()
    return accounted


def plot(path, labels, series, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    names = [n for n, _ in PATTERNS]
    x = np.arange(len(names))
    width = 0.8 / len(series)

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (label, vals) in enumerate(zip(labels, series)):
        off = (i - (len(series) - 1) / 2) * width
        bars = ax.bar(x + off, [vals[n] for n in names], width, label=label)
        ax.bar_label(bars, fmt="%.2f", padding=2, fontsize=8)

    ax.set_ylabel("wall-clock seconds per biological second")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.grid(axis="y", alpha=0.3)
    if len(series) > 1:
        ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print("plot written to {}".format(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("logs", nargs="+", help="run log(s); two are compared")
    ap.add_argument("--until", type=float, default=None,
                    help="only count up to this biological second (compare "
                         "runs of different lengths on equal footing)")
    ap.add_argument("--plot", default=None, help="write a bar chart here")
    args = ap.parse_args()

    for p in args.logs:
        if not os.path.exists(p):
            sys.exit("no such log: {}".format(p))

    per_second = []
    for p in args.logs:
        totals, counts, bio_time = parse(p, args.until)
        if not bio_time:
            sys.exit("{}: no 'Current biological time:' lines; is this a run log?"
                     .format(p))
        report(p, totals, counts, bio_time)
        per_second.append({n: totals[n] / bio_time for n, _ in PATTERNS})

    if len(per_second) == 2:
        a, b = per_second
        print("{} vs {}, per biological second".format(
            os.path.basename(args.logs[0]), os.path.basename(args.logs[1])))
        print("  " + "-" * 50)
        for name, _ in PATTERNS:
            if a[name]:
                print("  {:<16} {:>8.3f} -> {:>8.3f}  {:+7.1f}%".format(
                    name, a[name], b[name], 100 * (b[name] - a[name]) / a[name]))
            else:
                print("  {:<16} {:>8.3f} -> {:>8.3f}        -".format(
                    name, a[name], b[name]))
        print("  " + "-" * 50)
        print()

    if args.plot:
        plot(args.logs[0], [os.path.basename(p) for p in args.logs],
             per_second, args.plot)


if __name__ == "__main__":
    main()
