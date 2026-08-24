# Author: alfiap (added 2026-05) -- persistent CME solver worker
#
# Long-lived Python process that holds `lm` + GillespieDSolver loaded once.
# Reads CME .lm filenames from stdin (one per line), runs the solver in-process,
# writes "DONE <fname>" or "ERR <fname> <message>" to stdout.
#
# Replaces the per-call `os.system("python Run_CME.py %s")` shell-out, which
# pays ~0.5 s of `lm` import and CUDA-context init on every gCME (1 Hz).
#
# Protocol:
#   stdin:  "<lm_filename>\n"     -> run solver on that file
#           "EXIT\n"               -> shut down cleanly
#           ""                     -> EOF, shut down cleanly
#   stdout: "DONE <lm_filename>\n" on success
#           "ERR <lm_filename>: <message>\n" on failure (worker keeps running)

import sys
import traceback

import lm
from lm import GillespieDSolver


def main():
    print("CME_WORKER: started, lm/GillespieDSolver loaded", flush=True)

    cuda_devices = [0]

    for raw_line in sys.stdin:
        line = raw_line.strip()
        if not line:
            continue
        if line == "EXIT":
            print("CME_WORKER: received EXIT, shutting down", flush=True)
            return
        try:
            lm.runSolver(
                line,
                1,
                solver=GillespieDSolver(),
                cudaDevices=cuda_devices,
                checkpointInterval=0,
            )
            print("DONE " + line, flush=True)
        except Exception as exc:
            tb = traceback.format_exc()
            sys.stderr.write(tb)
            sys.stderr.flush()
            print(f"ERR {line}: {exc}", flush=True)


if __name__ == "__main__":
    main()
