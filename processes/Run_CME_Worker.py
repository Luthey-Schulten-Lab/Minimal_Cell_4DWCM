"""
Persistent CME solver worker process.

Authors
-------
Alfia Parvez — long-lived ``lm`` / GillespieDSolver worker for gCME hooks

Protocol
--------
stdin:  ``<lm_filename>\\n`` run; ``EXIT\\n`` or EOF shut down
stdout: ``DONE <fname>`` or ``ERR <fname>: <message>``

Used by ``MC_CME`` instead of per-call ``os.system(Run_CME.py ...)``.
"""

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
