# Run external tools (SV callers, Jasmine) and raise a readable error if one fails.
import subprocess


def run_command(cmd):
    result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} failed (exit {result.returncode})\n{result.stderr.strip()}"
        )
