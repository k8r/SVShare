# Run an external tool and raise a readable error if it fails.
import subprocess


# stdout, if given, is an open file the command's output is written to.
def run_command(cmd, stdout=None):
    result = subprocess.run(cmd, stdout=stdout, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} failed (exit {result.returncode})\n{result.stderr.strip()}"
        )
