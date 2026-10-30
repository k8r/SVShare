# Run an external tool and raise a readable error if it fails.
import subprocess


# stdout, if given, is an open file the command's output is written to.
def run_command(cmd, stdout=None):
    # With no keyboard input, a tool that asks a question, or runs one that
    # does, gets no answer instead of waiting forever for one.
    result = subprocess.run(cmd, stdin=subprocess.DEVNULL, stdout=stdout,
                            stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"{cmd[0]} failed (exit {result.returncode})\n{result.stderr.strip()}"
        )
