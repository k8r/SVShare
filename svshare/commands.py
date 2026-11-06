# Run outside tools, and record what was run.
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path


# A run of an outside tool, and the sample it was run on, if any.
@dataclass
class ToolRun:
    sample: str | None
    tool: str
    command: list[str]
    # The file the tool's output was written to, if any.
    stdout: Path | None = None

    def __str__(self):
        text = shlex.join(self.command)
        return f"{text} > {shlex.quote(str(self.stdout))}" if self.stdout else text


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
