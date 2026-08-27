# `report` subcommand: render a formatted report from analysis or filter results.

def run(args):
    output = args.output or (args.results / "report.html")
    print(f"Generating report from {args.results}")
    print(f"Output: {output}")
