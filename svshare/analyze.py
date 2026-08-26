from .reference import check_bam_reference_compatibility, ensure_reference_index

def run(args):
    ensure_reference_index(args.reference)

    print(f"Analyzing {len(args.samples)} sample(s):")
    for sample in args.samples:
        print(f"  - {sample}")
        for issue in check_bam_reference_compatibility(sample, args.reference):
            print(f"    ! {issue}")
    print(f"Reference: {args.reference}")
    print(f"Output directory: {args.output}")
