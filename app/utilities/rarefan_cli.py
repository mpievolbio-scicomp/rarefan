#!/usr/bin/env python3

import argparse
import os, sys
import subprocess
from io import StringIO
import shlex

if "REPIN_ecology.jar" in os.listdir():
    JAR = 'REPIN_ecology.jar'
else:
    JAR = os.path.join(os.environ["CONDA_PREFIX"], "lib", 'REPIN_ecology.jar')


MCL_THREADS = max(os.cpu_count()//2, 1)

def rarefan_command(**kwargs):
    cmd = " ".join(['java',
                            '-Xmx10g',
                            '-jar',
                            JAR,
                            kwargs['tmpdir'],
                            kwargs['outdir'],
                            kwargs['reference_strain'],
                            '{}'.format(kwargs['min_nmer_occurrence']),
                            '{}'.format(kwargs['nmer_length']),
                            kwargs['query_rayt_fname'],
                            kwargs['treefile'],
                            '{}'.format(kwargs['e_value_cutoff']),
                            {"y": "true", True: 'true', False: 'false', None: "false"}[kwargs['analyse_repins']],
                            '{}'.format(kwargs.get('num_threads', MCL_THREADS)),
                            '{}'.format(kwargs.get('distance_group_seeds', 15)),
                            ]
                           )

    return cmd


def tree_command(dir, outdir, treefile):

    # Check if treefile exists.
    in_treefile = os.path.join(dir, treefile)
    out_treefile = os.path.join(outdir, treefile)
    if os.path.isfile(in_treefile):
        return "cp {} {}".format(in_treefile, out_treefile)

    # List of genome sequence files in input directory.
    inputs = [os.path.join(dir, f) for f in os.listdir(dir) if f.split(".")[-1] in ["fas", "fna", "fn", "fasta", "fastn"]]

    # Only proceed if >= 4 genomes provided.
    if len(inputs) < 4:
        return "echo 'Tree generation skipped since < 4 genomes in the input data.'"

    # else (no treefile exists.)
    command = "andi -j {} | clustDist > {}".format(" ".join(inputs), out_treefile)

    return command


if __name__ == '__main__':

    # Local import needed here.
    from checkers import parse_results

    greeting = """
***************************************************************************
*                                                                         *
*                          WELCOME TO RAREFAN                             *
*                                                                         *
*                        RAyt And REpin ANalyzer                          *
*                                                                         *
* RAREFAN is released under the terms of the MIT License.                 *
* See LICENSE for details.                                                *
*                                                                         *
* Copyright (c) 2020 - 2022 Max Planck Institute for Evolutionary Biology *
*                                                                         *
***************************************************************************

"""
    print(greeting)

    parser = argparse.ArgumentParser()
    parser.add_argument("indir",
                        metavar="DIR",
                        help="Contains the genome DNA sequences and RAYT AA sequences to be analyzed.",
                        type=str,
    )

    parser.add_argument("-o", "--outdir",
                        dest="outdir",
                        metavar="OUTDIR",
                        default="./rarefan_out",
                        type=str,
                        help='Results will be written to OUTDIR. OUTDIR will be created if not existing (default: %(default)s).',
                        required=False,
    )

    parser.add_argument("-r", "--reference",
                        dest="reference",
                        metavar="REFERENCE",
                        help="Filename of the reference genome sequence",
                        required=True,
                        type=str,
    )

    parser.add_argument("-c", "--min_nmer_occurrence",
                        dest="min_nmer_occurrence",
                        metavar="MIN_NMER_OCCURRENCE",
                        help="Only Nmers of NMER_LENGTH that occur more frequently than MIN_NMER_OCCURRENCE will be taken into account (default: %(default)d). See RAREFAN manual for details.",
                        required=False,
                        default=55,
                        type=int,
    )

    parser.add_argument("-l", "--min_nmer_length",
                        dest="nmer_length",
                        metavar="NMER_LENGTH",
                        help="Only Nmers of NMER_LENGTH that occur more frequently than MIN_NMER_OCCURRENCE will be taken into account (default: %(default)d). See RAREFAN manual for details.",
                        required=False,
                        default=21,
                        type=int,
    )

    parser.add_argument("-d", "--distance_group_seeds",
                        dest="distance_group_seeds",
                        metavar="DISTANCE_GROUP_SEEDS",
                        help="Set the distance between the REPIN groups (???)",
                        required=False,
                        default=15,
                        type=int,
    )


    parser.add_argument("-q", "--query_rayt",
                        dest="query_rayt",
                        metavar="QUERY_RAYT",
                        help="Filename or path of the amino acid sequence file containing the RAYT protein sequence (default: %(default)s).",
                        required=True,
                        type=str,
    )

    parser.add_argument("-e", "--e_value_cutoff",
                        dest="e_value_cutoff",
                        metavar="E_VALUE_CUTOFF",
                        help="e-value cutoff for tblastn of the query rayt sequence against the submitted genomes (default: %(default)s).",
                        required=False,
                        default=1.e-30,
    )

    parser.add_argument("-R", "--no-repins",
                        dest="no_repins",
                        help="Do not analyse REPINS (default: %(default)r).",
                        required=False,
                        action='store_true',
                        default=False,
    )

    parser.add_argument("-j", "--num_threads",
                        dest="num_threads",
                        metavar="THREADS",
                        help="Number of threads for parallel cluster analysis with MCL (default: %(default)d).",
                        required=False,
                        default=MCL_THREADS,
                        type=str,
    )

    parser.add_argument("-t", "--treefile",
                        dest='treefile',
                        metavar="TREEFILE",
                        help="Filename or path of the phylogenetic tree of submitted genomes (newik format, '.nwk' extension). If none given and more than four genomes are submitted, the tree will be calculated and written to OUTDIR/tmptree.nwk (default: %(default)s).",
                        required=False,
                        type=str,
                        default='tmptree.nwk',
    )

    parser.add_argument("-i", "--interactive",
                        dest='interactive',
                        action='store_true',
                        required=False,
                        default=False,
                        help="Interactive mode. Ask for confirmation before starting the analysis run."
    )

    # Parse arguments.
    args = parser.parse_args()

    # Convert to absolute paths.
    args.outdir = os.path.abspath(args.outdir)
    args.query_rayt = os.path.abspath(args.query_rayt)

    # Ask for confirmation if interactive mode.
    if args.interactive:
        items_str = ["", "RAREFAN run parameters:", "=======================",]
        for k,v in args.__dict__.items():
            items_str.append(f"{k}: {v}")

        items_str.append("=======================", )
        items_str.append("\n")
        items_str.append("Please confirm with 'y' or abort with 'n'.  ")

        report = "\n".join(items_str)

        response = input(report)

        if response.lower() != 'y':
            print("Bye, bye!")
            sys.exit(1)

    # Assembla java command.
    java_command = rarefan_command(tmpdir=args.indir,
                                outdir=args.outdir,
                                reference_strain=args.reference,
                                min_nmer_occurrence=args.min_nmer_occurrence,
                                nmer_length=args.nmer_length,
                                query_rayt_fname=args.query_rayt,
                                treefile=args.treefile,
                                e_value_cutoff=args.e_value_cutoff,
                                analyse_repins=(not args.no_repins),
                                mcl_threads=args.num_threads,
                                distance_group_seeds=args.distance_group_seeds,
                                )

    # Log the command.
    print("RAYT/REPIN analysis")
    print("===================")
    print(java_command)

    # Start the subprocess.
    with subprocess.Popen(shlex.split(java_command),
                          shell=False,
                          stdout=subprocess.PIPE,
                          # stderr=subprocess.STDOUT,
                          bufsize=1,
                          universal_newlines=True,
                            ) as proc:

        # Print output directly to stdout.
        for line in proc.stdout:
            print(line, end="")

    if proc.returncode != 0:
        sys.exit(proc.returncode)

    ### Collect results.
    results = parse_results(args.outdir, os.path.basename(args.reference))
    counts = results['counts']

    print("")
    print("Results")
    print("-------")
    print(f"{counts['rayts']} RAYTs across all submitted genomes using tblastn with {args.query_rayt} at an e-value threshold of {args.e_value_cutoff}.")
    print(f"{counts['nmers']} {args.nmer_length}bp long sequences in the reference genome that occur more frequently than {args.min_nmer_occurrence} times.")
    print(f"{sum(counts['repins'].values())} REPINs in the reference genome.")
    print("")

    ### Generate tree.
    print("Tree generation")
    print("===============")
    tree_command = tree_command(args.indir, args.outdir, args.treefile)

    with subprocess.Popen(tree_command,
                          shell=True,
                          stdout=subprocess.PIPE,
                          bufsize=1,
                          universal_newlines=True,
                            ) as proc:

        # Print output directly to stdout.
        for line in proc.stdout:
            print(line, end="")

    if proc.returncode != 0:
        sys.exit(proc.returncode)

    print(f"""

Your RAREFAN run is complete. Data was written to {args.outdir} You can now visualize your results with the command

    $> rarefan_plots -d DIR -r QUERY_RAYT_ID -t TREEFILE -o OUTFILE

    """)
