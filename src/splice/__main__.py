"""
SPLICE scans the ends of contigs from bacterial genome assemblies looking
for prophage-like region.  Regions identified as prophage fragment candidates
are conjoined based on known genome organization data and extracted with DEPhT
"""
import argparse
import csv
from datetime import datetime
import pathlib
import shutil
import sys
from tempfile import mkdtemp

from Bio import SeqIO
from depht.functions.blastn import BLASTN_EVALUE, BLASTN_OUTFMT
from depht.functions.fasta import write_fasta
from depht.functions.sniff_format import sniff_format
from depht.functions.subprocess import run_command

OUTFMT = "10 qstart qend sseqid sstart send sstrand qcovs"


def parse_args(unparsed_args):
    """
    Parse command line arguments.

    :return: parsed_args:
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("infile", type=pathlib.Path, nargs="+",
                        help="path to genome file(s) to scan for prophages")
    parser.add_argument("outdir", type=pathlib.Path,
                        help="path where output files can be written")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="toggle pipeline verbosity")
    parser.add_argument("-db", "--database", type=pathlib.Path, required=True,
                        help="database of prophage-like reference genomes")
    parser.add_argument("-s", "--search_space", default=40000, type=int,
                        help=("the size of the search space applied over "
                              "ends of contigs"))
    parser.add_argument("-D", "--run_depht", action="store_true")

    return parser.parse_args(unparsed_args)


def main():
    """
    Main function that interfaces with command line args and the program
    workflow
    """
    args = parse_args(sys.argv[1:])

    # Mark program start time
    mark = datetime.now()

    outdir = pathlib.Path(args.outdir).resolve()
    if not outdir.is_dir():
        print(f"'{str(outdir)}' does not exist - creating it...")
        outdir.mkdir(parents=True)

    # TO BE REPLACED
    tmp_dir = pathlib.Path.cwd().joinpath("tmp")
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True)

    phage_db = args.database.resolve()

    splice(args.infile, outdir, phage_db, tmp_dir, verbose=args.verbose,
           run_depht=args.run_depht)

    # shutil.rmtree(tmp_dir)

    print(f"\nTotal runtime: {str(datetime.now() - mark)}")


def splice(infiles, outdir, phage_db, tmp_dir, verbose=False,
           search_space=15000, splice_space=50000, run_depht=False):
    for infile in infiles:
        if not infile.is_file():
            if verbose:
                print("'{str(infile)}' does not exist - skipping it...")
            continue

        if infile.name == ".DS_Store":
            if verbose:
                print("skipping .DS_Store file...")
            continue

        genome_tmp_dir = pathlib.Path(mkdtemp(dir=tmp_dir))
        if not genome_tmp_dir.is_dir():
            genome_tmp_dir.mkdir()

        fmt = sniff_format(infile)

        if not isinstance(fmt, str):
            print(f"Could not recognize file format of '{str(infile)}'.  "
                   "Skipping.")
            continue

        if verbose:
            print(f"\nparsing  '{str(infile)}' as {fmt}...")

        genome_id = infile.stem
        records = [x for x in SeqIO.parse(infile, fmt)
                   if len(x) >= search_space]

        fragment_map = search_contig_ends(
                            records, phage_db, search_space, genome_tmp_dir)

        splice_map = map_splice_events(fragment_map)
        splice_map = clean_splice_map(splice_map)

        if not splice_map:
            if verbose:
                print(f"No fragmented prophages in {genome_id} detected.")
            continue
        else:
            if verbose:
                print(f"{len(splice_map.keys())} fragmented prophages in "
                      f"{genome_id} detected:")

                clusters = list(splice_map.keys())
                clusters.sort()

                for cluster in clusters:
                    print(f"\t- {cluster}-like prophage detected")

        spliced_fragments = splice_fragments(splice_map, records, genome_id)

        genome_outdir = outdir.joinpath(genome_id)
        genome_outdir.mkdir(exist_ok=True)

        mine_spliced_fragments(spliced_fragments, genome_outdir,
                               run_depht=run_depht, verbose=verbose)


def search_contig_ends(records, phage_db, search_space, genome_tmp_dir):
    frag_hits_map = dict()
    for record in records:
        l_seq = record.seq[:search_space]
        r_seq = record.seq[-1*search_space:]

        l_fasta = genome_tmp_dir.joinpath(f"{record.id}_left.fasta")
        write_fasta([record.id], [str(l_seq)], l_fasta)
        r_fasta = genome_tmp_dir.joinpath(f"{record.id}_right.fasta")
        write_fasta([record.id], [str(r_seq)], r_fasta)

        l_results = blastn(l_fasta, phage_db, genome_tmp_dir, outfmt=OUTFMT)
        if l_results:
            load_frag_hits_map(record.id, l_results, "left", frag_hits_map)

        r_results = blastn(r_fasta, phage_db, genome_tmp_dir, outfmt=OUTFMT)
        if r_results:
            load_frag_hits_map(record.id, r_results, "right", frag_hits_map)

    return frag_hits_map


def load_frag_hits_map(contig_id, results, label, frag_hits_map, min_cov=30):
    for result in results:
        if int(result["qcovs"]) < min_cov:
            continue

        contig_hits = frag_hits_map.get(contig_id, dict())

        end_specific_contig_hits = contig_hits.get(label, dict())

        cluster = result["sseqid"].split("_")[0]
        cluster_specific_data = end_specific_contig_hits.get(cluster)

        if cluster_specific_data is None:
            sub_start = result["sstart"]
            sub_end = result["send"]
            if result["sstrand"] == "minus":
                sub_start = result["send"]
                sub_end = result["sstart"]

            end_specific_contig_hits[cluster] = [int(sub_start), int(sub_end),
                                                 result["sstrand"]]
        else:
            sub_start = result["sstart"]
            sub_end = result["send"]
            if result["sstrand"] == "minus":
                sub_start = result["send"]
                sub_end = result["sstart"]

            cluster_specific_data[0] = min((cluster_specific_data[0],
                                            int(sub_start)))
            cluster_specific_data[1] = max((cluster_specific_data[1],
                                            int(sub_end)))

            end_specific_contig_hits[cluster] = cluster_specific_data

        contig_hits[label] = end_specific_contig_hits
        frag_hits_map[contig_id] = contig_hits


def clean_splice_map(splice_map):
    cleaned_splice_map = dict()
    for cluster, fragments in splice_map.items():
        if len(fragments) <= 1:
            continue

        cleaned_fragments = list()
        for fragment in fragments:
            if not cleaned_fragments:
                cleaned_fragments.append(fragment)
                continue

            if cleaned_fragments[-1][0] == fragment[0]:
                if fragment[4] == "minus":
                    if cleaned_fragments[-1][1] == "right" and \
                                fragment[1] == "left":
                        cleaned_fragments[-1][1] = "whole"
                        cleaned_fragments[-1][2] = min((
                                                    cleaned_fragments[-1][2],
                                                    fragment[2]))
                        cleaned_fragments[-1][3] = max((
                                                    cleaned_fragments[-1][3],
                                                    fragment[3]))

                        continue

                if cleaned_fragments[-1][1] == "left" and \
                        fragment[1] == "right":
                    cleaned_fragments[-1][1] = "whole"
                    cleaned_fragments[-1][2] = min((cleaned_fragments[-1][2],
                                                    fragment[2]))
                    cleaned_fragments[-1][3] = max((cleaned_fragments[-1][3],
                                                    fragment[3]))

                    continue

            cleaned_fragments.append(fragment)

        cleaned_splice_map[cluster] = cleaned_fragments

    return cleaned_splice_map


def splice_fragments(splice_map, records, genome_id,
                     splice_space=50000, verbose=True):
    spliced_fragments = list()
    record_seq_map = {record.id: record.seq for record in records}
    for cluster, fragments in splice_map.items():
        joined_contigs = list()
        sequences = list()

        for fragment in fragments:
            contig_seq = record_seq_map[fragment[0]]
            if fragment[4] == "minus":
                contig_seq = contig_seq.reverse_complement()
                if fragment[1] == "left":
                    fragment[1] = "right"
                elif fragment[1] == "right":
                    fragment[1] = "left"

            seq_space = min((splice_space, len(contig_seq)))

            if fragment[1] == "left":
                sequences.append(str(contig_seq[:seq_space]))
            elif fragment[1] == "right":
                sequences.append(str(contig_seq[-1*seq_space:]))
            elif fragment[1] == "whole":
                sequences.append(str(contig_seq))

            joined_contigs.append(str(fragment[0]))

        spliced_name = "_".join([genome_id] + joined_contigs)
        spliced_sequence = "".join(sequences)

        spliced_fragments.append((spliced_name, spliced_sequence, cluster))

    return spliced_fragments


def mine_spliced_fragments(spliced_fragments, outdir, run_depht=False,
                           verbose=False):
    for frag_id, frag_seq, cluster in spliced_fragments:
        fragment_dir = outdir.joinpath(cluster)
        fragment_dir.mkdir(exist_ok=True)

        fragment_fasta = fragment_dir.joinpath(f"{frag_id}.fasta")
        write_fasta([frag_id], [frag_seq], fragment_fasta)

        if run_depht:
            if verbose:
                print(f"\n\tRunning DEPHT on {cluster}-like prophage...")
            depht_outdir = fragment_dir.joinpath("depht_output")

            depht(fragment_fasta, depht_outdir)

            if verbose:
                for prophage_dir in depht_outdir.iterdir():
                    if not prophage_dir.is_dir():
                        continue

                    print("\tDEPHT validated SPLICEd "
                          f"{cluster}-like prophage.")
                    break


def map_splice_events(fragment_map):
    splice_map = dict()

    for contig_id, contig_hits in fragment_map.items():
        for end_label, end_specific_data in contig_hits.items():
            for cluster, cluster_specific_data in end_specific_data.items():
                spliced_segments = splice_map.get(cluster, list())

                fragment_data = [contig_id, end_label] + cluster_specific_data
                spliced_segments.append(fragment_data)

                splice_map[cluster] = spliced_segments

    for spliced_segments in splice_map.values():
        spliced_segments.sort(key=lambda x: x[2])

    return splice_map


def depht(infile, outdir, p=3):
    infile = infile.resolve()
    outdir = outdir.resolve()

    run_command(f"depht '{infile}' '{outdir}' -p {p}")


def blastn(query, target, tmp_dir, mode="db", evalue=BLASTN_EVALUE,
           word_size=None, gapopen=None, gapextend=None, outfmt=BLASTN_OUTFMT):
    """
    Runs blastn in either query/subject mode or query/database mode, as
    indicated by `mode`. Returns hits better than `evalue`.

    NOTE: **kwargs will be interpreted as additional blastn parameters.

    :param query: the query FASTA file to use
    :type query: pathlib.Path
    :param target: the target to BLAST against
    :type target: pathlib.Path
    :param tmp_dir: the directory where temporary files can go
    :type tmp_dir: pathlib.Path
    :param mode: how to treat the target (db or subject sequence)
    :type mode: str
    :param evalue: the e-value cutoff to use
    :type evalue: float
    :param word_size: specify a word size (>=4) to use with blastn
    :type word_size: int or None
    :param gapopen: specify a gap open penalty to use with blastn
    :type gapopen: int or None
    :return: results
    """
    # Store any results here
    results = list()

    # Create output filepath
    outfile = tmp_dir.joinpath(f"{query.stem}_blastn_results.csv")

    # Prepare blastn command
    if mode == "db":
        command = f"blastn -query {query} -db {target}"
    elif mode == "subject":
        command = f"blastn -query {query} -subject {target}"
    else:
        raise ValueError("valid blastn modes are: 'db', 'subject'")
    command += f" -out {outfile} -outfmt '{outfmt}'"

    if evalue:
        command += f" -evalue {evalue}"

    if word_size:
        command += f" -word_size {word_size}"

    if gapopen:
        command += f" -gapopen {gapopen}"

    # TODO:
    #  Gap extend and gap open are dependant on each other, for the future,
    #  there needs to be some logic mandating the use of either none or both
    if gapextend:
        command += f" -gapextend {gapextend}"

    run_command(command)

    # Return parsed hits as list of dictionaries
    fields = outfmt.split()[1:]
    try:
        blastn_reader = open(outfile, "r")
        csv_reader = csv.DictReader(blastn_reader, fieldnames=fields)
        for row in csv_reader:
            results.append(row)
        blastn_reader.close()
    except FileNotFoundError:
        pass

    return results


if __name__ == "__main__":
    main()
