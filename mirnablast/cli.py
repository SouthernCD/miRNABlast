import os
import re
import argparse
from toolbiox.lib.common.genome.seq_base import read_fasta_big, sub_seq
from toolbiox.lib.common.os import cmd_run, have_file
from toolbiox.lib.common.fileIO import tsv_file_dict_parse

PATH_OF_THIS_SCRIPT = os.path.split(os.path.realpath(__file__))[0]
BLASTN_PATH = PATH_OF_THIS_SCRIPT + "/dep/blastn"
MAKEBLASTDB_PATH = PATH_OF_THIS_SCRIPT + "/dep/makeblastdb"
RNAFOLD_PATH = PATH_OF_THIS_SCRIPT + "/dep/RNAfold"


def get_lr(brax):
    hash = {}
    lefts = []
    i = 0
    for ch in brax:
        i = i + 1
        if ch == "(":
            lefts.append(i)
        elif ch == ")":
            left = lefts.pop()
            hash[left] = i
    return hash


def get_rl(brax):
    lr = get_lr(brax)
    rl = {}
    for left in lr:
        right = lr[left]
        rl[right] = left
    return rl


def check_mir_star(brax, strand, fold_start, fold_stop, mirkey_start, mirlen):
    if strand == "+":
        offset = mirkey_start - fold_start
    else:
        offset = fold_stop - mirkey_start - mirlen + 1

    if not (offset + mirlen) < len(brax):
        return ("N10",)

    mir_brax = brax[offset:offset + mirlen]
    n_mir_up = mir_brax.count(".")

    # print(n_mir_up)

    if n_mir_up > 5:
        return ("N11",)

    n_bulges = 0
    bulged_nts = 0
    last_right = 0
    last_left = 0

    if re.match(r'^[\.\(]+$', mir_brax):
        # find star of a 5p miRNA
        left_right = get_lr(brax)

        # get one-based start and stop of the mature miRNA, relative to brax
        mir_brax_start = offset + 1
        mir_brax_end = mir_brax_start + mirlen - 1

        # March through the mature miRNA until the last 2 bases (3' overhang),
        # tracking bulges, and getting the star 5p and 3p
        for left in range(mir_brax_start, mir_brax_end - 2):
            if left in left_right:
                if last_right and last_left:  # was it a bulge?
                    left_delta = left - last_left
                    right_delta = last_right - left_right[left]
                    if abs(left_delta - right_delta):
                        n_bulges = n_bulges + 1
                        bulged_nts = bulged_nts + abs(left_delta - right_delta)
                    else:
                        # This is the first pair found in the duplex
                        # Infer the 3p end of the miRNA* .. 2nt offset
                        star_brax_end = left_right[left] + \
                            2 + (left - mir_brax_start)
                last_left = left
                last_right = left_right[left]
        if n_bulges > 2 or bulged_nts > 3:
            return ("N13",)
        else:
            if not last_left:
                return ("N10",)
            # Calculate the star 5p position, relative to brax
            star_brax_start = left_right[last_left] - \
                ((mir_brax_end - 2) - last_left)
            # above, the right-hand term is 0 if the last position analyzed of the mature miRNA was paired.
        # encode and return (below)
    elif re.match(r'^[\.\)]+$', mir_brax):
        # find star of a 3p miRNA
        right_left = get_rl(brax)
        # get one-based start and stop of the mature miRNA, relative to brax
        mir_brax_start = offset + 1
        mir_brax_end = mir_brax_start + mirlen - 1

        # March through the mature miRNA until the last 2 bases (3' overhang),
        # tracking bulges, and getting the star 5p and 3p
        for right in range(mir_brax_start, mir_brax_end - 2):
            if right in right_left:
                if last_right and last_left:
                    # was it a bulge?
                    left_delta = right_left[right] - last_left
                    right_delta = last_right - right
                    if abs(left_delta - right_delta):
                        n_bulges = n_bulges + 1
                        bulged_nts = bulged_nts + abs(left_delta - right_delta)
                else:
                    # This is the first pair found in the duplex
                    # Infer the 3p end of the miRNA* .. 2nt offset
                    star_brax_end = right_left[right] + \
                        (right - mir_brax_start) + 2
                last_left = right_left[right]
                last_right = right
        if n_bulges > 2 or bulged_nts > 3:
            return ("N13",)
        else:
            if not last_right:
                return ("N10",)
            # Calculate the star 5p position, relative to brax
            star_brax_start = right_left[last_right] - \
                ((mir_brax_end - 2) - last_right)
    else:
        return ("N12",)

    if star_brax_start and star_brax_end:
        if star_brax_start >= star_brax_end:
            return ("N10",)
        else:
            if strand == "+":
                star_left = star_brax_start + fold_start - 1
            else:
                star_left = fold_stop - star_brax_end + 1
            mdz = star_brax_end - star_brax_start + 1
            return star_left, strand, mdz
    else:
        return ("N10",)


def blast_and_hit_parse(query_file, db_file, bls_file, aln_ratio, min_len, foldsize, num_threads):
    # load query fasta file
    record_dict = {}
    for record in read_fasta_big(query_file):
        record_dict[record.seqname_short()] = record

    # load db fasta file
    record_db_dict = {}
    for record in read_fasta_big(db_file):
        record_db_dict[record.seqname_short()] = record

    # makeblastdb
    if not have_file(db_file + ".nal"):
        cmd_string = MAKEBLASTDB_PATH + " -in %s -dbtype nucl" % (db_file)
        cmd_run(cmd_string, silence=True)

    # make blast
    cmd_string = BLASTN_PATH + " -query %s -db %s -out %s -word_size 4 -gapopen 5 -gapextend 2 -reward 1 -penalty -3 -evalue 10 -outfmt 6 -num_threads %d" % (
        query_file, db_file, bls_file, num_threads)
    cmd_run(cmd_string, silence=True)

    # load blast results
    outfmt6_title = (
        "query", "subject", "identity", "aln_len", "miss", "gap", "qstart", "qend", "sstart", "send", "evalue", "score")
    hit_dict = tsv_file_dict_parse(bls_file, fieldnames=outfmt6_title)

    # drop_short_hit
    def drop_short_aln(aln_len, query_len):
        aln_len_limit = max(min_len, query_len * aln_ratio)
        if aln_len >= aln_len_limit:
            return 1
        else:
            return 0

    candi_hit_id = [ID for ID in hit_dict if
                    drop_short_aln(abs(int(hit_dict[ID]["sstart"]) - int(hit_dict[ID]["send"])),
                                   len(record_dict[hit_dict[ID]["query"]].seqs))]

    output_list = []
    for ID in candi_hit_id:
        loc_range = (int(hit_dict[ID]["send"]), int(hit_dict[ID]["sstart"]))
        mirkey_start = min(loc_range)
        mirkey_end = max(loc_range)
        loc_delta = int(hit_dict[ID]["send"]) - int(hit_dict[ID]["sstart"])
        mirlen = abs(loc_delta) + 1

        if loc_delta > 0:
            strand = "+"
        else:
            strand = "-"

        loc_center = min(loc_range) + int(0.5 * mirlen)
        fold_start = max(1, loc_center - (int(0.5 * foldsize)))
        fold_stop = min(loc_center + (int(0.5 * foldsize)),
                        len(record_db_dict[hit_dict[ID]["subject"]].seqs))

        fold_seq = sub_seq(
            record_db_dict[hit_dict[ID]["subject"]].seqs, fold_start, fold_stop, strand, True)

        cmd_string = "echo %s | %s --noPS " % (fold_seq, RNAFOLD_PATH)
        flag_tmp, output, error = cmd_run(
            cmd_string, silence=True)
        output = output.decode()
        brax = output.split("\n")[1].split(" ")[0]

        check_mir_star_output = check_mir_star(
            brax, strand, fold_start, fold_stop, mirkey_start, mirlen)

        if len(check_mir_star_output) == 1:
            continue

        star_left, strand, mdz = check_mir_star_output

        query = hit_dict[ID]["query"]
        chr = hit_dict[ID]["subject"]
        mirkey_start = mirkey_start
        mirkey_end = mirkey_end
        strand = strand
        mir_seq = sub_seq(
            record_db_dict[hit_dict[ID]["subject"]].seqs, mirkey_start, mirkey_end, strand, True)
        mistar_start = star_left
        mistar_end = star_left + mdz - 1
        mistar_seq = sub_seq(
            record_db_dict[hit_dict[ID]["subject"]].seqs, mistar_start, mistar_end, strand, True)
        display_start = max(1, min(mistar_start, mistar_end,
                                   mirkey_start, mirkey_end) - 10)
        display_end = min(len(record_db_dict[hit_dict[ID]["subject"]].seqs),
                          max(mistar_start, mistar_end, mirkey_start, mirkey_end) + 10)
        display_fold_seq = sub_seq(record_db_dict[hit_dict[ID]["subject"]].seqs, display_start, display_end, strand,
                                   True)

        cmd_string = "echo %s | %s --noPS " % (display_fold_seq, RNAFOLD_PATH)
        flag_tmp, output, error = cmd_run(
            cmd_string, silence=True)
        output = output.decode()
        display_fold_brax = output.split("\n")[1].split(" ")[0]

        output_list.append((query, chr, strand, mirkey_start, mirkey_end, mir_seq, mistar_start, mistar_end, mistar_seq,
                            display_start, display_end, display_fold_seq, display_fold_brax))
    return output_list


def mirnablast(args):

    bls_file = os.path.dirname(os.path.abspath(
        args.output_file)) + "/" + str(os.getpid()) + ".tmp.bls"
    output = blast_and_hit_parse(args.query_file, args.genome_file, bls_file,
                                 args.aln_ratio, args.min_len, args.foldsize, args.num_threads)

    with open(args.output_file, "w") as f:
        f.write(
            "query\tchr\tstrand\tmirna_start\tmirna_end\tmir_seq\tmistar_start\tmistar_end\tmistar_seq\tdisplay_start\tdisplay_end\tdisplay_fold_seq\tdisplay_fold_brax\n")
        for i in output:
            printer = ""
            for j in i:
                printer = printer + str(j) + "\t"
            printer = printer.rstrip("\t")
            f.write(printer + "\n")


def main():
    parser = argparse.ArgumentParser(
        prog='miRNABlast', description='Use mature query miRNA sequences to BLAST the genome, fold the flanking sequences of the subjects to find possible homologous miRNA locus\n'
    )

    parser.add_argument(
        "query_file", help="Path for query mature miRNA fasta", type=str)
    parser.add_argument("genome_file", help="Path for genome fasta", type=str)
    parser.add_argument("output_file", help="Path for output file", type=str)
    parser.add_argument("-r", "--aln_ratio",
                        help="the blast hit which length small than aln_ratio * query length will not be use. (default as 0.8)",
                        default=0.8, type=float)
    parser.add_argument("-l", "--min_len",
                        help="the blast hit which length small than min_len will not be use. (default as 18)",
                        default=18, type=int)
    parser.add_argument("-s", "--foldsize", help="length of the flanking sequence which used in RNAfold (default as 300)",
                        default=300,
                        type=str)
    parser.add_argument("-p", "--num_threads",
                        help="num of threads for BLASTN", default=1, type=int)

    args = parser.parse_args()

    mirnablast(args)


if __name__ == "__main__":
    main()
