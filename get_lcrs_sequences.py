import os
from subprocess import Popen, PIPE


class LCRsSequences:
    def __init__(self, input_file,
                 output_file,
                 seg_parameters=None,
                 path_to_seg="segmasker"):
        self.input = input_file
        self.output_file = output_file
        self.output = None
        self.path_to_seg = path_to_seg
        self.seg_parameters = seg_parameters
        self.proteins = {}

    def prepare_params(self):
        txt = ""
        if self.seg_parameters is None:
            txt = "-locut 1.5 -hicut 1.8 -window 15"
        else:
            for value, key in self.seg_parameters.items():
                txt += "-" + value + " " + str(key) + " "
        return txt

    def run_seg(self):
        input_fasta = open(self.input, "r").read()
        FNULL = open(os.devnull, 'w')
        params = self.path_to_seg + " " + self.prepare_params() +" "
        p = Popen(params.split(), stdout=PIPE, stdin=PIPE, stderr=FNULL)

        stdout = p.communicate(input=input_fasta.encode("ascii"))[0]
        self.output = stdout.decode()
        parsed_output = self.parse_output()
        self.write_output(parsed_output)

    def parse_output(self):
        retval = {}
        protein_list = self.read_input()
        for line in self.output.splitlines():
            if line.startswith(">"):
                seq = protein_list[line]
                retval[line] = []
                last_header = line
            else:
                line_items = line.split("-")
                beg = int(line_items[0].strip())
                end = int(line_items[1].strip())
                region = seq[beg:end]
                if last_header in retval.keys():
                    retval[last_header].append({"region": region, "end": end, "beg": beg})
                else:
                    retval[last_header] = [{"region": region, "end": end, "beg": beg}]
        return retval

    def read_input(self):
        res = {}
        with open(self.input) as f:
            for line in f.readlines():
                if line.startswith(">"):
                    header = line.strip()
                    res[header] = ""
                elif line:
                    res[header] += line.strip()
        return res

    def write_output(self, parsed_output):
        plik = open(self.output_file, "w")
        for header, regions in parsed_output.items():
            if len(regions) > 0:
                for region in regions:
                    plik.write(header + " LCR:begin=" + str(region['beg']) + ", end=" + str(region['end']) + "\n")
                    plik.write(region['region'] + "\n")
        plik.close()

import click


@click.command()
@click.option('--sequence_db', default="./examples/LCRs.fasta", help='Path to fasta database with LCRs.')
@click.option('--output_file', default="./seg_tmp.fasta", help='Output with clusters of LCRs.')
@click.option('--seg_param', default=None,
              help='Parameters of seg in dictionary. For example:{"locut": 1.5, "hicut":1.8, "window": 15}')
def run(sequence_db, output_file, seg_param):
    if seg_param is not None:
        tmp = LCRsSequences(sequence_db, output_file, eval(seg_param))
    else:
        tmp = LCRsSequences(sequence_db, output_file)

    tmp.run_seg()
    print("poszlo")


if __name__ == "__main__":
    run()
