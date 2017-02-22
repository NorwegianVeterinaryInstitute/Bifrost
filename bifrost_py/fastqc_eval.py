"""
This script takes in a directory consisting of one of more fastqc zipped
results. These are unpacked, and the summary.txt files are analyzed to give an
overview of the quality of the fastqc files in question.

Information about the quality is divided into R1 and R2, and also on lanes, if
multiple lanes have been used.

Usage:

    $ python fastqc_eval.py -d directory_name -o summary_file.txt

"""
import os
import os.path
import argparse
import zipfile

from filetypes import SeqFileName


# figure out the zip files, unzip, get summary.txt file
# per txt file,
#     parse file, put contents into dict
# summarize dictionary

class FastqcSummaryFile(SeqFileName):
    def __init__(self, name, data):
        super().__init__(name)
        self.passfail = self.parse_fastqc_summary(data)

    @staticmethod
    def parse_fastqc_summary(data):
        """
        This staticmethod parses the string found in a fastqc summary.txt file
        :param data: string containing the data
        :return: dictionary with contents of file, key being type, i.e. col 2
        """
        passfail = {}
        lines = data.split("\n")
        for line in lines:
            if line == "":
                break
            fields = line.split("\t")
            key = fields[1]
            value = fields[0]
            passfail[key] = value
        return passfail


def get_summary_filecontents(directory):
    fastqc_sumfile_dict = {}
    for path, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".zip"):
                fullname = os.path.join(path, filename)
                archivename = filename.rstrip(".zip")
                summary_name = os.path.join(archivename, "summary.txt")
                with zipfile.ZipFile(fullname) as zf:
                    try:
                        data =zf.read(summary_name)
                    except KeyError:
                        print("ERROR: did not find summary.txt in zip file")
                data = data.decode("utf-8")
                fastqc_sumfile_dict[archivename] = FastqcSummaryFile(filename, data)
    return fastqc_sumfile_dict

def sort_files(fastqc_sumfile_dict):
    pairs = {fastqc_sumfile_dict[key].pairnumber \
             for key in fastqc_sumfile_dict.keys()}



def create_stats(content_dict):
    pass

def write_output(summary, output):
    pass

def main(directory, output):
    content_dict = get_summary_filecontents(directory)
    summary = create_stats(content_dict)
    write_output(summary, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", metavar="DIRECTORY",
                        help="Directory containing directories with " +
                             "fastq zip files in them")
    parser.add_argument("-o", "--output", metavar="STRING",
                        help="output summary results file name")
    args = parser.parse_args()

    main(args.directory, args.output)