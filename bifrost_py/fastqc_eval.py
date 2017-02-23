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


class FastqcSumFile(SeqFileName):
    passvalues = ["PASS", "WARN", "FAIL"]

    properties = ['Basic Statistics',
                  'Per base sequence quality',
                  'Per tile sequence quality',
                  'Per sequence quality scores',
                  'Per base sequence content',
                  'Per sequence GC content',
                  'Per base N content',
                  'Sequence Length Distribution',
                  'Sequence Duplication Levels',
                  'Overrepresented sequences',
                  'Adapter Content',
                  'Kmer Content']

    def __init__(self, name, data):
        super().__init__(name)
        self.passfail = self.parse_fastqc_summary(data)
        self.warn = self.get_val_count("WARN")
        self.fail = self.get_val_count("FAIL")

    def get_val_count(self, val):
        valcount = 0
        for prop in self.passfail:
            if self.passfail[prop] == val:
                valcount += 1
        return valcount

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
    fastqc_sum_file_dict = {}
    for path, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".zip"):
                fullname = os.path.join(path, filename)
                archive_name = filename.rstrip(".zip")
                summary_name = os.path.join(archive_name, "summary.txt")
                with zipfile.ZipFile(fullname) as zf:
                    try:
                        data =zf.read(summary_name)
                    except KeyError:
                        print("ERROR: did not find summary.txt in zip file")
                data = data.decode("utf-8")
                fastqc_sum_file_dict[archive_name] = FastqcSumFile(filename, data)
    return fastqc_sum_file_dict


def sort_files(fastqc_sumfile_dict):
    sorted_files = {}

    for entry in fastqc_sumfile_dict:
        fastqc_sumfile = fastqc_sumfile_dict[entry]
        lane = fastqc_sumfile.lane
        pair = fastqc_sumfile.pairnumber
        if not lane in sorted_files:
            sorted_files[lane] = {}
        if not pair in sorted_files[lane]:
           sorted_files[lane][pair] = []
        sorted_files[lane][pair].append(fastqc_sumfile)

    return sorted_files


def process_file_sets(sorted_files):
    output = ""
    for lane in sorted_files:
        for pair in sorted_files[lane]:
            output += lane + " " + pair + '\n'
            output += "Property\t" + "\t".join(FastqcSumFile.passvalues) + "\n"
            output += create_stats(sorted_files[lane][pair])
            output += "\n"
    return output


def get_bad_files(content_dict):
    failfiles = {}
    for entry in content_dict:
        fastqc_sum_file = content_dict[entry]
        if fastqc_sum_file.fail not in failfiles:
            failfiles[fastqc_sum_file.fail] = []
        failfiles[fastqc_sum_file.fail].append(fastqc_sum_file)

    failwarnfiles = {}
    for entry in content_dict:
        fastqc_sum_file = content_dict[entry]
        failwarn = fastqc_sum_file.warn + fastqc_sum_file.fail
        if failwarn not in failwarnfiles:
            failwarnfiles[failwarn] = []
        failwarnfiles[failwarn].append(fastqc_sum_file)

    mostfails = max(failfiles)
    output = "All files with max count of fails, here: {}\n".format(mostfails)
    output += "\n".join([x.filename for x in failfiles[mostfails]])
    output += "\n\n"

    mostfailwarns = max(failwarnfiles)
    output += "All files with max count of fails and warns, here: {}\n".format(mostfailwarns)
    output += "\n".join([x.filename for x in failwarnfiles[mostfailwarns]]) + "\n"

    return output


def create_stats(file_set):
    properties = FastqcSumFile.properties
    property_count = {key: {} for key in properties}
    for prop in property_count:
        property_count[prop] = {key: 0 for key in FastqcSumFile.passvalues}
    for entry in file_set:
        for prop in properties:
            passvalue = entry.passfail[prop]
            property_count[prop][passvalue] += 1

    output = ""
    for prop in property_count:
        output += prop
        for value in FastqcSumFile.passvalues:
            output += "\t" + str(property_count[prop][value])
        output += "\n"

    return output

def write_output(summary, outputfile):
    fo = open(outputfile, "w")
    fo.write(summary)
    fo.close()


def main(directory, outputfile):
    content_dict = get_summary_filecontents(directory)
    sorted_files = sort_files(content_dict)
    output = process_file_sets(sorted_files)
    output += get_bad_files(content_dict)
    write_output(output, outputfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", metavar="DIRECTORY",
                        help="Directory containing directories with " +
                             "fastq zip files in them")
    parser.add_argument("-o", "--output", metavar="STRING",
                        help="output summary results file name")
    args = parser.parse_args()

    main(args.directory, args.output)