"""
This module contains classes and info relating to filetypes handled in this
project
"""

class SeqFileName:
    def __init__(self, filename):
        self.filename = filename
        self.pairnumber = self.pairnumber(filename)
        self.lane = self.lane(filename)

    @staticmethod
    def pairnumber(filename):
        val = "NA"
        if "R1" in filename:
            val = "R1"
        elif "R2" in filename:
            val = "R2"
        return val

    @staticmethod
    def lane(filename):
        lane = "NA"
        fields = filename.split("_")
        for field in fields:
            if field.startswith("L00"):
                lane = field
        return lane


