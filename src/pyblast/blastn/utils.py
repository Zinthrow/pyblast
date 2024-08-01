from .wrapper import Blastn
import uuid
from io import StringIO


class Localdb:
    def __init__(self) -> None:
        self.id = None
        self.taxid = 0
        self.dbpath = None
        self.info = None
        self.blastwrapper = Blastn()

    def get_blastdb_info(self, dbpath):
        text = self.blastwrapper.get_blast_db_info(dbpath)
        info = self.parse_blast_db_info(text)
        self.info = info
        return info

    def make_from_fasta(self, fasta_path: str, dbname: str = None):
        if dbname is None:
            self.id = uuid.uuid4()
        else:
            self.id = dbname
        self.dbpath = self.blastwrapper.makeblastdb(fasta_path, self.id, self.taxid)
        self.get_blastdb_info(self.dbpath)

    def parse_blast_db_info(self, text: str):
        """
        Parses information about a BLAST database from a formatted string.

        :param text: Multi-line string containing formatted BLAST database information.
        :type text: str
        :returns: A dictionary with parsed data including database name, sequence count,
                    total bases, date, longest sequence, BLASTDB version, and volume paths.
        :rtype: dict
        """
        lines = text.split('\n')
        if 'Error' in text:
            raise SyntaxError(text + '\n Check to make sure you have the valid blastdb path')
        data = {}

        # Parse the database name, sequence count, and total bases
        first_line = lines[0].split('\s+')
        data['Database'] = first_line[0].split(': ')[1]
        seq_info = lines[1].split('; ')
        data['Sequence Count'] = int(seq_info[0].split(None)[0].replace(',', ''))
        data['Total Bases'] = int(seq_info[1].split(None)[0].replace(',', ''))

        # Parse the date and longest sequence
        date_line = lines[3].split('\t')
        data['Date'] = date_line[0].split(': ')[1]
        data['Longest Sequence'] = int(date_line[1].split(' ')[2].replace(',', ''))

        # Parse BLASTDB version
        data['BLASTDB Version'] = int(lines[5].split(': ')[1])

        # Parse volumes
        data['Volumes'] = lines[8].strip()

        return data


class Parser:
    def __init__(self) -> None:
        self.wrapper = None
        self.parameters = {}


    

    

