from .wrapper import Blastn as blast_wrapper
import uuid
import pandas as pd


class Localdb:
    def __init__(self) -> None:
        self.id = None
        self.taxid = 0
        self.dbpath = None
        self.info = None
        self.blast_wrapper = blast_wrapper()

    def get_blastdb_info(self, dbpath):
        text = self.blast_wrapper.get_blast_db_info(dbpath)
        info = self.parse_blast_db_info(text)
        self.info = info
        return info

    def make_from_fasta(self, fasta_path: str, dbname: str = None):
        if dbname is None:
            self.id = str(uuid.uuid4())
        else:
            self.id = dbname
        self.dbpath = self.blast_wrapper.makeblastdb(fasta_path, self.id, self.taxid)
        self.get_blastdb_info(self.dbpath)

    def pull_fasta_from_db(self, 
                           subjects: list, 
                           out_fasta: str, 
                           buffer: bool = False):
        
        subjects = ','.join(subjects)
        buf = self.blast_wrapper.write_db_hits_to_fasta(self.dbpath, subjects, out_fasta, buffer=buffer)
        if buffer:
            return buf

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


class BlastParser:
    def __init__(self) -> None:
        self.wrapper = blast_wrapper()
        self.parameters = {}

    def read_tab_file(self, filepath, outfmt=None):
        if outfmt is None:
            outfmt = self.wrapper.get_default_parameters()['outfmt']
        headers = outfmt.split(' ')[1:]
        out_df = pd.read_csv(filepath, sep='\t')
        out_df.columns = headers

        return out_df


    

    

