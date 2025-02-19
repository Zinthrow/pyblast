from .wrapper import Blastn as blast_wrapper
from .utils import Localdb, BlastParser

class PyBlastnAPI:
    def __init__(self):
        self.helper = None
        self.wrapper = blast_wrapper()
        self.local_db = Localdb()
        self.parser = BlastParser()

    def blast_fasta_against_fasta(self, 
                                  input_fasta, 
                                  database_fasta, 
                                  out_file='out.tab',
                                  **blast_parameters):
        local_db = self.local_db

        local_db.make_from_fasta(database_fasta)

        dbpath = local_db.dbpath

        cmd = self.wrapper
        cmd.blastn_from_database(input_fasta, out_file, db=dbpath, **blast_parameters)

        out_df = self.parser.read_tab_file(out_file, **blast_parameters)

        return out_df
    
    def blast_fasta_against_db(self, 
                              input_fasta,
                              dbpath, 
                              out_file='out.tab',
                              **blast_parameters):

        cmd = self.wrapper
        cmd.blastn_from_database(input_fasta, out_file, db=dbpath, **blast_parameters)

        out_df = self.parser.read_tab_file(out_file, **blast_parameters)

        return out_df