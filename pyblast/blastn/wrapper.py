import os
from io import StringIO
from subprocess import call, run

from sh import blastn, makeblastdb, blastdbcmd


def blastn_from_database(self, query: str, out: str, **kwargs) -> StringIO:
    """
    Used for specificity analysis using a blast database.

    :param query: input fasta path
    :param out: .tab output file
    :param kwargs: dictionary of additional optional parameters
    """
    # Default parameters
    params = {
        'query': query,
        'db': 'nt',
        'task': 'blastn',
        'perc_identity': 92,
        'qcov_hsp_perc': 92,
        'max_target_seqs': 5000,
        'threads': 30,
        'buffer': False,
        'word_size': 7,
        'gapopen': 8,
        'gapextend': 6,
        'reward': 5,
        'penalty': -4,
        'max_hsps': 1,
        'evalue': 10000000000000,
        'remote': None,
        'outfmt': '6',
        'no_scientific_names': None,
        'negative_taxid_path': None,
        'dust': 'no'
    }

    # Update default parameters with any user-provided ones
    params.update(kwargs)

    cline = []
    # Building the command line
    for option in params:
        val = params[option]
        if val is not None:
            cline.append('-' + option)
            cline.append(val)

    if params['remote']:
        cline.append("-remote")

    if params['buffer']:
        buf = StringIO()
        blastn(cline, _out=buf)
        buf.seek(0)
        return buf
    else:
        cline.extend(["-out", out])
        blastn(cline)

    def fetch_taxdb(self):
        """
        This file is necessary to run against nt database. A copy must be held locally.
        """
        cline = ["wget", "https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz"]
        call(cline)

        cline = ["tar", "-xvf", "taxdb.tar.gz"]
        call(cline)

    def makeblastdb(
        self,
        input_fasta: str,
        dbname: str,
        taxid: int,
        dbtype: str = "nucl",
        blast_dir="blastdb/",
    ) -> str:
        """create a local blast database

        NOTE Will NOT work if fasta path has a space in it.

        :param input_fasta: input records to create database from
        :type input_fasta: str
        :param dbname: name of blast files
        :type dbname: str
        :param taxid: taxid for fasta
        :type taxid: int
        :param dbtype: default dbtype, defaults to "nucl"
        :type dbtype: str, optional
        :param blast_dir: local blast dir, defaults to "blastdb/"
        :type blast_dir: str, optional
        :return: path to blast paths
        :rtype: str
        """
        if isinstance(blast_dir, str):
            if not os.path.exists(blast_dir):
                os.mkdir(blast_dir)
            if not blast_dir.endswith("/"):
                blast_dir += "/"
            dbname = blast_dir + dbname

        makeblastdb(
            "-dbtype", dbtype, "-in", input_fasta, "-taxid", taxid, "-out", dbname
        )

        return dbname

    def write_db_hits_to_fasta(self, db_path, subjects, out_fasta, buffer=False):
        blastdbcmd("-db", db_path, "-entry", subjects, "-out", out_fasta)

    def get_blast_db_info(self, db_path):
        """
        Extracts information about a BLAST database using blastdbcmd.
        
        :param db_path: Path to the BLAST database (excluding the extension).
        :type db_path: str
        :returns: Information about the BLAST database.
        :rtype: str
        
        This function utilizes the `blastdbcmd` utility from the BLAST+ suite
        to extract database details such as number of sequences, total length,
        database version, and more. The function requires the `blastdbcmd`
        command to be available in the system's PATH, or the full path to the
        executable needs to be specified.
        """
        try:
            # Command to get database information
            result = run(
                ['blastdbcmd', '-db', db_path, '-info'],
                capture_output=True, 
                text=True
            )
            if result.returncode == 0:
                return result.stdout
            else:
                return f"Error: {result.stderr}"
        except Exception as e:
            return str(e)


