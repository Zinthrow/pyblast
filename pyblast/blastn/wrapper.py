import os
from io import StringIO
from subprocess import call, run

from sh import blastn, makeblastdb, blastdbcmd

class Blastn(object):
    def blastn_subject(
        self,
        query,
        subject,
        out,
        task="blastn-short",
        perc_identity=92,
        qcov_hsp_perc=92,
        max_target_seqs=10000,
        buffer=False,
        threads=30,
        word_size=7,
    ):
        """
        Input:
        query: str, path to fasta file of universal regions
        subject: str, path to curated fasta file list
        out: str, path to output blast tab file
        perc_identity: int, percentage identical to query entries
        qcov_hsp_perc: int, percentage of query entry covered by hits
        max_target_seqs: int, number of hits to show in tab file

        Used for sensitivity analysis of final universal regions.
        """
        outfmt = "6 qseqid sseqid pident qcovhsp length qlen mismatch gapopen qstart qend sstart send evalue bitscore"
        cline = [
            "-query",
            query,
            "-subject",
            subject,
            "-task",
            task,
            "-word_size",
            word_size,
            "-perc_identity",
            perc_identity,
            "-qcov_hsp_perc",
            qcov_hsp_perc,
            "-max_target_seqs",
            max_target_seqs,
            "-outfmt",
            outfmt,
            "-num_threads",
            threads,
        ]
        cline = [str(c) for c in cline]
        if buffer is True:
            buf = StringIO()
            blastn(cline, _out=buf)
            buf.seek(0)
            return buf
        else:
            cline.extend(["-out", out])
            blastn(cline)

    def blastn_from_database(
        self,
        query: str,
        out: str,
        negative_taxid_path: str,
        db: str = "nt",
        task: str = "blastn",
        perc_identity: int = 92,
        qcov_hsp_perc: int = 92,
        max_target_seqs: int = 5000,
        threads: int = 30,
        buffer: bool = False,
        word_size: int = 7,
        gapopen: int = 8,
        gapextend: int = 6,
        reward: int = 5,
        penalty: int = -4,
        max_hsps: int = 1,
        evalue: int = 10000000000000,
        remote: bool = False,
        outfmt: str = "6",
        no_scientific_names: bool = False,
    ) -> StringIO:  # megablast
        """Used for specificity analysis using a blast database.


        :param query: input fasta path
        :type query: str
        :param out: .tab output file
        :type out: str
        :param negative_taxid_path: path to negative taxids, can be set to none to skip
        :type negative_taxid_path: str
        :param db: blastdb path, defaults to 'nt'
        :type db: str, optional
        :param task: blastn type, defaults to 'blastn-short'
        :type task: str, optional
        :param perc_identity:  defaults to 92
        :type perc_identity: int, optional
        :param qcov_hsp_perc: defaults to 100
        :type qcov_hsp_perc: int, optional
        :param max_target_seqs: defaults to 100
        :type max_target_seqs: int, optional
        :param threads: defaults to 30
        :type threads: int, optional
        :param buffer: defaults to False
        :type buffer: bool, optional
        :param word_size:  defaults to 7
        :type word_size: int, optional
        :param reward: defaults to 1
        :type reward: int, optional
        :param penalty: defaults to -3
        :type penalty: int, optional
        :param gapopen: defaults to 5
        :type gapopen: int, optional
        :param gapextend: defaults to 2
        :type gapextend: int, optional
        :param evalue: defaults to 100000
        :type evalue: int, optional
        :param remote: defaults to False
        :type remote: bool, optional
        :param no_scientific_names: defaults to False
        :type no_scientific_names: bool, optional
        :return: i/o buffer or none
        :rtype: StringIO
        """
        if not db:
            db = "nt"
        if outfmt == "6":
            outfmt = "6 qseqid sseqid pident qcovhsp length qlen mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq staxids sscinames"
        if no_scientific_names is True:
            outfmt = outfmt.replace(" staxids sscinames", "")
        cline = [
            "-db",
            db,
            "-perc_identity",
            perc_identity,
            "-task",
            task,
            "-qcov_hsp_perc",
            qcov_hsp_perc,
            "-dust",
            "no",
            "-max_target_seqs",
            max_target_seqs,
            "-query",
            query,
            "-outfmt",
            outfmt,
            "-word_size",
            word_size,
            "-reward",
            reward,
            "-penalty",
            penalty,
            "-gapopen",
            gapopen,
            "-gapextend",
            gapextend,
            "-max_hsps",
            max_hsps,
            "-num_threads",
            threads,
            "-evalue",
            evalue,
            "-negative_taxidlist",
            negative_taxid_path,
        ]
        cline = [str(c) for c in cline]

        if not negative_taxid_path:
            cline = cline[:-2]

        if remote is True:
            cline.append("-remote")
        if buffer is True:
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


