#!/usr/bin/env python

from contextlib import contextmanager
from psycopg2 import connect, ProgrammingError, OperationalError
from gzip import open as gzopen
from itertools import izip

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "BSD"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

FULL_RECORD_DUMP = "SELECT gg_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,"\
                   "decision,prokmsaname,isolation_source,clone,organism,"\
                   "strain,specific_host,authors,title,journal,pubmed,"\
                   "submit_date,country,nt.tax_string AS ncbi_tax_string,"\
                   "st.tax_string AS silva_tax_string,"\
                   "gg.tax_string AS greengenes_tax_string,"\
                   "h.tax_string AS hugenholtz_tax_id,non_acgt_percent,"\
                   "perc_ident_to_invariant_core,aseq.sequence AS aligned_seq "\
                   "FROM record g "\
                   "LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id "\
                   "LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id "\
                   "LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id "\
                   "LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id "\
                   "LEFT JOIN sequence aseq ON aseq.seq_id=g.%s "\
                   "WHERE g.gg_id IN (%s)"

SINGLE_RECORD = "SELECT gg_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,"\
                "decision,prokmsaname,isolation_source,clone,organism,"\
                "strain,specific_host,authors,title,journal,pubmed,"\
                "submit_date,country,nt.tax_string AS ncbi_tax_string,"\
                "st.tax_string AS silva_tax_string,"\
                "gg.tax_string AS greengenes_tax_string,"\
                "h.tax_string AS hugenholtz_tax_id,non_acgt_percent,"\
                "perc_ident_to_invariant_core,ssu.sequence AS ssualign_seq,"\
                "pn.sequence AS pynast_seq,"\
                "ua.sequence AS unaligned_seq "\
                "FROM record g "\
                "LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id "\
                "LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id "\
                "LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id "\
                "LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id "\
                "LEFT JOIN sequence ssu ON ssu.seq_id=g.aligned_seq_id "\
                "LEFT JOIN sequence pn on pn.seq_id=g.pynast_aligned_seq_id "\
                "LEFT JOIN sequence ua on ua.seq_id=g.unaligned_seq_id "\
                "WHERE g.gg_id='%s' or g.ncbi_acc_w_ver='%s'"

SINGLE_RECORD_ORDER = ['gg_id', 'ncbi_acc_w_ver', 'ncbi_gi', 'db_name',
                       'gold_id', 'decision', 'prokmsaname',
                       'isolation_source', 'clone', 'organism', 'strain',
                       'specific_host', 'authors', 'title', 'journal',
                       'pubmed', 'submit_date', 'country', 'ncbi_tax_string',
                       'silva_tax_string', 'greengenes_tax_string',
                       'hugenholtz_tax_id', 'non_acgt_percent',
                       'perc_ident_to_invariant_core', 'ssualign_seq',
                       'pynast_seq', 'unaligned_seq']

FULL_RECORD_ORDER = ["gg_id", "ncbi_acc_w_ver", "ncbi_gi", "db_name",
                     "gold_id", "decision", "prokmsaname", "isolation_source",
                     "clone", "organism", "strain", "specific_host", "authors",
                     "title", "journal", "pubmed", "submit_date", "country",
                     "ncbi_tax_string", "silva_tax_string",
                     "greengenes_tax_string", "hugenholtz_tax_id",
                     "non_acgt_percent", "perc_ident_to_invariant_core",
                     "aligned_seq"]

FULL_GG_ORDER = ["gg_id", "ncbi_acc_w_ver", "ncbi_gi", "db_name",
                 "gold_id", "decision", "prokmsaname", "isolation_source",
                 "clone", "organism", "strain", "specific_host", "authors",
                 "title", "journal", "pubmed", "submit_date", "country",
                 "ncbi_tax_id", "silva_tax_id",
                 "greengenes_tax_id", "hugenholtz_tax_id",
                 "non_acgt_percent", "perc_ident_to_invariant_core",
                 "unaligned_seq_id", "aligned_seq_id",
                 "pynast_aligned_seq_id"]

_sql_insert_tax = """INSERT INTO taxonomy(tax_id, tax_version, tax_string)
                     VALUES (%d, '%s', '%s')"""
_sql_insert_seq = """INSERT INTO sequence(seq_id, sequence)
                     VALUES (%d, '%s')"""
_sql_update_rec = """UPDATE record
                     SET %s=%s
                     WHERE gg_id=%d"""
_sql_insert_otu_cluster = """
               INSERT INTO otu_cluster(cluster_id, rep_id, rel_id, similarity,
                                       method)
               VALUES (%d, %d, %d, %f, '%s')"""
_sql_insert_otu = """INSERT INTO otu(cluster_id, gg_id)
                     VALUES (%d, %d)"""
_sql_create_tmp = "CREATE TEMPORARY TABLE %s LIKE %s"""
_sql_drop = "DROP TABLE %s"
_sql_insert_rec = "INSERT INTO record (%s) VALUES (%s)"
_sql_insert_rel = "INSERT INTO gg_release (gg_id,name) VALUES (%d, '%s')"
_sql_select_relids = "SELECT gg_id FROM gg_release WHERE name='%s'"
_sql_select_relid = """SELECT rel_id
                       FROM gg_release
                       WHERE gg_id=%d AND name='%s'"""
_sql_set_search_path = "SET search_path TO %s"

_sql_record_exists_ncbi =  """SELECT ncbi_acc_w_ver
                              FROM record
                              WHERE ncbi_acc_w_ver='%s'"""
_sql_record_exists_ggid = "SELECT gg_id FROM record WHERE gg_id=%d"

_sql_select_multiple_tax = """SELECT g.gg_id, t.tax_string
                              FROM record g INNER JOIN
                                   taxonomy t ON g.%s=t.tax_id
                              WHERE g.gg_id in (%s)"""

_sql_select_single_tax = """SELECT t.tax_string
                            FROM record g INNER JOIN
                                 taxonomy t ON g.%s=t.tax_id
                            WHERE g.gg_id=%d"""
_sql_select_max = "SELECT MAX(%s) FROM %s"

class GreengenesDB(object):
    def __init__(self, host='localhost', user='ggadmin', passwd='',
                 debug=False, database='greengenes'):
        self.con = connect(host=host, user=user, password=passwd,
                           database=database)

        if debug:
            self._set_development_schema()
            self._create_db()
            self._populate_debug_db()
        else:
            self._set_production_schema()

    def __del__(self):
        self.con.close()
        del self.con

    def _set_development_schema(self):
        """Set development schema"""
        self._execute(_sql_set_search_path % "development")

    def _set_production_schema(self):
        """Set production schema"""
        self._execute(_sql_set_search_path % "production")

    def to_arb(self, ids, aln_seq_field, directio_basename=None, size=10000):
        """Fetch ARB records

        If direct IO, data written direct to file(s). Data are written gzip'd,
        and spread over multiple files. directio_basename is the base
        filename, and that name is tagged with a unique number
        """
        bin_ids = (ids[i:i+size] for i in xrange(0, len(ids), size))

        if directio_basename is not None:
            file_count = 0
        else:
            out = []

        cursor = self.con.cursor()
        for chunk in bin_ids:
            if directio_basename is not None:
                out = gzopen(directio_basename + '_%d.txt.gz' % file_count,
                             'w')

            joined_ids = ','.join(map(str, chunk))
            cursor.execute(FULL_RECORD_DUMP % (aln_seq_field, joined_ids))

            for rec in cursor.fetchall():
                rec_lines = []
                rec_lines.append("BEGIN\n")
                for o, x in zip(FULL_RECORD_ORDER, rec):
                    if o == 'aligned_seq':
                        rec_lines.append("warning=\n")

                    if x is not None:
                        rec_lines.append("%s=%s\n" % (o, str(x)))
                    else:
                        rec_lines.append("%s=\n" % o)
                rec_lines.append("END\n\n")

                if directio_basename is not None:
                    out.write(''.join(rec_lines))
                else:
                    out.extend(rec_lines)

            if directio_basename is not None:
                out.close()
                file_count += 1

        cursor.close()
        if directio_basename is None:
            return out
        else:
            return []

    def update_greengenes_tax(self, tax_map, version):
        """Update Greengenes taxonomy fields"""
        self._update_tax('greengenes_tax_id', tax_map, version)

    def update_ncbi_tax(self, tax_map, version='NA'):
        """Update NCBI taxonomy fields"""
        self._update_tax('ncbi_tax_id', tax_map, version)

    def _update_tax(self, col, taxmap, version):
        """Insert taxonomy records, update greengenes records"""
        tax_id = self._get_max_taxid() + 1
        for gg_id, tax in taxmap.iteritems():
            gg_id = int(gg_id)

            if tax is None:
                self._execute(_sql_update_rec % (col, "NULL", gg_id))
            else:
                self._execute(_sql_insert_tax % (tax_id, version, tax))
                self._execute(_sql_update_rec % (col, tax_id, gg_id))
                tax_id += 1

        self.con.commit()

    def get_release(self, name):
        """Return the GG IDs associated with a release name"""
        with self._execute_and_more(_sql_select_relids % name) as cur:
            return [i[0] for i in cur.fetchall()]

    def get_ncbi_tax_multiple(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_tax('ncbi_tax_id', ggids)

    def get_greengenes_tax_multiple(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_tax('greengenes_tax_id', ggids)

    def _get_multiple_tax(self, field, ggids):
        """Get multiple taxonomy strings by GGIDs"""
        res = {int(i): None for i in ggids}
        ggids = ",".join(map(str, ggids))
        sql = _sql_select_multiple_tax % (field, ggids)

        with self._execute_and_more(sql) as cur:
            res.update(dict(cur.fetchall()))

        return res

    def get_ncbi_tax(self, ggid):
        """Get a taxonomy string by GGID"""
        return self._get_single_tax("ncbi_tax_id", ggid)

    def get_greengenes_tax(self, ggid):
        """Get a taxonomy string by GGID"""
        return self._get_single_tax("greengenes_tax_id", ggid)

    def _get_single_tax(self, field, ggid):
        """Get a single taxonomy string by ggid"""
        sql = _sql_select_single_tax % (field, int(ggid))

        with self._execute_and_more(sql) as cur:
            res = cur.fetchone()

        if res is None:
            return None
        else:
            return res[0]

    def get_pynast_seq(self, gg_id):
        """Get a single PyNAST sequence by GG_ID"""
        return self._get_seq("pynast_aligned_seq_id", gg_id)

    def get_ssualign_seq(self, gg_id):
        """Get a single SSU Align sequence by GG_ID"""
        return self._get_seq("aligned_seq_id", gg_id)

    def get_unaligned_seq(self, gg_id):
        """Get a single unaligned sequence by GG_ID"""
        return self._get_seq("unaligned_seq_id", gg_id)

    def _get_seq(self, field, gg_id):
        """Get a sequence by GG ID

        Returns None if not found or GG_ID doesn't exist
        """
        sql = """SELECT s.sequence
                 FROM record g INNER JOIN
                      sequence s ON g.%s=s.seq_id
                 WHERE g.gg_id=%d""" % (field, int(gg_id))

        with self._execute_and_more(sql) as cur:
            if cur.rowcount == 0:
                return None
            else:
                return cur.fetchone()[0]

    def _update_seq(self, seqs, col):
        """update greengenes record"""
        seq_id = self._get_max_seqid() + 1
        for gg_id, seq in seqs.items():
            gg_id = int(gg_id)

            self._execute(_sql_insert_seq % (seq_id, seq))
            self._execute(_sql_update_rec % (col, seq_id, gg_id))
            seq_id += 1

        self.con.commit()

    def update_pynast_seq(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq(seqs, "pynast_aligned_seq_id")

    def update_ssualign_seq(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq(seqs, "aligned_seq_id")

    def update_unaligned_seq(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq(seqs, "unaligned_seq_id")

    def _get_max_otu_cluster_id(self):
        """Returns the max observed OTU cluster ID"""
        sql = _sql_select_max % ("cluster_id", "otu_cluster")
        with self._execute_and_more(sql) as cur:
            id_ = cur.fetchone()[0]
            return id_ if id_ is not None else 0

    def _get_max_ggid(self):
        """Returns the max observed gg id"""
        sql = _sql_select_max % ("gg_id", "record")
        with self._execute_and_more(sql) as cur:
            return cur.fetchone()[0]

    def _get_max_taxid(self):
        """Returns the max observed taxonomy id"""
        sql = _sql_select_max % ("tax_id", "taxonomy")
        with self._execute_and_more(sql) as cur:
            return cur.fetchone()[0]

    def _get_max_seqid(self):
        """Returns the max observed sequence id"""
        sql = _sql_select_max % ("seq_id", "sequence")
        with self._execute_and_more(sql) as cur:
            return cur.fetchone()[0]

    def __contains__(self, item):
        ncbi = _sql_record_exists_ncbi % str(item)

        try:
            item = int(item)
            ggid = _sql_record_exists_ggid % item
        except ValueError:
            ggid = None

        if ggid is not None:
            with self._execute_and_more(ggid) as cur:
                ggid_cnt = cur.rowcount
        else:
            ggid_cnt = 0

        with self._execute_and_more(ncbi) as cur:
            ncbi_cnt = cur.rowcount

        if ggid_cnt != 0 or ncbi_cnt != 0:
            return True
        else:
            return False

    @contextmanager
    def _execute_and_more(self, sql):
        """Execute, rollback if we hit an error, otherwise get a cursor"""
        with self.con.cursor() as cursor:
            try:
                _ = cursor.execute(sql)
            except ProgrammingError:
                self.con.rollback()
                raise ValueError("Unable to execute:\n%s!" % sql)
            except OperationalError:
                self.con.rollback()
                raise ValueError("Bad value in:\n%s!" % sql)
            yield cursor

    def _execute(self, sql):
        """Execute and rollback if we hit an error"""
        with self._execute_and_more(sql):
            pass

    @staticmethod
    def _build_rec(dbrec):
        """Rebuild a record from db select results"""
        return {k: v for k, v in izip(SINGLE_RECORD_ORDER, dbrec)}

    def select_record(self, id_):
        """Return a record from the db by greengenes id"""
        if id_ not in self:
            raise ValueError("%d doesn't appear in the db!" % id_)

        with self._execute_and_more(SINGLE_RECORD % (id_, id_)) as cur:
            return self._build_rec(cur.fetchone())

    def insert_sequence(self, seq):
        """Load a sequence, return seq_id"""
        seq_id = self._get_max_seqid() + 1

        self._execute(_sql_insert_seq % (seq_id, seq))

        self.con.commit()
        return seq_id

    def insert_taxonomy(self, tax, tax_version):
        """Load a taxonomy string, return tax_id"""
        tax_id = self._get_max_taxid() + 1

        self._execute(_sql_insert_tax % (tax_id, tax_version, tax))

        self.con.commit()
        return tax_id

    def insert_record(self, record, releasename="in_holding"):
        """Insert a Greengenes record"""
        if record['ncbi_acc_w_ver'] in self:
            raise ValueError("record %s exists!" % record['ncbi_acc_w_ver'])

        ggid = self._get_max_ggid() + 1

        record['gg_id'] = ggid

        colnames = ','.join(FULL_GG_ORDER)

        vals = []
        for c in FULL_GG_ORDER:
            val = record.get(c, "NULL")
            if not val:
                val = "NULL"
            if val != "NULL":
                val = "'%s'" % str(val).replace("'", "\\'")
            vals.append(val)
        vals = ','.join(vals)

        self._execute(_sql_insert_rec % (colnames, vals))
        self._execute(_sql_insert_rel % (ggid, releasename))

        self.con.commit()

        return ggid

    def insert_otu(self, rep_id, members, method, similarity, rel_name):
        """Insert an OTU

        rep_id : a gg_id
        rel_name : a release name
        members : a list of gg_id
        method : a string < 16 bytes
        similarity : a float
        """
        if rep_id not in members:
            members.append(rep_id)

        sql = _sql_select_relid % (rep_id, rel_name)
        with self._execute_and_more(sql) as cur:
            rel_id = cur.fetchone()[0]

        if rel_id is None:
            raise ValueError("%d doesn't appear to be in %s" %
                             (rep_id, rel_name))

        c_id = self._get_max_otu_cluster_id() + 1

        sql = _sql_insert_otu_cluster % (c_id, rep_id, rel_id, similarity,
                                         method)
        self._execute(sql)

        for member in members:
            self._execute(_sql_insert_otu % (c_id, member))

        self.con.commit()

    def _create_db(self, schema='development'):
        """Create a small test database"""
        cursor = self.con.cursor()
        cursor.execute("DROP TABLE IF EXISTS %s.otu" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.otu_cluster" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.gg_release" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.chimera" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.record" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.taxonomy" % schema)
        cursor.execute("DROP TABLE IF EXISTS %s.sequence" % schema)

        cursor.execute("""
            CREATE TABLE %s.taxonomy(
            tax_id INT NOT NULL,
            tax_version VARCHAR(16) NOT NULL,
            tax_string VARCHAR(500) NOT NULL,
            PRIMARY KEY(tax_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""
            CREATE TABLE %s.sequence(
            seq_id INT NOT NULL,
            sequence VARCHAR(20000),
            PRIMARY KEY(seq_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""
            CREATE TABLE %s.record (
                gg_id INT NOT NULL,
                ncbi_acc_w_ver VARCHAR(20) NOT NULL,
                ncbi_gi INT NULL,
                db_name VARCHAR(20) NULL,
                gold_id VARCHAR(16) NULL,
                decision VARCHAR(15) NOT NULL,
                prokmsaname VARCHAR(1000) NULL,
                isolation_source VARCHAR(1000) NULL,
                clone VARCHAR(300) NULL,
                organism VARCHAR(100) NULL,
                strain VARCHAR(300) NULL,
                specific_host VARCHAR(1200) NULL,
                authors VARCHAR(10000) NULL,
                title VARCHAR(1000) NULL,
                journal VARCHAR(400) NULL,
                pubmed INT NULL,
                /*STUDY_ID INT NOT NULL,*/
                submit_date VARCHAR(25) NULL,
                country VARCHAR(400) NULL,
                ncbi_tax_id INT NULL,
                silva_tax_id INT NULL,
                greengenes_tax_id INT NULL,
                hugenholtz_tax_id INT NULL,
                non_acgt_percent FLOAT NULL,
                perc_ident_to_invariant_core FLOAT NULL,
                unaligned_seq_id INT NULL,
                aligned_seq_id INT NULL,
                pynast_aligned_seq_id int null,
                max_non_acgt_streak int null,
                PRIMARY KEY(gg_id),
                FOREIGN KEY(ncbi_tax_id) REFERENCES taxonomy(tax_id),
                FOREIGN KEY(hugenholtz_tax_id) REFERENCES taxonomy(tax_id),
                FOREIGN KEY(silva_tax_id) REFERENCES taxonomy(tax_id),
                FOREIGN KEY(greengenes_tax_id) REFERENCES taxonomy(tax_id),
                FOREIGN KEY(unaligned_seq_id) REFERENCES sequence(seq_id),
                FOREIGN KEY(aligned_seq_id) REFERENCES sequence(seq_id),
                FOREIGN KEY(pynast_aligned_seq_id) REFERENCES sequence(seq_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""
            CREATE TABLE %s.gg_release(
            rel_id SERIAL NOT NULL,
            name VARCHAR(20) NOT NULL,
            gg_id INT NOT NULL,
            PRIMARY KEY(rel_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""CREATE TABLE %s.otu_cluster (
            cluster_id INT NOT NULL,
            rep_id INT NOT NULL,
            rel_id INT NOT NULL,
            similarity FLOAT NOT NULL,
            method VARCHAR(16),
            PRIMARY KEY(cluster_id),
            FOREIGN KEY(rep_id) REFERENCES record(gg_id),
            FOREIGN KEY(rel_id) REFERENCES gg_release(rel_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""CREATE TABLE %s.otu (
            otu_id SERIAL NOT NULL,
            cluster_id INT NOT NULL,
            gg_id INT NOT NULL,
            PRIMARY KEY(otu_id),
            FOREIGN KEY(cluster_id) REFERENCES otu_cluster(cluster_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

        cursor.execute("""
            CREATE TABLE %s.chimera(
            chim_id SERIAL NOT NULL,
            gg_id INT NOT NULL,
            reason VARCHAR(1000),
            PRIMARY KEY(chim_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

    def _populate_debug_db(self):
        """Source a subset from production db"""
        cursor = self.con.cursor()
        cursor.execute("""INSERT INTO development.sequence
                          SELECT s.seq_id,s.sequence
                          FROM production.record g
                           INNER JOIN production.sequence s
                           ON g.pynast_aligned_seq_id=s.seq_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.sequence
                          SELECT s.seq_id,s.sequence
                          FROM production.record g
                           INNER JOIN production.sequence s
                           ON g.unaligned_seq_id=s.seq_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.sequence
                          SELECT s.seq_id,s.sequence
                          FROM production.record g
                           INNER JOIN production.sequence s
                           ON g.aligned_seq_id=s.seq_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.taxonomy
                          SELECT t.tax_id, t.tax_version, t.tax_string
                          FROM production.record g
                           INNER JOIN production.taxonomy t
                           ON g.ncbi_tax_id=t.tax_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.taxonomy
                          SELECT t.tax_id, t.tax_version, t.tax_string
                          FROM production.record g
                           INNER JOIN production.taxonomy t
                           ON g.hugenholtz_tax_id=t.tax_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.taxonomy
                          SELECT t.tax_id, t.tax_version, t.tax_string
                          FROM production.record g
                           INNER JOIN production.taxonomy t
                           ON g.greengenes_tax_id=t.tax_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.taxonomy
                          SELECT t.tax_id, t.tax_version, t.tax_string
                          FROM production.record g
                           INNER JOIN production.taxonomy t
                           ON g.silva_tax_id=t.tax_id
                          WHERE g.gg_id < 100""")
        cursor.execute("""INSERT INTO development.record
                          SELECT   GG_ID,
                                   NCBI_ACC_W_VER,
                                   NCBI_GI,
                                   DB_NAME ,
                                   GOLD_ID,
                                   DECISION,
                                   PROKMSANAME,
                                   ISOLATION_SOURCE,
                                   CLONE,
                                   ORGANISM,
                                   STRAIN,
                                   SPECIFIC_HOST,
                                   AUTHORS,
                                   TITLE,
                                   JOURNAL,
                                   PUBMED,
                                   SUBMIT_DATE,
                                   COUNTRY,
                                   NCBI_TAX_ID,
                                   SILVA_TAX_ID,
                                   GREENGENES_TAX_ID,
                                   HUGENHOLTZ_TAX_ID,
                                   NON_ACGT_PERCENT,
                                   PERC_IDENT_TO_INVARIANT_CORE,
                                   UNALIGNED_SEQ_ID,
                                   ALIGNED_SEQ_ID,
                                   pynast_aligned_seq_id,
                                   max_non_acgt_streak
                          FROM production.record
                          WHERE gg_id < 100""")
        cursor.execute("""INSERT INTO development.gg_release
                          SELECT r.rel_id, r.name, r.gg_id
                          FROM production.record g
                           INNER JOIN production.gg_release r
                           ON g.gg_id=r.gg_id
                          WHERE g.gg_id < 100""")
        self.con.commit()
