#!/usr/bin/env python

from psycopg2 import connect, ProgrammingError, OperationalError
from random import choice
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

FULL_RECORD_DUMP = """SELECT gg_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,decision,prokmsaname,isolation_source,clone,organism,strain,specific_host,authors,title,journal,pubmed,submit_date,country,nt.tax_string AS ncbi_tax_string,st.tax_string AS silva_tax_string,gg.tax_string AS greengenes_tax_string,h.tax_string AS hugenholtz_tax_id,non_acgt_percent,perc_ident_to_invariant_core,aseq.sequence AS aligned_seq
FROM record g
LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id
LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id
LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id
LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id
LEFT JOIN sequence aseq ON aseq.seq_id=g.%s
WHERE g.gg_id IN (%s)"""

SINGLE_RECORD = """SELECT gg_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,decision,prokmsaname,isolation_source,clone,organism,strain,specific_host,authors,title,journal,pubmed,submit_date,country,nt.tax_string AS ncbi_tax_string,st.tax_string AS silva_tax_string,gg.tax_string AS greengenes_tax_string,h.tax_string AS hugenholtz_tax_id,non_acgt_percent,perc_ident_to_invariant_core,ssu.sequence AS ssualign_seq, pynast.sequence AS pynast_seq, unaligned.sequence AS unaligned_seq
FROM record g
LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id
LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id
LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id
LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id
LEFT JOIN sequence ssu ON ssu.seq_id=g.aligned_seq_id
LEFT JOIN sequence pynast on pynast.seq_id=g.pynast_aligned_seq_id
LEFT JOIN sequence unaligned on unaligned.seq_id=g.unaligned_seq_id
WHERE g.gg_id='%s' or g.ncbi_acc_w_ver='%s'"""

SINGLE_RECORD_ORDER = ['gg_id','ncbi_acc_w_ver','ncbi_gi','db_name','gold_id','decision','prokmsaname','isolation_source','clone','organism','strain','specific_host','authors','title','journal','pubmed','submit_date','country','ncbi_tax_string','silva_tax_string','greengenes_tax_string','hugenholtz_tax_id','non_acgt_percent','perc_ident_to_invariant_core','ssualign_seq','pynast_seq','unaligned_seq']

FULL_RECORD_ORDER = ["gg_id","ncbi_acc_w_ver","ncbi_gi","db_name",
                     "gold_id","decision","prokmsaname","isolation_source",
                     "clone","organism","strain","specific_host","authors",
                     "title","journal","pubmed","submit_date","country",
                     "ncbi_tax_string","silva_tax_string",
                     "greengenes_tax_string","hugenholtz_tax_id",
                     "non_acgt_percent","perc_ident_to_invariant_core",
                     "aligned_seq"]

FULL_GG_ORDER = ["gg_id","ncbi_acc_w_ver","ncbi_gi","db_name",
                 "gold_id","decision","prokmsaname","isolation_source",
                 "clone","organism","strain","specific_host","authors",
                 "title","journal","pubmed","submit_date","country",
                 "ncbi_tax_id","silva_tax_id",
                 "greengenes_tax_id","hugenholtz_tax_id",
                 "non_acgt_percent","perc_ident_to_invariant_core",
                 "unaligned_seq_id", "aligned_seq_id",
                 "pynast_aligned_seq_id"]

ALPHA = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

_sql_insert_tax = """INSERT INTO taxonomy(tax_id, tax_version, tax_string)
                     VALUES (%d, '%s', '%s')"""
_sql_insert_seq = """INSERT INTO sequence(seq_id, sequence)
                     VALUES (%d, '%s')"""
_sql_update_rec = """UPDATE record
                     SET %s=%s
                     WHERE gg_id=%d"""
_sql_insert_otu_cluster = """
                 INSERT INTO otu_cluster(cluster_id, rep_id, similarity, method)
                 VALUES (%d, %d, %f, '%s')"""
_sql_insert_otu = """INSERT INTO otu(cluster_id, gg_id)
                     VALUES (%d, %d)"""
_sql_create_tmp = "CREATE TEMPORARY TABLE %s LIKE %s"""
_sql_drop =       "DROP TABLE %s"
_sql_insert_rec = "INSERT INTO record (%s) VALUES (%s)"
_sql_insert_rel = "INSERT INTO gg_release (gg_id,name) VALUES (%d, '%s')"
_sql_select_relids = "SELECT gg_id FROM gg_release WHERE name='%s'"
_sql_set_search_path = "SET search_path TO %s"

#__sql_select_rec = "SELECT %s FROM greengenes" % (",".join(FULL_GG_ORDER))
#_sql_select_rec_ggid = __sql_select_rec + " WHERE gg_id=%d"
#_sql_select_rec_ncbi = __sql_select_rec + " WHERE ncbi_acc_w_ver=%d"
#_sql_select_rec_both = __sql_select_rec + """ WHERE gg_id='%s' OR
#                                                ncbi_acc_w_ver='%s'"""

class GreengenesDB(object):

    def __init__(self, host='localhost', user='ggadmin', passwd='',
                 debug=False, database='greengenes'):
        self.con = connect(host=host, user=user, password=passwd,
                           database=database)
        self.cursor = self.con.cursor()

        if debug:
            self._set_development_schema()
            self._create_db()
            self._populate_debug_db()

    def __del__(self):
        del self.cursor
        del self.con

    def _set_development_schema(self):
        self._execute_safe(_sql_set_search_path % "development")

    def bulkFetchARBRecords(self, ids, aln_seq_field, directio_basename=None,
                            size=10000):
        """Bulk fetch ARB records

        If direct IO, data written direct to file(s). Data are written gzip'd,
        and spread over multiple files. directio_basename is the base
        filename, and that name is tagged with a unique number
        """
        bin_ids = (ids[i:i+size] for i in xrange(0, len(ids), size))

        if directio_basename is not None:
            file_count = 0
        else:
            out = []

        for chunk in bin_ids:
            if directio_basename is not None:
                out = gzopen(directio_basename + '_%d.txt.gz' % file_count,
                             'w')

            joined_ids = ','.join(map(str, chunk))
            self.cursor.execute(FULL_RECORD_DUMP % (aln_seq_field, joined_ids))

            for rec in self.cursor.fetchall():
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

        if directio_basename is None:
            return out
        else:
            return []

    def updateGreengenesTaxonomy(self, tax_map, version):
        """Update Greengenes taxonomy fields"""
        self._update_tax_field('greengenes_tax_id', tax_map, version)

    def updateNCBITaxonomy(self, tax_map, version='NA'):
        """Update NCBI taxonomy fields"""
        self._update_tax_field('ncbi_tax_id', tax_map, version)

    def _update_tax_field(self, col, taxmap, version):
        """Insert taxonomy records, update greengenes records"""
        tax_id = self._get_max_taxid() + 1
        for gg_id, tax in taxmap.iteritems():
            gg_id = int(gg_id)

            if tax is None:
                self._execute_safe(_sql_update_rec % (col, "NULL", gg_id))
            else:
                self._execute_safe(_sql_insert_tax % (tax_id, version, tax))
                self._execute_safe(_sql_update_rec % (col, tax_id, gg_id))
                tax_id += 1

        self.con.commit()

    def getReleaseGGIDs(self, name):
        """Return the IDs associated with a release name"""
        self._execute_safe(_sql_select_relids % name)
        return [i[0] for i in self.cursor.fetchall()]

    def getNCBITaxonomyMultipleGGID(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_taxonomy_strings_ggid('ncbi_tax_id', ggids)

    def getGreengenesTaxonomyMultipleGGID(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_taxonomy_strings_ggid('greengenes_tax_id',
                                                        ggids)

    def _get_multiple_taxonomy_strings_ggid(self, field, ggids):
        """Get multiple taxonomy strings by GGIDs"""
        res = dict([(int(i), None) for i in ggids])
        ggids = ",".join(map(str, ggids))
        try:
            self.cursor.execute("""
                SELECT g.gg_id, t.tax_string
                FROM record g INNER JOIN
                     taxonomy t ON g.%s=t.tax_id
                WHERE g.gg_id in (%s)""" % (field, ggids))
        except ProgrammingError:
            raise ValueError("Failed query!")

        res.update(dict(self.cursor.fetchall()))
        return res

    def getNCBITaxonomySingleGGID(self, ggid):
        """Get a taxonomy string by GGID"""
        return self._get_single_taxonomy_string_ggid("ncbi_tax_id", ggid)

    def getGreengenesTaxonomySingleGGID(self, ggid):
        """Get a taxonomy string by GGID"""
        return self._get_single_taxonomy_string_ggid("greengenes_tax_id", ggid)

    def _get_single_taxonomy_string_ggid(self, field, ggid):
        """Get a single taxonomy string by ggid"""
        try:
            self.cursor.execute("""
                SELECT t.tax_string
                FROM record g INNER JOIN
                     taxonomy t ON g.%s=t.tax_id
                WHERE g.gg_id=%d""" % (field, int(ggid)))
        except ProgrammingError:
            raise ValueError("Failed query!")

        res = self.cursor.fetchone()
        if res is None:
            return None
        else:
            return res[0]

    def getPyNASTSequenceGGID(self, gg_id):
        """Get a single PyNAST sequence by GG_ID"""
        return self._get_sequence_gg_id("pynast_aligned_seq_id", gg_id)

    def getSSUAlignSequenceGGID(self, gg_id):
        """Get a single SSU Align sequence by GG_ID"""
        return self._get_sequence_gg_id("aligned_seq_id", gg_id)

    def getUnalignedSequenceGGID(self, gg_id):
        """Get a single unaligned sequence by GG_ID"""
        return self._get_sequence_gg_id("unaligned_seq_id", gg_id)

    def _get_sequence_gg_id(self, field, gg_id):
        """Get a sequence by GG ID

        Returns None if not found or GG_ID doesn't exist
        """
        try:
            n = self.cursor.execute("""SELECT s.sequence
                                       FROM record g INNER JOIN
                                            sequence s ON g.%s=s.seq_id
                                       WHERE g.gg_id=%d""" %(field,int(gg_id)))
        except ProgrammingError:
            return None

        if n == 0:
            return None
        else:
            return self.cursor.fetchone()[0]

    def _update_seq_field(self, seqs, col):
        """update greengenes record"""
        seq_id = self._get_max_seqid() + 1
        for gg_id, seq in seqs.items():
            gg_id = int(gg_id)

            self._execute_safe(_sql_insert_seq % (seq_id, seq))
            self._execute_safe(_sql_update_rec % (col, seq_id, gg_id))
            seq_id += 1

        self.con.commit()

    def updatePyNASTSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "pynast_aligned_seq_id")

    def updateSSUAlignSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "aligned_seq_id")

    def updateUnalignedSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "unaligned_seq_id")

    def _get_max_otu_cluster_id(self):
        """Returns the max observed OTU cluster ID"""
        query = "SELECT MAX(cluster_id) FROM otu_cluster"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        id_ = self.cursor.fetchone()[0]
        return id_ if id_ is not None else 0

    def _get_max_ggid(self):
        """Returns the max observed gg id"""
        query = "SELECT MAX(gg_id) FROM record"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    def _get_max_taxid(self):
        """Returns the max observed taxonomy id"""
        query = "SELECT MAX(tax_id) FROM taxonomy"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    def _get_max_seqid(self):
        """Returns the max observed sequence id"""
        query = "SELECT MAX(seq_id) FROM sequence"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    ### can we drop?
    def _create_tmp_sequence_table(self):
        """ """
        name = self._get_tmp_table_name("tmp_seq")
        self._execute_safe(_sql_create_tmp % (name, "SEQUENCE"))
        self.con.commit()

        return name

    ### can we drop?
    def _drop_tmp_sequence_table(self, name):
        """ """
        if not name.startswith('tmp_'):
            raise ValueError("%s doesn't appear to be temporary!" % name)

        self._execute_safe(_sql_drop % name)

        self.con.commit()

    ### can we drop?
    def _get_tmp_table_name(self, tag='tmp_'):
        """Returns a random table name"""
        if not tag.startswith('tmp'):
            tag = 'tmp_' + tag

        name = None
        while name is None:
            name = "%s_%s" % (tag, ''.join([choice(ALPHA) for i in range(5)]))
            self.cursor.execute("""
                SELECT table_name
                FROM information_schema.tables
                WHERE table_name='%s'""" % name)
            if self.cursor.fetchall():
                name = None
        return name

    def __contains__(self, item):
        ncbi = """select ncbi_acc_w_ver
                  from record
                  where ncbi_acc_w_ver='%s'""" % str(item)
        try:
            item = int(item)
            ggid = "select gg_id from record where gg_id=%d" % item
        except ValueError:
            ggid = None

        if ggid is not None:
            self.cursor.execute(ggid)
            ggid_cnt = self.cursor.rowcount
        else:
            ggid_cnt = 0

        self.cursor.execute(ncbi)
        ncbi_cnt = self.cursor.rowcount

        if ggid_cnt != 0 or ncbi_cnt != 0:
            return True
        else:
            return False

    def _execute_safe(self, sql):
        try:
            self.cursor.execute(sql)
        except ProgrammingError:
            self.con.rollback()
            raise ValueError("Unable to execute:\n%s!" % sql)
        except OperationalError:
            self.con.rollback()
            raise ValueError("Bad value in:\n%s!" % sql)

    def _build_rec(self, dbrec):
        """Rebuild a record from db select results"""
        return {k: v for k, v in izip(SINGLE_RECORD_ORDER, dbrec)}

    def selectRecord(self, id_):
        """Return a record from the db"""
        if id_ not in self:
            raise ValueError("%d doesn't appear in the db!" % id_)

        self._execute_safe(SINGLE_RECORD % (id_, id_))
        return self._build_rec(self.cursor.fetchone())

    def insertSequence(self, seq):
        """Load a sequence, return seq_id"""
        seq_id = self._get_max_seqid() + 1

        self._execute_safe(_sql_insert_seq % (seq_id, seq))

        self.con.commit()
        return seq_id

    def insertTaxonomy(self, tax, tax_version):
        """Load a taxonomy string, return tax_id"""

        tax_id = self._get_max_taxid() + 1

        self._execute_safe(_sql_insert_tax % (tax_id, tax_version, tax))

        self.con.commit()
        return tax_id

    def insertRecord(self, record, releasename="in_holding"):
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

        self._execute_safe(_sql_insert_rec % (colnames, vals))
        self._execute_safe(_sql_insert_rel % (ggid, releasename))

        self.con.commit()

        return ggid

    def insertOTU(self, rep, members, method, similarity):
        """Insert an OTU

        rep : a gg_id
        members : a list of gg_id
        method : a string < 16 bytes
        similarity : a float
        """
        if rep not in members:
            members.append(rep)

        c_id = self._get_max_otu_cluster_id() + 1

        sql = _sql_insert_otu_cluster % (c_id, rep, similarity, method)
        self._execute_safe(sql)

        for member in members:
            self._execute_safe(_sql_insert_otu % (c_id, member))

        self.con.commit()

    def release(self, prior_version, tmp_table='in_holding', thresholds=None):
        """Yield sequences tagged by decision type

        """
        if thresholds is not None:
            # expects [(column_name, op, value)]
            criteria = [op.join([c,v]) for c,op,v in thresholds]
            filter_criteria = " AND %s " % ' AND '.join(criteria)
        else:
            filter_criteria = ""

        #self._execute_safe(_sql_release_ni % (tmp_table, filter_criteria))
        #return (i for i in self.db.cursor.fetchall()):

    def _create_db(self, schema='development'):
        """Create a small test database"""
        self.cursor.execute("DROP TABLE IF EXISTS %s.otu" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.otu_cluster" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.gg_release" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.chimera" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.record" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.taxonomy" % schema)
        self.cursor.execute("DROP TABLE IF EXISTS %s.sequence" % schema)


        self.cursor.execute("""
            CREATE TABLE %s.taxonomy(
            tax_id INT NOT NULL,
            tax_version VARCHAR(16) NOT NULL,
            tax_string VARCHAR(500) NOT NULL,
            PRIMARY KEY(tax_id)
            )""" % schema)

        self.cursor.execute("""
            CREATE TABLE %s.sequence(
            seq_id INT NOT NULL,
            sequence VARCHAR(20000),
            PRIMARY KEY(seq_id)
            )""" % schema)
        self.con.commit()

        self.cursor.execute("""
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

        self.cursor.execute("""
            CREATE TABLE %s.gg_release(
            rel_id SERIAL NOT NULL,
            name VARCHAR(20) NOT NULL,
            gg_id INT NOT NULL,
            PRIMARY KEY(rel_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

        self.cursor.execute("""CREATE TABLE %s.otu_cluster (
            cluster_id INT NOT NULL,
            rep_id INT NOT NULL,
            similarity FLOAT NOT NULL,
            method VARCHAR(16),
            PRIMARY KEY(cluster_id),
            FOREIGN KEY(rep_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

        self.cursor.execute("""CREATE TABLE %s.otu (
            otu_id SERIAL NOT NULL,
            cluster_id INT NOT NULL,
            gg_id INT NOT NULL,
            PRIMARY KEY(otu_id),
            FOREIGN KEY(cluster_id) REFERENCES otu_cluster(cluster_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

        self.cursor.execute("""
            CREATE TABLE %s.chimera(
            chim_id SERIAL NOT NULL,
            gg_id INT NOT NULL,
            reason VARCHAR(1000),
            PRIMARY KEY(chim_id),
            FOREIGN KEY(gg_id) REFERENCES record(gg_id)
            )""" % schema)
        self.con.commit()

    def _populate_debug_db(self):
        self.cursor.execute("""INSERT INTO development.sequence
                               SELECT s.seq_id,s.sequence
                               FROM production.record g
                                INNER JOIN production.sequence s
                                ON g.pynast_aligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.sequence
                               SELECT s.seq_id,s.sequence
                               FROM production.record g
                                INNER JOIN production.sequence s
                                ON g.unaligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.sequence
                               SELECT s.seq_id,s.sequence
                               FROM production.record g
                                INNER JOIN production.sequence s
                                ON g.aligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM production.record g
                                INNER JOIN production.taxonomy t
                                ON g.ncbi_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM production.record g
                                INNER JOIN production.taxonomy t
                                ON g.hugenholtz_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM production.record g
                                INNER JOIN production.taxonomy t
                                ON g.greengenes_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM production.record g
                                INNER JOIN production.taxonomy t
                                ON g.silva_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO development.record
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
        self.cursor.execute("""INSERT INTO development.gg_release
                               SELECT r.rel_id, r.name, r.gg_id
                               FROM production.record g
                                INNER JOIN production.gg_release r
                                ON g.gg_id=r.gg_id
                               WHERE g.gg_id < 100""")
        self.con.commit()
