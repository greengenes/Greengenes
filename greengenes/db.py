#!/usr/bin/env python

from MySQLdb import (connect, ProgrammingError, OperationalError, DataError)
from random import choice
from gzip import open as gzopen

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2013, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "0.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

FULL_RECORD_DUMP = """SELECT gg_id,prokmsa_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,decision,prokmsaname,isolation_source,clone,organism,strain,specific_host,authors,title,journal,pubmed,submit_date,country,nt.tax_string AS ncbi_tax_string,st.tax_string AS silva_tax_string,gg.tax_string AS greengenes_tax_string,h.tax_string AS hugenholtz_tax_id,non_acgt_percent,perc_ident_to_invariant_core,small_gap_intrusions,bel3_div_ratio,bel3_a,bel3_b,b3a.tax_string AS bel3_a_tax,b3b.tax_string AS bel3_b_tax,chim_slyr_a,chim_slyr_b,csa.tax_string AS chim_slyr_a_tax,csb.tax_string AS chim_slyr_b_tax,aseq.sequence AS aligned_seq
FROM greengenes g
LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id
LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id
LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id
LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id
LEFT JOIN taxonomy b3a ON b3a.tax_id=g.bel3_a_tax_id
LEFT JOIN taxonomy b3b ON b3b.tax_id=g.bel3_b_tax_id
LEFT JOIN taxonomy csa ON csa.tax_id=g.chim_slyr_a_tax_id
LEFT JOIN taxonomy csb ON csb.tax_id=g.chim_slyr_b_tax_id
LEFT JOIN sequence aseq ON aseq.seq_id=g.%s
WHERE g.gg_id IN (%s)"""

SINGLE_RECORD = """SELECT gg_id,prokmsa_id,ncbi_acc_w_ver,ncbi_gi,db_name,gold_id,decision,prokmsaname,isolation_source,clone,organism,strain,specific_host,authors,title,journal,pubmed,submit_date,country,nt.tax_string AS ncbi_tax_string,st.tax_string AS silva_tax_string,gg.tax_string AS greengenes_tax_string,h.tax_string AS hugenholtz_tax_id,non_acgt_percent,perc_ident_to_invariant_core,small_gap_intrusions,ssu.sequence AS ssualign_seq, pynast.sequence AS pynast_seq, unaligned.sequence AS unaligned_seq
FROM greengenes g
LEFT JOIN taxonomy st ON st.tax_id=g.silva_tax_id
LEFT JOIN taxonomy nt ON nt.tax_id=g.ncbi_tax_id
LEFT JOIN taxonomy h ON h.tax_id=g.hugenholtz_tax_id
LEFT JOIN taxonomy gg ON gg.tax_id=g.greengenes_tax_id
LEFT JOIN sequence ssu ON ssu.seq_id=g.aligned_seq_id
LEFT JOIN sequence pynast on pynast.seq_id=g.pynast_aligned_seq_id
LEFT JOIN sequence unaligned on unaligned.seq_id=g.unaligned_seq_id
WHERE g.gg_id='%s' or g.ncbi_acc_w_ver='%s'"""

SINGLE_RECORD_ORDER = ['gg_id','prokmsa_id','ncbi_acc_w_ver','ncbi_gi','db_name','gold_id','decision','prokmsaname','isolation_source','clone','organism','strain','specific_host','authors','title','journal','pubmed','submit_date','country','ncbi_tax_string','silva_tax_string','greengenes_tax_string','hugenholtz_tax_id','non_acgt_percent','perc_ident_to_invariant_core','small_gap_intrusions','ssualign_seq','pynast_seq','unaligned_seq']

FULL_RECORD_ORDER = ["gg_id","prokmsa_id","ncbi_acc_w_ver","ncbi_gi","db_name",
                     "gold_id","decision","prokmsaname","isolation_source",
                     "clone","organism","strain","specific_host","authors",
                     "title","journal","pubmed","submit_date","country",
                     "ncbi_tax_string","silva_tax_string",
                     "greengenes_tax_string","hugenholtz_tax_id",
                     "non_acgt_percent","perc_ident_to_invariant_core",
                     "small_gap_intrusions","bel3_div_ratio","bel3_a","bel3_b",
                     "bel3_a_tax","bel3_b_tax","chim_slyr_a","chim_slyr_b",
                     "chim_slyr_a_tax","chim_slyr_b_tax","aligned_seq"]

FULL_GG_ORDER = ["gg_id","prokmsa_id","ncbi_acc_w_ver","ncbi_gi","db_name",
                 "gold_id","decision","prokmsaname","isolation_source",
                 "clone","organism","strain","specific_host","authors",
                 "title","journal","pubmed","submit_date","country",
                 "ncbi_tax_id","silva_tax_id",
                 "greengenes_tax_id","hugenholtz_tax_id",
                 "non_acgt_percent","perc_ident_to_invariant_core",
                 "small_gap_intrusions", "unaligned_seq_id", "aligned_seq_id",
                 "pynast_aligned_seq_id"]

ALPHA = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'

_sql_insert_tax = """INSERT INTO taxonomy(tax_id, tax_version, tax_string)
                     VALUES (%d, '%s', '%s')"""
_sql_insert_seq = """INSERT INTO sequence(seq_id, sequence)
                     VALUES (%d, '%s')"""
_sql_update_rec = """UPDATE greengenes
                     SET %s=%s
                     WHERE gg_id=%d"""
_sql_insert_otu_cluster = """
                 INSERT INTO otu_cluster(cluster_id, rep_id, similarity, method)
                 VALUES (%d, %d, %f, '%s')"""
_sql_insert_otu = """INSERT INTO otu(cluster_id, gg_id)
                     VALUES (%d, %d)"""
_sql_create_tmp = "CREATE TEMPORARY TABLE %s LIKE %s"""
_sql_drop =       "DROP TABLE %s"
_sql_insert_rec = "INSERT INTO greengenes (%s) VALUES (%s)"
_sql_insert_rel = "INSERT INTO gg_release (gg_id,name) VALUES (%d, '%s')"

#__sql_select_rec = "SELECT %s FROM greengenes" % (",".join(FULL_GG_ORDER))
#_sql_select_rec_ggid = __sql_select_rec + " WHERE gg_id=%d"
#_sql_select_rec_ncbi = __sql_select_rec + " WHERE ncbi_acc_w_ver=%d"
#_sql_select_rec_both = __sql_select_rec + """ WHERE gg_id='%s' OR
#                                                ncbi_acc_w_ver='%s'"""

class GreengenesMySQL(object):

    def __init__(self, host='localhost',user='greengenes',passwd='',
                 debug=False):
        self.con = connect(host,user=user,passwd=passwd)
        self.cursor = self.con.cursor()

        if debug:
            self._create_debug_db()
        else:
            self.cursor.execute('USE greengenes')

        self._have_locks = False
        self._lock_aliases = {}

    def __del__(self):
        if self._have_locks:
            self._unlock()

        del self.cursor
        del self.con

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
                out = gzopen(directio_basename + '_%d.txt.gz' % file_count, 'w')

            joined_ids = ','.join(map(str, chunk))
            self.cursor.execute(FULL_RECORD_DUMP % (aln_seq_field, joined_ids))

            for rec in self.cursor.fetchall():
                rec_lines = []
                rec_lines.append("BEGIN\n")
                for o,x in zip(FULL_RECORD_ORDER, rec):
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
        self._lock([('greengenes', None, 'write'),
                    ('taxonomy', None, 'write')])

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
        self._unlock()

    def getNCBITaxonomyMultipleGGID(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_taxonomy_strings_ggid('ncbi_tax_id', ggids)

    def getGreengenesTaxonomyMultipleGGID(self, ggids):
        """Query multiple GGIDs at a time"""
        return self._get_multiple_taxonomy_strings_ggid('greengenes_tax_id', \
                ggids)

    def _get_multiple_taxonomy_strings_ggid(self, field, ggids):
        """Get multiple taxonomy strings by GGIDs"""
        res = dict([(int(i), None) for i in ggids])
        ggids = ",".join(map(str, ggids))
        try:
            n = self.cursor.execute("""
                SELECT g.gg_id, t.tax_string
                FROM greengenes g INNER JOIN
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
            n = self.cursor.execute("""
                SELECT t.tax_string
                FROM greengenes g INNER JOIN
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
                                       FROM greengenes g INNER JOIN
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
        self._lock([('greengenes', None, 'write'),
                    ('sequence', None, 'write')])

        seq_id = self._get_max_seqid() + 1
        for gg_id, seq in seqs.items():
            gg_id = int(gg_id)

            self._execute_safe(_sql_insert_seq % (seq_id, seq))
            self._execute_safe(_sql_update_rec % (col, seq_id, gg_id))
            seq_id += 1

        self.con.commit()
        self._unlock()

    def updatePyNASTSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "pynast_aligned_seq_id")

    def updateSSUAlignSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "aligned_seq_id")

    def updateUnalignedSequences(self, seqs):
        """seqs -> {gg_id:sequence}"""
        self._update_seq_field(seqs, "unaligned_seq_id")

    def _lock(self, to_lock):
        """to_lock is [(tablename, alias, locktype)]"""
        if self._have_locks:
            raise ValueError("Already have locks!")

        query_builder = []
        lock_aliases = {}
        for table, alias, locktype in to_lock:
            if alias is None:
                query_builder.append("%s %s" % (table, locktype))
            else:
                query_builder.append("%s as %s %s" % (table, alias, locktype))
            lock_aliases[table] = alias
        query = "LOCK TABLES %s" % ', '.join(query_builder)

        try:
            self.cursor.execute(query)
        except ProgrammingError:
            raise ValueError("Unable to lock %s!" %
                             ', '.join([i[0] for i in to_lock]))

        self._have_locks = True
        self._lock_aliases = lock_aliases

    def _unlock(self):
        """Unlock locked tables"""
        if self._have_locks:
            self.cursor.execute("UNLOCK TABLES")
            self._have_locks = False
            self._lock_aliases = {}

    def _get_max_otu_cluster_id(self):
        """Returns the max observed OTU cluster ID"""
        if self._have_locks:
            if 'otu_cluster' not in self._lock_aliases:
                raise ValueError("Cannot access table 'otu_cluster'")
            else:
                if self._lock_aliases['otu_cluster'] is None:
                    query = "SELECT MAX(cluster_id) FROM otu_cluster"
                else:
                    query = """SELECT MAX(cluster_id)
                               FROM otu_cluster %s""" % \
                                       self._lock_aliases['otu_cluster']
        else:
            query = "SELECT MAX(cluster_id) FROM otu_cluster"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        id_ = self.cursor.fetchone()[0]
        return id_ if id_ is not None else 0

    def _get_max_ggid(self):
        """Returns the max observed gg id"""
        if self._have_locks:
            if 'greengenes' not in self._lock_aliases:
                raise ValueError("Cannot access table 'greengenes'")
            else:
                if self._lock_aliases['greengenes'] is None:
                    query = "SELECT MAX(gg_id) FROM greengenes"
                else:
                    query = """SELECT MAX(gg_id)
                           FROM greengenes %s""" % self._lock_aliases['greengenes']
        else:
            query = "SELECT MAX(gg_id) FROM greengenes"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    def _get_max_taxid(self):
        """Returns the max observed taxonomy id"""
        if self._have_locks:
            if 'taxonomy' not in self._lock_aliases:
                raise ValueError("Cannot access table 'taxonomy'")
            else:
                if self._lock_aliases['taxonomy'] is None:
                    query = "SELECT MAX(tax_id) FROM taxonomy"
                else:
                    query = """SELECT MAX(tax_id)
                           FROM taxonomy %s""" % self._lock_aliases['taxonomy']
        else:
            query = "SELECT MAX(tax_id) FROM taxonomy"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    def _get_max_seqid(self):
        """Returns the max observed sequence id"""
        if self._have_locks:
            if 'sequence' not in self._lock_aliases:
                raise ValueError("Cannot access table 'sequence'")
            else:
                if self._lock_aliases['sequence'] is None:
                    query = "SELECT MAX(seq_id) FROM sequence"
                else:
                    query = """SELECT MAX(seq_id)
                           FROM sequence %s""" % self._lock_aliases['sequence']
        else:
            query = "SELECT MAX(seq_id) FROM sequence"

        try:
            self.cursor.execute(query)
        except:
            raise ValueError("Unable to query")

        return self.cursor.fetchone()[0]

    def _create_tmp_sequence_table(self):
        """ """
        name = self._get_tmp_table_name("tmp_seq")
        self._execute_safe(_sql_create_tmp % (name, "SEQUENCE"))
        self.con.commit()

        return name

    def _drop_tmp_sequence_table(self, name):
        """ """
        if not name.startswith('tmp_'):
            raise ValueError("%s doesn't appear to be temporary!" % name)

        self._execute_safe(_sql_drop % name)

        self.con.commit()

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
                  from greengenes
                  where ncbi_acc_w_ver='%s'""" % str(item)
        try:
            item = int(item)
            ggid = "select gg_id from greengenes where gg_id=%d" % item
        except ValueError:
            ggid = None

        if ggid is not None:
            ggid_cnt = self.cursor.execute(ggid)
        else:
            ggid_cnt = 0

        ncbi_cnt = self.cursor.execute(ncbi)

        if ggid_cnt != 0 or ncbi_cnt != 0:
            return True
        else:
            return False

    def _execute_safe(self, sql):
        try:
            self.cursor.execute(sql)
        except ProgrammingError:
            self.con.rollback()
            self._unlock()
            raise ValueError("Unable to execute:\n%s!" % sql)
        except OperationalError:
            self.con.rollback()
            self._unlock()
            raise ValueError("Bad value in:\n%s!" % sql)

    def _build_rec(self, dbrec):
        """Rebuild a record from db select results"""
        return {k:v for k,v in zip(SINGLE_RECORD_ORDER, dbrec)}

    def selectRecord(self, id_):
        """Return a record from the db"""
        if id_ not in self:
            raise ValueError("%d doesn't appear in the db!" % id_)

        self._execute_safe(SINGLE_RECORD % (id_, id_))
        return self._build_rec(self.cursor.fetchone())

    def insertSequence(self, seq):
        """Load a sequence, return seq_id"""
        self._lock([('sequence', None, 'write')])
        seq_id = self._get_max_seqid() + 1

        self._execute_safe(_sql_insert_seq % (seq_id, seq))

        self.con.commit()
        self._unlock()
        return seq_id

    def insertTaxonomy(self, tax, tax_version):
        """Load a taxonomy string, return tax_id"""
        self._lock([('taxonomy', None, 'write')])

        tax_id = self._get_max_taxid() + 1

        self._execute_safe(_sql_insert_tax % (tax_id, tax_version, tax))

        self.con.commit()
        self._unlock()
        return tax_id

    def insertRecord(self, record, releasename="in_holding"):
        """Insert a Greengenes record"""
        if record['ncbi_acc_w_ver'] in self:
            raise ValueError("record %s exists!" % record['ncbi_acc_w_ver'])

        self._lock([('greengenes', None, 'write'),
                    ('gg_release', None, 'write')])

        ggid = self._get_max_ggid() + 1

        record['gg_id'] = ggid
        record['prokmsa_id'] = ggid

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
        self._unlock()

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

        self._lock([('otu_cluster', None, 'write'),
                    ('otu', None, 'write')])

        c_id = self._get_max_otu_cluster_id() + 1

        sql = _sql_insert_otu_cluster % (c_id, rep, similarity, method)
        self._execute_safe(sql)

        for member in members:
            self._execute_safe(_sql_insert_otu % (c_id, member))

        self.con.commit()
        self._unlock()

    def _create_debug_db(self):
        """Create a small test database"""
        self.cursor.execute('USE test_greengenes')

        self.cursor.execute("DROP TABLE IF EXISTS otu")
        self.cursor.execute("DROP TABLE IF EXISTS otu_cluster")
        self.cursor.execute("DROP TABLE IF EXISTS gg_release")
        self.cursor.execute("DROP TABLE IF EXISTS greengenes")
        self.cursor.execute("DROP TABLE IF EXISTS taxonomy")
        self.cursor.execute("DROP TABLE IF EXISTS sequence")

        self.cursor.execute("""
            CREATE TABLE TAXONOMY(
            TAX_ID INT NOT NULL,
            TAX_VERSION VARCHAR(16) NOT NULL,
            TAX_STRING VARCHAR(500) NOT NULL,
            PRIMARY KEY(TAX_ID)
            )""")

        self.cursor.execute("""
            CREATE TABLE SEQUENCE(
            SEQ_ID INT NOT NULL,
            SEQUENCE VARCHAR(20000),
            PRIMARY KEY(SEQ_ID)
            )""")
        self.con.commit()

        self.cursor.execute("""
            CREATE TABLE GREENGENES (
            GG_ID INT NOT NULL,
            PROKMSA_ID INT NOT NULL,
            NCBI_ACC_W_VER VARCHAR(20) NOT NULL,
            NCBI_GI INT NULL,
            DB_NAME VARCHAR(20) NULL,
            GOLD_ID VARCHAR(16) NULL,
            DECISION VARCHAR(15) NOT NULL,
            PROKMSANAME VARCHAR(1000) NULL,
            ISOLATION_SOURCE VARCHAR(1000) NULL,
            CLONE VARCHAR(300) NULL,
            ORGANISM VARCHAR(100) NULL,
            STRAIN VARCHAR(300) NULL,
            SPECIFIC_HOST VARCHAR(1200) NULL,
            AUTHORS VARCHAR(10000) NULL,
            TITLE VARCHAR(1000) NULL,
            JOURNAL VARCHAR(400) NULL,
            PUBMED INT NULL,
            /*STUDY_ID INT NOT NULL,*/
            SUBMIT_DATE VARCHAR(25) NULL,
            COUNTRY VARCHAR(400) NULL,
            NCBI_TAX_ID INT NULL,
            SILVA_TAX_ID INT NULL,
            GREENGENES_TAX_ID INT NULL,
            HUGENHOLTZ_TAX_ID INT NULL,
            NON_ACGT_PERCENT FLOAT NULL,
            PERC_IDENT_TO_INVARIANT_CORE FLOAT NULL,
            SMALL_GAP_INTRUSIONS FLOAT NULL,
            BEL3_DIV_RATIO FLOAT NULL,
            BEL3_A VARCHAR(25) NULL,
            BEL3_B VARCHAR(25) NULL,
            BEL3_A_TAX_ID INT NULL,
            BEL3_B_TAX_ID INT NULL,
            CHIM_SLYR_A VARCHAR(25) NULL,
            CHIM_SLYR_B VARCHAR(25) NULL,
            CHIM_SLYR_A_TAX_ID INT NULL,
            CHIM_SLYR_B_TAX_ID INT NULL,
            UNALIGNED_SEQ_ID INT NULL,
            ALIGNED_SEQ_ID INT NULL,
            pynast_aligned_seq_id int null,
            PRIMARY KEY(GG_ID),
            FOREIGN KEY(NCBI_TAX_ID) REFERENCES TAXONOMY(TAX_ID),
            FOREIGN KEY(HUGENHOLTZ_TAX_ID) REFERENCES TAXONOMY(TAX_ID),
            FOREIGN KEY(SILVA_TAX_ID) REFERENCES TAXONOMY(TAX_ID),
            FOREIGN KEY(GREENGENES_TAX_ID) REFERENCES TAXONOMY(TAX_ID),
            FOREIGN KEY(UNALIGNED_SEQ_ID) REFERENCES SEQUENCE(SEQ_ID),
            FOREIGN KEY(ALIGNED_SEQ_ID) REFERENCES SEQUENCE(SEQ_ID),
            FOREIGN KEY(pynast_aligned_seq_id) REFERENCES SEQUENCE(SEQ_ID)
            )""")
        self.con.commit()

        self.cursor.execute("""
            CREATE TABLE GG_RELEASE(
            REL_ID INT NOT NULL AUTO_INCREMENT,
            NAME VARCHAR(20) NOT NULL,
            GG_ID INT NOT NULL,
            PRIMARY KEY(REL_ID),
            FOREIGN KEY(GG_ID) REFERENCES GREENGENES(GG_ID)
            )""")
        self.con.commit()

        self.cursor.execute("""CREATE TABLE otu_cluster (
            cluster_id INT NOT NULL,
            rep_id INT NOT NULL,
            similarity FLOAT NOT NULL,
            method VARCHAR(16),
            PRIMARY KEY(cluster_id),
            FOREIGN KEY(rep_id) REFERENCES greengenes(gg_id)
            )""")
        self.con.commit()

        self.cursor.execute("""CREATE TABLE otu (
            otu_id int NOT NULL AUTO_INCREMENT,
            cluster_id INT NOT NULL,
            gg_id INT NOT NULL,
            PRIMARY KEY(otu_id),
            FOREIGN KEY(cluster_id) REFERENCES otu_cluster(cluster_id),
            FOREIGN KEY(gg_id) REFERENCES greengenes(gg_id)
            )""")
        self.con.commit()

        self.cursor.execute("""INSERT INTO test_greengenes.sequence
                               SELECT s.seq_id,s.sequence
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.sequence s
                                ON g.pynast_aligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.sequence
                               SELECT s.seq_id,s.sequence
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.sequence s
                                ON g.unaligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.sequence
                               SELECT s.seq_id,s.sequence
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.sequence s
                                ON g.aligned_seq_id=s.seq_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.taxonomy t
                                ON g.ncbi_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.taxonomy t
                                ON g.hugenholtz_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.taxonomy t
                                ON g.greengenes_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.taxonomy
                               SELECT t.tax_id, t.tax_version, t.tax_string
                               FROM greengenes.greengenes g
                                INNER JOIN greengenes.taxonomy t
                                ON g.silva_tax_id=t.tax_id
                               WHERE g.gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.greengenes
                               SELECT *
                               FROM greengenes.greengenes
                               WHERE gg_id < 100""")
        self.cursor.execute("""INSERT INTO test_greengenes.gg_release
                               SELECT r.rel_id, r.name, r.gg_id
                               FROM test_greengenes.greengenes g
                                INNER JOIN greengenes.gg_release r
                                ON g.gg_id=r.gg_id""")
        self.con.commit()
