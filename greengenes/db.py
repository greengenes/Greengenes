#!/usr/bin/env python

from MySQLdb import connect, ProgrammingError
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

ALPHA = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
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

    def _update_tax_field(self, field_name, taxmap, version):
        """Insert taxonomy records, update greengenes records"""
        self._lock([('greengenes', None, 'write'), 
                    ('taxonomy', None, 'write')])

        next_tax_id = self._get_max_taxid() + 1
        for gg_id, tax_string in taxmap.items():
            if tax_string is None:
                try:
                    self.cursor.execute("""
                        UPDATE greengenes
                        SET %s=NULL
                        WHERE gg_id=%d""" % (field_name, int(gg_id)))
                    self.con.commit()
                except ProgrammingError:
                    self.con.rollback()
                    self._unlock()
                    raise ValueError, "Unable to insert taxonomy for gg_id %s!" \
                            % gg_id

            else:
                try:
                    self.cursor.execute("""
                        INSERT INTO taxonomy(tax_id, tax_version, tax_string)
                        VALUES (%d, '%s', '%s')""" % (next_tax_id, version, 
                                                      tax_string))
                    self.cursor.execute("""
                        UPDATE greengenes
                        SET %s=%d
                        WHERE gg_id=%d
                        """ % (field_name, next_tax_id, int(gg_id)))
                    self.con.commit()
                    next_tax_id += 1
                except ProgrammingError:
                    self.con.rollback()
                    self._unlock()
                    raise ValueError, "Unable to insert taxonomy for gg_id %s!" \
                            % gg_id

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
            raise ValueError, "Failed query!"
        
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
            raise ValueError, "Failed query!"

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

    def _update_seq_field(self, seqs, field_name):
        """update greengenes record"""
        self._lock([('greengenes', None, 'write'), 
                    ('sequence', None, 'write')])

        next_seq_id = self._get_max_seqid() + 1
        for gg_id, seq in seqs.items():
            try:
                self.cursor.execute("""
                    INSERT INTO sequence(seq_id, sequence)
                    VALUES (%d, '%s')""" % (next_seq_id, seq))
                self.cursor.execute("""
                    UPDATE greengenes
                    SET %s=%d
                    WHERE gg_id=%s
                    """ % (field_name, next_seq_id, gg_id))
                self.con.commit()
                next_seq_id += 1
            except:
                self.con.rollback()
                self._unlock()
                raise ValueError, "Unable to insert sequence for gg_id %s!" \
                        % gg_id

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
            raise ValueError, "Already have locks!"

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
            raise ValueError, "Unable to lock %s!" % \
                       ', '.join([i[0] for i in to_lock])
        
        self._have_locks = True
        self._lock_aliases = lock_aliases

    def _unlock(self):
        """Unlock locked tables"""
        if self._have_locks:
            self.cursor.execute("UNLOCK TABLES")
            self._have_locks = False
            self._lock_aliases = {}

    def _get_max_taxid(self):
        """Returns the max observed taxonomy id"""
        if self._have_locks:
            if 'taxonomy' not in self._lock_aliases:
                raise ValueError, "We're locked and cannot access table 'taxonomy'"
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
            raise ValueError, "Unable to query"

        return self.cursor.fetchone()[0]

    def _get_max_seqid(self):
        """Returns the max observed sequence id"""
        if self._have_locks:
            if 'sequence' not in self._lock_aliases:
                raise ValueError, "We're locked and cannot access table 'sequence'"
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
            raise ValueError, "Unable to query"

        return self.cursor.fetchone()[0]

    def _create_tmp_sequence_table(self):
        """ """
        name = self._get_tmp_table_name("tmp_seq")
        
        try:
            self.cursor.execute("""CREATE TEMPORARY TABLE %s 
                                   LIKE SEQUENCE""" % name)
            self.con.commit()
        except ProgrammingError:
            self.con.rollback()
            raise ValueError, "Unable to create table %s" % name

        return name

    def _drop_tmp_sequence_table(self, name):
        """ """
        if not name.startswith('tmp_'):
            raise ValueError, "%s doesn't appear to be temporary!" % name

        try:
            self.cursor.execute("DROP TABLE %s" % name)
            self.con.commit()
        except ProgrammingError:
            self.con.rollback()
            raise ValueError, "Unable to drop table %s" % name

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
    
    def _create_debug_db(self):
        """Create a small test database"""
        self.cursor.execute('USE test_greengenes')

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
