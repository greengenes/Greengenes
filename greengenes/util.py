#!/usr/bin/env python

from gzip import open as gzopen
from datetime import datetime
import os

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2012, Greengenes"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Development"

class NoSequenceError(Exception):
    pass

def log_f(line):
    start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
    return "%s\t%s\n" % (start_time, line)

def greengenes_open(file_fp, permission='U'):
    """Read or write the contents of a file
    
    file_fp : file path
    permission : either 'U','r','w','a'
    
    NOTE: univeral line breaks are always used, so 'r' is automatically changed
    into 'U'
    """
    if permission not in ['U','r','w','a']:
        raise IOError, "Unknown permission: %s" % permission

    if file_fp.endswith('gz'):
        # gzip doesn't support Ub
        if permission == 'U':
            permission = 'r'
        return gzopen(file_fp, permission)
    else:
        if permission == 'r':
            permission = 'U'
        return open(file_fp, permission)

def generate_log_fp(output_dir, basefile_name='log', suffix='txt',
                    timestamp_pattern='%Y%m%d%H%M%S'):
    """Originally from the QIIME project under qiime.workflow"""
    timestamp = datetime.now().strftime(timestamp_pattern)
    filename = '%s_%s.%s' % (basefile_name,timestamp,suffix)
    return os.path.join(output_dir,filename)

class WorkflowLogger(object):
    """Originally from the QIIME project under qiime.workflow"""
    def __init__(self,log_fp=None,open_mode='w',script_name="Not specified"):
        if log_fp:
            self._f = open(log_fp,open_mode)
        else:
            self._f = None 
        start_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('Logging started at %s\n' % start_time)
        self.write('Script: %s\n' % script_name)
        self.write('Greengenes workflow version: %s\n\n' % __version__)

    def write(self,s):
        if self._f:
            self._f.write(s)
            # Flush here so users can see what step they're
            # on after each write, since some steps can take
            # a long time, and a relatively small amount of 
            # data is being written to the log files.
            self._f.flush()
        else:
            pass 
    
    def close(self):
        end_time = datetime.now().strftime('%H:%M:%S on %d %b %Y')
        self.write('\nLogging stopped at %s\n' % end_time)
        if self._f:
            self._f.close()
        else:
            pass

    def __del__(self):
        """Destructor"""
        self.close()
