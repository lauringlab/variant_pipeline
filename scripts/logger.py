import datetime
import resource
import sys
import traceback


class StdErrLogger():
    """Writes basic utilization data to stderr"""
    def __init__(self, verbose = False):
        self._verbose = verbose

    def log(self, message, verbose=None):
        """Logs message to std err with optional stats on peak utilization"""
        verbose = self._verbose if verbose is None else verbose
        if (verbose):
            usage = resource.getrusage(resource.RUSAGE_SELF)
            memory_used = usage.ru_maxrss/1024
            function_name = traceback.extract_stack()[-2:-1][0][2]
            message = "usertime(s)={0:.0f}|systime(s)={1:.0f}|"\
                "peak_memory_used(mb)={2}|{3}|{4}". \
                format(usage.ru_utime, usage.ru_stime, memory_used, 
                    function_name, message)
        sys.stderr.write("{0}|{1}\n".format(datetime.datetime.today(), message))

