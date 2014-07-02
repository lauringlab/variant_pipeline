import argparse
import datetime
import errno
import multiprocessing
import sys
import os
import csv
import zlib

from ngs_plumbing.dna import PackedDNABytes

from logger import StdErrLogger
import trimmer
import demultiplexer

#    <-- 5 ---><---- 6 --><----- 7 ---->
#    R R R R R B B B B B B C C C C C C C

TAG_LENGTH_PCR_COUNTER = 5
TAG_LENGTH_SAMPLE_BARCODE = 6
TAG_LENGTH_KLENOW_FRAGMENT = 7
LEFT_PRIMER_SEQUENCE = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
RIGHT_PRIMER_SEQUENCE = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
PRIMER_TAIL_SEARCH_LENGTH = 8
TRIM_READS = True
DEMULTIPLEX = True


class AmplificationPrimer():
    def __init__(self, sequence, search_length, start_position):
        self.sequence = sequence.upper()
        self.search_excerpt = self.sequence[(-1 * search_length):]
        self.start_position = start_position
        self.primer_length = len(self.sequence)
    
    def length_to_trim(self, read_sequence):
        '''
        Returns number of bases to be trimmed from start to remove primer.
        We only look for the last few bases of the primer starting from end of barcode.
        '''
        if self._primer_found(read_sequence):
            return self.primer_length + 1
        else:
            return 0
    
    def _primer_found(self, read_sequence):
        return read_sequence.upper().find(
                self.search_excerpt, self.start_position, self.primer_length + 2) > 0

class BCodes():
    def __init__(self):
        self.rand_tag_len = TAG_LENGTH_PCR_COUNTER
        self.bcode_len = TAG_LENGTH_SAMPLE_BARCODE
        self.context_len = TAG_LENGTH_KLENOW_FRAGMENT
        self.name = []
        self.seq = []
        self.bcid = []
        self.lookup_seq = {}
    def get_bc_name(self,seq):
        return self.name[self.lookup_seq[seq]]
    def get_total_len(self):
        return self.rand_tag_len + self.bcode_len + self.context_len
    def get_bc_start(self):
        return self.rand_tag_len
    def get_bc_end(self):
        return self.rand_tag_len + self.bcode_len
    def read(self, barcode_file):
        with open(barcode_file,'rb') as bf:
            for r in csv.reader(bf, delimiter=','):
                name = r[0]
                if name in self.name:
                    raise  ValueError("barcode name [{0}] is not unique".format(name))
                self.name.append(name)
                self.lookup_seq[r[1]] = len(self.seq)
                self.seq.append(r[1])
                self.bcid.append(r[2])
class Dupes():
    def __init__(self, common_id,
            l_seq, l_qual, r_seq, r_qual):
        self.ids = zlib.compress(common_id)
        self.nominal = zlib.compress('%s %s %s %s' % (
            l_seq, l_qual, r_seq, r_qual))
    def merge(self, common_id,
            l_seq, l_qual, r_seq, r_qual):
        self.ids = zlib.compress('%s %s' % (zlib.decompress( self.ids ),
            common_id))
    def get_id_list(self):
        return zlib.decompress(self.ids).split(' ')
    def get_read_list(self):
        return zlib.decompress(self.nominal).split(' ')

###############################################################################

def make_key(flag,dna):
    '''
    This is lossy: N is converted to A.
    This is only used for hashing the fragmentary keys for
    each read
    '''
    return flag + PackedDNABytes(dna)

def get_int_flag(key):
    return int(key[0])
        
def condense(interval):
    LOGGER.log("condensing {0}".format(interval))
    global procs
    global bcodes
    global args
    bc_start = bcodes.get_bc_start()
    bc_end = bcodes.get_bc_end()
    bc_total_len = bcodes.get_total_len()
    res = {}
    for b in bcodes.seq:
        res[b] = {}
    start = int(interval.split(':')[0])
    end = int(interval.split(':')[1])
    with open(args.left,'r') as lf, open(args.right,'r') as rf:
        lf.seek(start)
        rf.seek(start)
        while lf.tell() < end and rf.tell() < end:
            # left-side read
            l_id = lf.readline().strip()
            l_seq = lf.readline().strip()
            l_plus = lf.readline().strip()
            l_qual = lf.readline().strip()
            assert(l_plus=='+')
            # right-side read
            r_id = rf.readline().strip()
            r_seq = rf.readline().strip()
            r_plus = rf.readline().strip()
            r_qual = rf.readline().strip()
            assert(r_plus=='+')
            
            for b in bcodes.seq:
                flag = 0
                if b==l_seq[bc_start:bc_end]:
                    flag = 1
                if b==r_seq[bc_start:bc_end]:
                    flag = flag ^ 2
                if flag:
                    assert(l_id.split(' ')[0]==r_id.split(' ')[0])
                    common_id = l_id.split(' ')[0][1:]
                    # only store if there's a match on one side or the other
                    t = make_key(str(flag),
                        l_seq[0:bc_total_len] + r_seq[0:bc_total_len])
                    if t in res[b]:
                        res[b][t].merge( common_id,
                            l_seq, l_qual, r_seq, r_qual )
                    else:
                        res[b][t] = Dupes( common_id,
                            l_seq, l_qual, r_seq, r_qual )
    return res

def search_primer(ps, primer):
    '''
    returns number of bases to be trimmed from start to remove primer
    we only look for the last few bases of the primer starting from end of barcode
    through full length of  
    '''
    primer_tail = primer[(-1 * PRIMER_TAIL_SEARCH_LENGTH):]
    start = TAG_LENGTH_PCR_COUNTER + TAG_LENGTH_SAMPLE_BARCODE
    end = primer.length() + 2
    if ps.upper().find(primer_tail, start, end) > 0:
        return primer.length()
    else:
        return 0

    
def print_results(results,left_out,right_out,summary_out):
    LOGGER.log("saving multiplexed results")
    left_primer = AmplificationPrimer(
            LEFT_PRIMER_SEQUENCE, PRIMER_TAIL_SEARCH_LENGTH, \
            TAG_LENGTH_PCR_COUNTER + TAG_LENGTH_SAMPLE_BARCODE)
    right_primer = AmplificationPrimer(
            RIGHT_PRIMER_SEQUENCE, PRIMER_TAIL_SEARCH_LENGTH, \
            TAG_LENGTH_PCR_COUNTER + TAG_LENGTH_SAMPLE_BARCODE)
    
    global bcodes
    summary_out.write('\t'.join(('pydmx_id','barcode','left_key','right_key',
        'num_dupes','barcode_flag',
        'left_seq','left_qual','right_seq','right_qual',
        'duplicate_ids\n')))
    n = -1
    for r in results:
        for b in r:
            for d in r[b]:
                n = n + 1
                f = get_int_flag(d)
                read_list = r[b][d].get_read_list()
                id_list = r[b][d].get_id_list()
                left_trim = TAG_LENGTH_PCR_COUNTER + TAG_LENGTH_SAMPLE_BARCODE
                right_trim = TAG_LENGTH_PCR_COUNTER + TAG_LENGTH_SAMPLE_BARCODE
                key_len = bcodes.get_total_len()
                keys = (read_list[0][:key_len],read_list[2][:key_len])
                if not 1 & f:
                    left_trim = left_trim + left_primer.length_to_trim(read_list[1])
                if not 2 & f:
                    right_trim = right_trim + right_primer.length_to_trim(read_list[3])
                read_list[0] = read_list[0][left_trim:]
                read_list[1] = read_list[1][left_trim:]
                read_list[2] = read_list[2][right_trim:]
                read_list[3] = read_list[3][right_trim:]
                shared_id = '@%s:DUPE_%d:ID_%d:FLAG_%d' % (b,len(id_list),n,f)
                left_out.write('%s 1\n%s\n+\n%s\n' % (shared_id, read_list[0], read_list[1]))
                right_out.write('%s 2\n%s\n+\n%s\n' % (shared_id, read_list[2], read_list[3]))
                summary_out.write('%d\t%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t[%s]\n' % (
                    n,b,keys[0],keys[1],len(id_list),f,
                    read_list[0],read_list[1],read_list[2],read_list[3],
                    ','.join(id_list)))


def get_breaks(left_file,right_file,processes):
    p = processes * 2
    assert(os.path.getsize(left_file)==os.path.getsize(right_file))
    total_bytes = os.path.getsize(left_file)
    chunk_bytes = int(total_bytes / p)
    #print('cpus %s total bytes = %d chunk bytes %d' % (processes,total_bytes,chunk_bytes))
    startpoints = [0]
    with open(left_file,'r') as lf, open(right_file,'r') as rf:
        start_byte = chunk_bytes
        while start_byte < total_bytes:
            lf.seek(start_byte)
            while lf.readline().strip() != '+': 
                pass
            lf.readline()
            start_byte = lf.tell()
            if start_byte < total_bytes:
                startpoints.append(start_byte)
            start_byte = start_byte + chunk_bytes
        #print('startpoints (len %d)' % len(startpoints))
        #print(startpoints)
        intervals = []
        for i in range(0,len(startpoints),2):
            #print('i %d' % i)
            s = startpoints[i]
            e = startpoints[i+1]
            if e == startpoints[-1]:
                e = total_bytes
            lf.seek(s)
            rf.seek(s)
            lline = lf.readline().split(' ')[0]
            rline = rf.readline().split(' ')[0]
            #print(lline,rline)
            assert(lline==rline)
            assert(s<e)
            intervals.append('%d:%d' % (s,e))
        return intervals


def trim_reads(left, right, output_dir):
    if not TRIM_READS:
        return (left, right)

    LOGGER.log("trimming input")
    trimmed_dir = os.path.join(output_dir,"trimmed")
    makepath(trimmed_dir)
    trimmed_left = os.path.join(trimmed_dir, os.path.basename(left))
    trimmed_right = os.path.join(trimmed_dir, os.path.basename(right))
    with open(left, 'r') as left_in, \
        open(right, 'r') as right_in, \
        open(trimmed_left, 'w') as left_out, \
        open(trimmed_right, 'w') as right_out:
        trimmer.trim(left_in, right_in, left_out, right_out)
    return (trimmed_left, trimmed_right)


def demultiplex(left_deduped_filename, right_deduped_filename, barcodes, demultiplexed_dir):
    if not DEMULTIPLEX:
        return
    
    LOGGER.log("saving demultiplexed reuslts")
    makepath(demultiplexed_dir)

    barcode_files = {}
    for barcode in barcodes.seq:
        sample_name = barcodes.get_bc_name(barcode) 
        left_barcode_filename = os.path.join(demultiplexed_dir, sample_name + "." + os.path.basename(left_deduped_filename))
        right_barcode_filename = os.path.join(demultiplexed_dir, sample_name + "." + os.path.basename(right_deduped_filename))
        barcode_files[barcode] = [open(left_barcode_filename, 'w'), open(right_barcode_filename, 'w')]

    with open(left_deduped_filename, 'r') as left_deduped, \
            open(right_deduped_filename, 'r') as right_deduped:
        demultiplexer.demultiplex(left_deduped, right_deduped, barcode_files)
    
    for barcode in barcodes.seq:
        for file in barcode_files[barcode]:
            file.close()

def makepath(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise


##############################################################################

LOGGER = StdErrLogger()
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument( '-l', '--left', required=True,
            help='The left side mates.')
    parser.add_argument( '-r', '--right', required=True,
            help='The right side mates.')
    parser.add_argument( '-b', '--barcodes', required=True,
            help='The barcodes file; csv with three '
                    'columns: NAME SEQUENCE NUMBER (no header).')
#    parser.add_argument( '-p', '--processors', type=int, default=2,
#            help='Number of processors; twice the number of processes '
#            'will actually be used; this is intentional. You should give '
#            'the number of available cpus as the value for "-p" '
#            'for maximum performance.')
    timestamp = datetime.datetime.today().strftime('%Y%m%d_%H%M')
    parser.add_argument( '-o', '--output_dir', default="out-" + timestamp,
            help='Output dir')
    
    args = parser.parse_args()
    LOGGER.log("{0} begins".format(os.path.basename(sys.argv[0])))
    output_dir = os.path.realpath(args.output_dir)
    LOGGER.log("will write results to [{0}]".format(output_dir))
    
    makepath(output_dir)

    # global
    bcodes = BCodes()
    bcodes.read(args.barcodes)

    left_filename, right_filename = trim_reads(args.left, args.right, output_dir)

    # global
    #procs = 2 * args.processors #cgates/pulintz: This was raising exceptions in bpipe execution
    procs = 1
    pool = multiprocessing.Pool(processes=procs)
    # map_async to allow ctrl-c
    results = pool.map_async(condense, get_breaks(
        left_filename,right_filename,procs)).get(9999999)
    multiplexed_dir = os.path.join(output_dir, 'multiplexed')
    makepath(multiplexed_dir)
    left_out_filename = os.path.join(multiplexed_dir, '1.fastq')
    right_out_filename = os.path.join(multiplexed_dir, '2.fastq')
    summary_out_filename = os.path.join(multiplexed_dir, 'summary.txt')
    with open(left_out_filename,'w') as left_out, \
            open(right_out_filename, 'w') as right_out, \
            open(summary_out_filename,'w') as summary_out:
        print_results(results, left_out, right_out, summary_out)

    demultiplex(left_out_filename, right_out_filename, bcodes, os.path.join(output_dir,"demultiplexed"))

    LOGGER.log("results saved to [{0}]".format(output_dir))
    LOGGER.log("done")




















