import argparse
import multiprocessing
import sys
import os
import csv
import zlib
from ngs_plumbing.dna import PackedDNABytes

class BCodes():
    '''
    <-- 5 ---><---- 6 --><----- 7 ---->
    R R R R R B B B B B B C C C C C C C
    '''
    def __init__(self):
        self.rand_tag_len = 5
        self.bcode_len = 6
        self.context_len = 7
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
                self.lookup_seq[r[1]] = len(self.seq)
                self.name.append(r[0])
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

def search_primer(ps):
    '''
    returns number of bases to be trimmed from start to remove primer
    full primer sequence is GTATTGCCAGTCACGACGATG, but we only look for
    the last 8 bases: "CGACGATG".  We look for these in bases [11:23] 
    '''
    if ps.upper().find('CGACGATG',11,23) > 0:
        return 22
    else:
        return 0
    
def print_results(results,left_out,right_out,summary_out):
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
                left_trim = 11
                right_trim = 11
                key_len = bcodes.get_total_len()
                keys = (read_list[0][:key_len],read_list[2][:key_len])
                if not 1 & f:
                    left_trim = left_trim + search_primer(read_list[0])
                if not 2 & f:
                    right_trim = right_trim + search_primer(read_list[3]) 
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




##############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument( '-o', '--output_prefix', default='OUT',
            help='Prefix for outputfiles.')
    parser.add_argument( '-l', '--left', required=True,
            help='The left side mates.')
    parser.add_argument( '-r', '--right', required=True,
            help='The right side mates.')
    parser.add_argument( '-b', '--barcodes', required=True,
            help='The barcodes file; csv with three '
                    'columns: NAME SEQUENCE NUMBER (no header).')
    parser.add_argument( '-p', '--processors', type=int, default=2,
            help='Number of processors; twice the number of processes '
            'will actually be used; this is intentional. You should give '
            'the number of available cpus as the value for "-p" '
            'for maximum performance.')
    args = parser.parse_args()

    # global
    bcodes = BCodes()
    bcodes.read(args.barcodes)
    # global
    procs = 2 * args.processors
    pool = multiprocessing.Pool(processes=procs)
    # map_async to allow ctrl-c
    results = pool.map_async(condense, get_breaks(
        args.left,args.right,procs)).get(9999999)
    with open('%s.1.fastq' % args.output_prefix,'w') as left_out, \
            open('%s.2.fastq' % args.output_prefix,'w') as right_out, \
            open('%s.summary.txt' % args.output_prefix,'w') as summary_out:
        print_results(results, left_out, right_out, summary_out)




















