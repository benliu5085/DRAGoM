import subprocess
import os
import sys
import pickle
import getopt
import time

FF_SINGLE_END = 1
FF_PAIR_END_2 = 2
FF_PAIR_END_1 = 3

def analyzeCIAGR(ss):
    ans = 0
    p = 0
    q = 1
    char = ""
    cntM = 0
    while q < len(ss):
        while ss[q] >= '0' and ss[q] <= '9':
            q += 1
        ans += int(ss[p:q])
        if ss[q] == 'M':
            cntM += int(ss[p:q])
        char += ss[q]
        p = q + 1
        q += 1

    if char == 'M':
        return(cntM, ans)
    if len(char) == 2:
        return(cntM, ans)
    if len(char) == 3 and char[1] == 'M':
        return(cntM, ans)
    return(-1, ans)


def predictSAM(sam, name_l):
    path_pred = {}
    # e = float(sam[sam.find('_')+1:sam.rfind('.')])
    with open(sam, 'r') as fin:
        for line in fin:
            if line[0] != '@':
                cont = line.split()

                read = name_l[int(cont[0])]
                (overlap, read_len) = analyzeCIAGR(cont[5])

                header = cont[2]
                path_id = header[ header.find('_')+1 : header.find(',') ]

                cc = header.split('_')
                RNA = cc[-1]
                e = float(cc[-2])
                a = int(cc[-4])
                b = int(cc[-3])

                if overlap >= 0.6*(b-a) or overlap >= 0.6*read_len:
                    if read not in path_pred:
                        path_pred[read] = e
                    if e < path_pred[read]:
                        path_pred[read] = e
                    # if e not in path_pred:
                    #     path_pred[e] = {}
                    # if RNA not in path_pred[e]:
                    #     path_pred[e][RNA] = set([])
                    # path_pred[e][RNA].add(read)
    return path_pred

def mapping(path):
    bs = 3000
    if os.path.isfile(path):
        cnt_path = 0
        with open(path,'r') as fin:
            for line in fin:
                if line[0] == '>':
                    cnt_path += 1

        if cnt_path > bs:
            cmd = CDHIT + " -M 0 -c 0.90 -T " + CPU + " -o 90." + path + " -i " + path
            subprocess.call(cmd, shell=True)
            cmd = "mkdir tmp/"
            subprocess.call(cmd, shell=True)

            cls = []
            this_cls = []
            with open("90." + path + ".clstr", 'r') as fin:
                for l in fin:
                    if l[0] == '>':
                        if len(this_cls) != 0:
                            cls.append(this_cls)
                        this_cls = []
                    else:
                        this_cls.append( getPathID( l.split()[2] ) )
            cls.append(this_cls)

            path2file = {}
            for test_cls in cls:
                for (i,it) in enumerate(test_cls):
                    path2file[it] = 'tmp/' + path[:path.rfind('.')] + "_set_" + str(int(i/bs)) + '.fa'

            with open(path,'r') as fin:
                for l in fin:
                    if l[0] == '>':
                        header = l
                    else:
                        pid = getPathID(header)
                        if pid in path2file:
                            fout = open(path2file[pid], 'a')
                            fout.write(header)
                            fout.write(l)
                            fout.close()

            # sam = path + '.sam'
            sam = path[:path.rfind('.')] + ".sam"
            for map_path in set(list(path2file.values())):
                cmd = BWA + " index " + map_path
                subprocess.call(cmd, shell=True)
                cmd = BWA + " mem -t 16 -a " + map_path + ' '+ reads + ' | ' + SAMTOOL + ' view -F 4 >> ' + sam
                subprocess.call(cmd, shell=True)

            cmd = "rm -rf tmp/"
            subprocess.call(cmd, shell=True)
        else:
            sam = path[:path.rfind('.')] + ".sam"
            cmd = BWA + " index " + path
            subprocess.call(cmd, shell=True)
            cmd = BWA + " mem -t 16 -a " + path + ' '+ reads + ' | ' + SAMTOOL + ' view -F 4 > ' + sam
            subprocess.call(cmd, shell=True)
        return sam
    else:
        return ""

def cutPath(tblout, fq, sg, anchor, bit_score):
    cm_l = {}
    with open(anchor, 'r') as fin:
        for l in fin:
            if l[0] != '#':
                cont = l.split()
                header = cont[0]
                RNA = cont[3]
                a = min(int(cont[7]), int(cont[8]))
                b = max(int(cont[7]), int(cont[8]))
                a = a - 1
                e = cont[15]
                s = float(cont[14])
                if s > bit_score:
                    if header not in cm_l:
                        cm_l[header] = []
                    cm_l[header].append((a,b,RNA,e))
    fout = open('cut.'+fq, 'a')
    seq = ""
    header = ""
    with open(sg, 'r') as fin:
        for l in fin:
            if l[0] == '>':
                header = l[1:l.find(' ')]
            else:
                seq = l.strip()
                pos = cm_l.get(header,[])
                for it in pos:
                    (a, b, RNA, e) = it
                    fout.write('>')
                    fout.write(header)
                    fout.write('_')
                    fout.write(str(a))
                    fout.write('_')
                    fout.write(str(b))
                    fout.write('_')
                    fout.write(e)
                    fout.write('_')
                    fout.write(RNA)
                    fout.write('\n')
                    fout.write(seq[a:b])
                    fout.write('\n')
    fout.close()

    try:
        cm_l = {}
        with open(tblout, 'r') as fin:
            for l in fin:
                if l[0] != '#':
                    cont = l.split()
                    header = cont[0]
                    RNA = cont[3]
                    a = min(int(cont[7]), int(cont[8]))
                    b = max(int(cont[7]), int(cont[8]))
                    a = a - 1
                    e = cont[15]
                    if header not in cm_l:
                        cm_l[header] = []
                    cm_l[header].append((a,b,RNA,e))

        fout = open('cut.'+fq, 'a')
        seq = ""
        header = ""
        with open(fq, 'r') as fin:
            for l in fin:
                if l[0] == '>':
                    header = l[1:l.find(' ')]
                else:
                    seq = l.strip()
                    pos = cm_l.get(header,[])
                    for it in pos:
                        (a, b, RNA, e) = it
                        fout.write('>')
                        fout.write(header)
                        fout.write('_')
                        fout.write(str(a))
                        fout.write('_')
                        fout.write(str(b))
                        fout.write('_')
                        fout.write(e)
                        fout.write('_')
                        fout.write(RNA)
                        fout.write('\n')
                        fout.write(seq[a:b])
                        fout.write('\n')
        fout.close()
    except:
        print(tblout + " is empty!")

def splitFile(filename):
    max_char_count = 500*1024*1024
    cnt_file = 1
    char_count = 0
    idx = 0
    name_s = set([])
    with open(filename, 'r') as fin:
        for l in fin:
            char_count += len(l)
            newName = filename + '_' + str(cnt_file)
            name_s.add(newName)
            fout = open(newName,'a')
            fout.write(l)
            fout.close()
            if idx % 2 == 1:
                if char_count > max_char_count:
                    cnt_file += 1
                    char_count = 0

            idx += 1
    return list(name_s)

def mergeFile(outname, name_l):
    fout = open(outname,'w')
    for f in name_l:
        with open(f,'r') as fin:
            for line in fin:
                fout.write(line)
        os.remove(f)
    fout.close()

def compress(seqs, CPU):
    cnt_file_before = sys.maxsize
    while os.path.getsize(seqs) > 1024*1024*1024:
        name_l = splitFile(seqs)

        # compress each file
        cnt_file = len(name_l)
        name_l.sort()
        for f in name_l:
            out =  f + ".cd"
            cmd = CDHIT + " -M 0 -c 0.99 -T " + CPU + " -o " + out + " -i " + f
            subprocess.call(cmd, shell=True)
            os.remove(f)
            os.remove(out+".clstr")
            os.rename(out, f)

        mergeFile(seqs, name_l)

        # stop split if no significant reduction
        if cnt_file_before - cnt_file < 3:
            break
        cnt_file_before = cnt_file

    # one last time of compress
    out =  seqs + ".cd"
    cmd = CDHIT + " -M 0 -c 0.99 -T " + CPU + " -o " + out + " -i " + seqs
    subprocess.call(cmd, shell=True)
    os.remove(seqs)
    os.remove(out+".clstr")
    os.rename(out, seqs)

def rename(read_state, READ1, READ2, READS):
    if read_state == FF_PAIR_END_1:
        read_prefix = READ1[0 : READ1[:READ1.rfind('.')].rfind('.')]
    elif read_state == FF_PAIR_END_2:
        read_prefix = READ1[0 : READ1[:READ1.rfind('.')].rfind('.')]
    elif read_state == FF_SINGLE_END:
        read_prefix = READS[0 : READS[:READS.rfind('.')].rfind('.')]

    reads = read_prefix + '.reads.fa'
    namemap = read_prefix + '.reads.namemap.txt'

    fout1 = open(reads, 'w')
    fout2 = open(namemap, 'w')
    if read_state == FF_PAIR_END_1:
        with open(READ1, 'r') as fin:
            idx = 0
            for l in fin:
                if idx % 4 == 0:
                    fout1.write('>')
                    fout1.write(str(idx/4))
                    fout1.write('\n')
                    fout2.write(l[1:])
                elif idx % 4 == 1:
                    fout1.write(l)
                idx += 1
    elif read_state == FF_PAIR_END_2:
        with open(READ1, 'r') as fin:
            idx = 0
            for l in fin:
                if idx % 4 == 0:
                    fout1.write('>')
                    fout1.write(str(idx/4))
                    fout1.write('\n')
                    fout2.write(l[1:])
                elif idx % 4 == 1:
                    fout1.write(l)
                idx += 1
        with open(READ2, 'r') as fin:
            for l in fin:
                if idx % 4 == 0:
                    fout1.write('>')
                    fout1.write(str(idx/4))
                    fout1.write('\n')
                    fout2.write(l[1:])
                elif idx % 4 == 1:
                    fout1.write(l)
                idx += 1
    elif read_state == FF_SINGLE_END:
        with open(READS, 'r') as fin:
            idx = 0
            for l in fin:
                if idx % 4 == 0:
                    fout1.write('>')
                    fout1.write(str(idx/4))
                    fout1.write('\n')
                    fout2.write(l[1:])
                elif idx % 4 == 1:
                    fout1.write(l)
                idx += 1
    fout1.close()
    fout2.close()

    return (reads, namemap)

def usage():
    print("USAGE: python runDragom.py [options]")
    print("\n\noptions:\n\n")
    print(" -1   <filename>    : fastq file with forward paired-end reads\n")
    print(" -2   <filename>    : fastq file with reverse paired-end reads\n")
    print("                      only use -1 <filename> indicates that input reads are interleaved paired-end reads\n")
    print(" -s   <filename>    : fastq file with single-end reads\n")
    print("                      if both pair-end or single-end input are provided, the single-end input will be ignored\n")
    print(" -p   <filename>    : parameter file\n")
    print(" -m   <int>         : [optional] maximum extension length for anchors [default: 100] \n")
    print(" -d                 : [optional] anchor masking flag [use this to disable anchor masking] \n")
    print(" -h                 : print help message\n")


""" find all dependency """
home_dir = os.path.join(os.getcwd(), sys.argv[0][0:-12])
CMPRESS  = home_dir + "lib/cmpress "
CMSEARCH = home_dir + "lib/cmsearch "
GRARNA   = home_dir + "bin/dragom.exe "
BWA      = home_dir + "lib/bwa "
CDHIT    = home_dir + "lib/cd-hit-est "
SAMTOOL  = home_dir + "lib/samtools "
SGA      = home_dir + "lib/sga "
SPADES   = home_dir + "lib/SPAdes/bin/spades.py "
""" ---------------------- """

""" load parameters """
try:
    opts, args = getopt.getopt(sys.argv[1:], '1:2:s:m:p:dhH', ["help"])
except getopt.GetoptError as err:
        print str(err)
        usage()
        exit()

READ1 = ""
READ2 = ""
READS = ""
CONFIG = ""
MAX_EXTEND_LENGTH = "100"
AGGRESSIVE = True

for o, a in opts:
    if o == "-1":
        READ1 = a
    elif o == "-2":
        READ2 = a
    elif o == "-s":
        READS = a
    elif o == "-p":
        CONFIG = a
    elif o == "-m":
        MAX_EXTEND_LENGTH = a
    elif o in ("-d"):
        AGGRESSIVE = False
    else:
        usage()
        exit()

read_state = 0
if   (READ1 != "") and (READ2 != ""):
    read_state = FF_PAIR_END_2
elif (READ1 == "") and (READ2 == "") and (READS != ""):
    read_state = FF_SINGLE_END
elif (READ1 != "") and (READ2 == ""):
    read_state = FF_PAIR_END_1
elif (READ1 == "") and (READ2 != ""):
    read_state = FF_PAIR_END_1
    READ1 = READ2
else:
    print("[ERROR]: wrong setup for input read!\n\n\n")
    usage()
    exit()

if CONFIG == "":
    print("[ERROR]: you must provide config file!\n\n\n")
    usage()
    exit()

""" ---------------------- """

CPU = ""
SGA_CORRECR_K = ""
SGA_CORRECR_Learn = ""
SGA_CORRECR_d = ""
SGA_OVERLAP_M = ""

SPADES_META = ""
SPADES_M = ""

cm_l = []
with open(CONFIG, 'r') as fin:
    for line in fin:
        if line[0] != "#":

            cont = line.strip().split()
            if len(cont) > 0:
                if cont[0] == "CPU":

                    CPU = cont[1]
                if cont[0] == "SGA_CORRECR_K":
                    SGA_CORRECR_K = cont[1]
                if cont[0] == "SGA_CORRECR_Learn":
                    SGA_CORRECR_Learn = "-learn"
                if cont[0] == "SGA_CORRECR_d":
                    SGA_CORRECR_d = cont[1]
                if cont[0] == "SGA_OVERLAP_M":
                    SGA_OVERLAP_M = cont[1]

                if cont[0] == "SPADES_META":
                    SPADES_META = "--meta"
                if cont[0] == "SPADES_M":
                    SPADES_M = cont[1]

                if cont[0] == "in_cm":
                    cm_l.append(cont[1])
""" ---------------------- """

fout_log = open("Dragom.run.log", 'w')
""" stage 1: build graph """

fout_log.write("========== Welcome to DRAGoM ==========\n\n\n")
fout_log.write("========== stage 1: building graph ==========\n")
fout_log.write("========== stage 1a: building sga graph ==========\n")

time_begin = time.time()
fout_log.write("========== stage 1a: building sga graph ==========\n")
# preprocess
if read_state == FF_PAIR_END_1:
    cmd = SGA + " preprocess -o reads.pp.fastq --pe-mode 2 " + READ1
elif read_state == FF_PAIR_END_2:
    cmd = SGA + " preprocess -o reads.pp.fastq --pe-mode 1 " + READ1 + " " + READ2
elif read_state == FF_SINGLE_END:
    cmd = SGA + " preprocess -o reads.pp.fastq --pe-mode 0 " + READS
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

# ec
cmd = SGA + " index --no-reverse -t " + CPU + " reads.pp.fastq"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

cmd = SGA + " correct -k " + SGA_CORRECR_K + " -d " + SGA_CORRECR_d + " " + SGA_CORRECR_Learn + " -t " + CPU + " reads.pp.fastq"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

# rmdup
cmd = SGA + " index -t " + CPU + " reads.pp.ec.fa"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

cmd = SGA + " rmdup -t " + CPU + " reads.pp.ec.fa"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

# overlap
cmd = SGA + " overlap -t " + CPU + " -m " + SGA_OVERLAP_M + " reads.pp.ec.rmdup.fa"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to build sga graph!\n\n")

fout_log.write("========== stage 1b: building SPAdes graph ==========\n")
time_begin = time.time()

if read_state == FF_PAIR_END_1:
    cmd = SPADES + " " + SPADES_META + " -t " + CPU + " -m " + SPADES_M + " --12 " + READ1 + " -o spades"
elif read_state == FF_PAIR_END_2:
    cmd = SPADES + " " + SPADES_META + " -t " + CPU + " -m " + SPADES_M + " -1 " + READ1 + " -2 " + READ2 + " -o spades"
elif read_state == FF_SINGLE_END:
    cmd = SPADES + " " + SPADES_META + " -t " + CPU + " -m " + SPADES_M + " -s " + READS + " -o spades"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to build SPAdes graph!\n\n")


fout_log.write("========== stage 1c: merging graph ==========\n")
time_begin = time.time()
subprocess.call("gzip -d reads.pp.ec.rmdup.asqg.gz", shell=True)

# getSG
cmd = GRARNA + " getSG reads.pp.ec.rmdup.asqg"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)
subprocess.call("mv reads.pp.ec.rmdup.rename.StringGraph.fq og.fq", shell=True)

# mapping1
cmd = BWA + " index spades/contigs.fasta"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

# mapping2
cmd = BWA + " mem -a -T 45 -t " + CPU + " spades/contigs.fasta og.fq > og.sam"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

# merge
cmd = GRARNA + " mergeSG og.fq og.sam spades/contigs.fasta"
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)
string_graph = "og.merged.fq"
time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to build Hybrid graph!\n\n")

fout_log.write("========== stage 2: finding anchors ==========\n")
time_begin = time.time()
(reads, namemap) = rename(read_state, READ1, READ2, READS)
anchors = ""
paths = []
for cm in cm_l:
    # compress CM #
    acc = cm[cm.rfind('/')+1 : cm.rfind('.')]
    subprocess.call(CMPRESS + ' ' + cm, shell=True)

    # find anchor #
    tblout = string_graph[:string_graph.rfind('.')+1] + acc + ".tblout"
    anchors = anchors + tblout + ' '
    bit_score = ""
    with open(cm, 'r') as fin:
        for line in fin:
            if line[0:2] == "GA":
                cont = line.split()
                bit_score = float(cont[1])/10
    cmd = CMSEARCH + " --rfam -Z 5000000 -T " + str(bit_score) + " --nohmmonly --cpu " + CPU + " --tblout " + tblout + ' ' + cm + ' ' + string_graph
    fout_log.write("Command line:\n    "+cmd+"\n")
    subprocess.call(cmd, shell=True)

time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to find anchors!\n\n")

fout_log.write("========== stage 3: extend anchors ==========\n")
time_begin = time.time()
cmd = GRARNA + " extendAnchor -m " + MAX_EXTEND_LENGTH + " -f " + CDHIT + " -t " + CPU
if AGGRESSIVE:
    cmd += " -d "
cmd += string_graph + ' ' + anchors
fout_log.write("Command line:\n    "+cmd+"\n")
subprocess.call(cmd, shell=True)

for cm in cm_l:
    path = string_graph[:string_graph.rfind('.')+1] + acc + ".path.fa"
    compress(path, CPU)

time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to extend anchors!\n\n")

fout_log.write("========== stage 4: predicting ==========\n")
time_begin = time.time()

readname_l = []
with open(namemap,'r') as fin:
    for line in fin:
        readname_l.append(line[:-1])

subprocess.call("rm -f summary.csv", shell=True)
with open("summary.csv",'a') as fout:
    fout.write("acc,predicted reads,name\n")

for cm in cm_l:
    acc = cm[cm.rfind('/')+1 : cm.rfind('.')]
    path = string_graph[:string_graph.rfind('.')+1] + acc + ".path.fa"
    anchor = string_graph[:string_graph.rfind('.')+1] + acc + ".tblout"
    tblout = string_graph[:string_graph.rfind('.')+1] + acc + ".path.tblout"
    cmd = CMSEARCH + " --rfam -Z 5000000 --cut_ga --nohmmonly --cpu " + CPU + " --tblout " + tblout + ' ' + cm + ' ' + path
    fout_log.write("Command line:\n    "+cmd+"\n")
    subprocess.call(cmd, shell=True)

    bit_score = 0
    rna_name = ""
    with open(cm, 'r') as fin:
        for line in fin:
            if line[0:2] == "GA":
                cont = line.split()
                bit_score = float(cont[1])
            if line[0:4] == "NAME":
                cont = line.split()
                rna_name = cont[1]

    cutPath(tblout, path, string_graph, anchor, bit_score)
    cpath = "cut." + path
    compress(cpath, CPU)

    sam = mapping(cpath)
    if sam != "":
        this_pred = predictSAM(sam, readname_l)

    with open(acc+".prediction.csv",'w') as fout:
        fout.write("Read_Header,E-value\n")
        for read in this_pred:
            fout.write(read)
            fout.write(',')
            fout.write(str(this_pred[read]))
            fout.write('\n')

    with open("summary.csv",'a') as fout:
        fout.write(acc)
        fout.write(',')
        fout.write(str(len(this_pred)))
        fout.write(',')
        fout.write(rna_name)
        fout.write('\n')

    cmd = "mv " + cpath + " " + acc + ".assembled_rna.fa"
    subprocess.call(cmd, shell=True)

time_end   = time.time()
fout_log.write("\n\n\nIt takes " + str(time_end - time_begin) + " s to predict!\n\n")
fout_log.write("\n\n\n==================================")
fout_log.write("\nThanks for using DRAGoM!\n")

# clean intermeidate files
cmd = "rm -f reads.pp* "
subprocess.call(cmd, shell=True)
cmd = "rm -rf spades/"
subprocess.call(cmd, shell=True)
cmd = "rm -f og.fq og.sam"
subprocess.call(cmd, shell=True)
cmd = "rm -f og.merge* "
subprocess.call(cmd, shell=True)
cmd = "rm -f cut.og.merge* "
subprocess.call(cmd, shell=True)
cmd = "rm -f " + reads + " " + namemap
subprocess.call(cmd, shell=True)

""" ---------------------- """

fout_log.close()
