from Bio.Blast import NCBIXML
from Bio import SeqIO
from sys import argv
import operator

blast_result=argv[1]
reference=argv[2]
result_upper=argv[3]
result_down=argv[4]
result_sum=argv[5]
result_medium=argv[6]

wu=open(result_upper,'w')
wd=open(result_down,'w')
ws=open(result_sum,'w')
wm=open(result_medium,'w')
refs = SeqIO.to_dict(SeqIO.parse(reference, "fasta"))

E_VALUE_THRESH = 10
query_list=[]
global count
for record in NCBIXML.parse(open(blast_result)): 
    if record.alignments: 
#        print("\n") 
#        print("query: %s" % record.query)
        count=[]
        n=0
        spacer=record.query
        count.append(spacer)
        length=record.query_length
        print(length)
	for align in record.alignments:
            for hsp in align.hsps:
              qs=hsp.query_start
              qe=hsp.query_end
              mismatches=hsp.identities-hsp.positives
              if hsp.expect < E_VALUE_THRESH and mismatches <= 3:
                  n+=1
                  subj=align.title.split(' ')[1]
                  ss=hsp.sbjct_start
                  se=hsp.sbjct_end
                  if se > ss:
                     direction="p"
                     fss=ss-qs-10
                     fss=0 if fss<0 else fss
                     fse=se+(length-qe)+10
                    # print("%s %s %s," % (direction,fss,fse))
                     if refs[subj][fss:fss+10].seq:
                         wu.write("{up10}\t{spacer}\n".format(up10=str(refs[subj][fss:fss+10].seq),spacer=spacer))
                     if refs[subj][fse-10:fse].seq:
                         wd.write("{do10}\t{spacer}\n".format(do10=str(refs[subj][fse-10:fse].seq),spacer=spacer))
                     if refs[subj][(ss-qs):(se+(length-qe))].seq:
                         wm.write("{me}\t{spacer}\n".format(me=str(refs[subj][(ss-qs):(se+(length-qe))].seq),spacer=spacer))
                  else:
                     direction="n" 
                     fss=ss+qs+10
                     fse=se-(length-qe)-10
                     fse=0 if fse<0 else fse
                     #print("%s %s %s," % (direction,fss,fse))
                     if refs[subj][fss-10:fss].seq:
                         wu.write("{up10}\t{spacer}\n".format(up10=str(refs[subj][fss-10:fss].seq.complement())[::-1],spacer=spacer))
                     if refs[subj][fse:fse+10].seq:
                         wd.write("{do10}\t{spacer}\n".format(do10=str(refs[subj][fse:fse+10].seq.complement())[::-1],spacer=spacer))
                     if refs[subj][(se-length+qe):(ss+qs)].seq:
                         wm.write("{me}\t{spacer}\n".format(me=str(refs[subj][(se-length+qe):(ss+qs)].seq.complement())[::-1],spacer=spacer))
        count.append(n)
        query_list.append(count)
print(query_list)
for x in sorted(query_list,key=operator.itemgetter(1), reverse=True):
    ws.write("{spacer}\t{hits_number}\n".format(spacer=str(x[0]),hits_number=x[1]))
ws.close()
wu.close()
wd.close()
