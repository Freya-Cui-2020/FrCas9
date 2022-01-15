from Bio.Blast import NCBIXML
from sys import argv
import operator

r1=argv[1]
r2=argv[2]
blast_result=argv[3]
cds=argv[4]

r1=r1.strip().split('\t')[0]
r2=r2.strip().split('\t')[0]

r1=int(float(r1))
r2=int(float(r2))
rc=[]
for record in NCBIXML.parse(open(blast_result)): 
  if record.alignments: 
#    print("I am here~~~")
#    dr=record.query
    dr_len=record.query_length
    for align in record.alignments:
      for hsp in align.hsps:
#        qs=hsp.query_start
#        qe=hsp.query_end
        ql=hsp.identities
        t1=int(dr_len/5)+1
#        t2=int(float(dr_len)/3*2)+1
        mismatches=hsp.identities-hsp.positives
        ss=hsp.sbjct_start
        se=hsp.sbjct_end
        sr=int(((align.title.split(':'))[1]).split('-')[0])
#        print(str(sr))
#        print("test",str(ql),str(mismatches),str(ss+sr),str(se+sr))
        if ql>=t1 and ql<dr_len and mismatches<=3:
#          print("found",str(ql),str(mismatches))
          if ss>se:
            rr=se
          else:
            rr=ss
          rc.append((sr+rr-150))
          rc.append((sr+rr+150))
  else:
    print("not find the tracrRNA hits in 500bp range")
    exit()
if rc=="":
  print("No suitable potantial tracrRNA range!!!!!")
  exit()
print("-----------------------------------------------------------")
print("Potantial tracrRNA position is HERE.")
print(rc)

if min(rc)<(r1+500):
  r1=min(rc)
if max(rc)>(r2-500):
  r2=max(rc)
#print(str(r1),str(r2))

cds_dic=[]
cds_input=open(cds,'r')
for line in cds_input:
  line=line.strip()
  if line[0]=="#":
    continue
  else:
    lines=line.split('\t')
#    print(str(lines[3:5]))
    cds_dic.append([int(lines[3]),int(lines[4])])
cds_input.close()

#print(str(r1),str(r2))
b1=[]
change=0
newcds=sorted(cds_dic,key=operator.itemgetter(0))
if r1<int(newcds[0][0]):
  r1=0
  change=1
for k in newcds:
#  print(type(k))
#  print(str(k))
  if r1> int(k[0]) and int(r1<k[1]):
#    print(str(k[0]))
    r1=k[0]
    change=1
  elif r1>int(k[1]):
#    print(str(k[1]))
    b1.append(k[1])
  if r2> int(k[0]) and r2<int(k[1]):
    #print(str(k[0]))
    r2=k[1]
    break
  elif r2<k[0]:
    r2=k[0]
    break
#print(b1)    
if change==0:
  r1=max(b1)

print("-----------------------------------------------------------")
print("FINAL RANGE IS HERE")
print("final_range\t{0}\t{1}".format(r1,r2))

    
        