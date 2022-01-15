import os
from sys import argv
import operator


drfile=argv[1]
rpfile=argv[2]
trans_dic={"A":"T","T":"A","C":"G","G":"C"}
def make_trans(myseq):
    trans_seq="".join([trans_dic[x] for x in myseq[::-1]])
    return trans_seq

drinput=open(drfile,'r')
mydic=[]
for line in drinput:
    line=line.strip()
    if line.startswith(">"):
        start=(line[1:].split("_"))[1]
    else:
        mydic.append([start,line])
drinput.close()

mydrdic=sorted(mydic,key=operator.itemgetter(0))
myfirst=mydrdic[0][1]
mylast=mydrdic[-1][1]

sum_first=0
sum_last=0
for i in range(1,len(mydrdic)-1,1):
#    print(mydrdic[i][1])
    if myfirst==mydrdic[i][1]:
        sum_first+=1
    elif mylast==mydrdic[i][1]:
        sum_last+=1
#print(str(len(mydrdic)),str(sum_first),str(sum_last))

output=open(rpfile,'w')
if sum_first==(len(mydrdic)-2):
    print("Confirmed Repeat\t+\t{0}\t{1}".format(mydrdic[0][0],myfirst))
    output.write(">repeat_+_{0}\n{1}".format(mydrdic[0][0],myfirst))
elif sum_last==(len(mydrdic)-2):
    print("Confirmed Repeat\t-\t{0}\t{1}".format(mydrdic[-1][0],make_trans(mylast)))
    output.write(">repeat_-_{0}\n{1}".format(mydrdic[-1][0],make_trans(mylast)))
elif sum_first>sum_last:
    print("Doubtful Repeat\t+\t{0}\t{1}".format(mydrdic[0][0],myfirst))
    output.write(">repeat_+_{0}\n{1}".format(mydrdic[0][0],myfirst))
elif sum_first<sum_last:
    print("Doubtful Repeat\t-\t{0}\t{1}".format(mydrdic[-1][0],make_trans(mylast)))
    output.write(">repeat_-_{0}\n{1}".format(mydrdic[-1][0],make_trans(mylast)))
output.close()
