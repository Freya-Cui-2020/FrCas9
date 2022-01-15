import json
from sys import argv

dir=argv[1]

with open(dir+'/result.json','r') as f:
     data=json.load(f)
array=data['Sequences'][0]['Crisprs']
cas=data['Sequences'][0]['Cas']
if array:
    level=data['Sequences'][0]['Crisprs'][0]['Evidence_Level']
    direction=data['Sequences'][0]['Crisprs'][0]['Potential_Orientation']
    if cas:
        print("find crisprs and cas")
    else:
        print("find crisprs and no cas,further candidates")
else:
    if cas:
        print("find cas,no cirsprs,further candidates")
    else:
        print("excluded")
