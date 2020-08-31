#!/usr/bin/env python
import json

with open("Original_Biochemistry.json", "r") as read_file:
    data = json.load(read_file)

data['compounds']=[]
data['reactions']=[]
data['compound_aliases']={}
data['reaction_aliases']={}

data['id']="MSD_v1.0"
del(data['__VERSION__'])

#NB compartments and cues currently not being handled
with open("Base_Biochemistry.json", "w") as write_file:
    write_file.write(json.dumps(data, indent=4))

