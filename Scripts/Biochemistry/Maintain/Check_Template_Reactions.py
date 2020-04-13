#!/usr/bin/env python
import os, sys, time
from csv import DictReader
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Template_Dir='../../Templates/'
Unbalanced_Reactions=dict()
for template in os.listdir(Template_Dir):
        if not os.path.isdir(os.path.join(Template_Dir, template)):
            continue
        Unbalanced_Reactions[template]=dict()
        with open(os.path.join(Template_Dir, template, 'Reactions.tsv')) as infile:
            for line in DictReader(infile, dialect='excel-tab'):
                if(line['id'] in Reactions_Dict and "OK" not in Reactions_Dict[line['id']]['status']):
                    Unbalanced_Reactions[template][line['id']]=1

time_str = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
file = open("Template_Update.txt","a")
file.write("====================================\n")
file.write(time_str+"\n")
for template in sorted(Unbalanced_Reactions.keys()):
    file.write(template+": "+str(len(Unbalanced_Reactions[template].keys()))+"\n")
file.close()
