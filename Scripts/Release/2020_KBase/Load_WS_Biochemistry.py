#!/usr/bin/env python
from biokbase.workspace.client import Workspace
import os, sys, json
Token = os.environ['KB_AUTH_TOKEN']
Workspace_URL = 'https://appdev.kbase.us/services/ws'
WSClient = Workspace(url = Workspace_URL, token = Token)
print(WSClient.ver())

if(len(sys.argv)<2 or os.path.isfile(sys.argv[1]) is False):
    print("Takes one argument, the path to and including biochemistry object json")
    sys.exit()

with open(sys.argv[1], "r") as read_file:
    data = json.load(read_file)

results = WSClient.save_objects({'workspace':'kbase',
                                 'objects':[{'type':'KBaseBiochem.Biochemistry',
                                             'data':data,
                                             'name':'default'}]})
