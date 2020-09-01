#!/usr/bin/env python
from biokbase.workspace.client import Workspace
import os, sys, json
Token = os.environ['KB_AUTH_TOKEN']
Workspace_URL = 'https://appdev.kbase.us/services/ws'
WSClient = Workspace(url = Workspace_URL, token = Token)
print(WSClient.ver())

with open('MSD_v1.0_Biochem.json', "r") as read_file:
    data = json.load(read_file)

results = WSClient.save_objects({'workspace':'kbase',
                                 'objects':[{'type':'KBaseBiochem.Biochemistry',
                                             'data':data,
                                             'name':'default'}]})
