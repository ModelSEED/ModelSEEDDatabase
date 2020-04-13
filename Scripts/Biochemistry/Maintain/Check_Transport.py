#!/usr/bin/env python
from BiochemPy import Reactions

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()

Update_Reactions=0
for rxn in sorted(Reactions_Dict.keys()):

    if(Reactions_Dict[rxn]["status"] == "EMPTY"):
        continue

    old_stoichiometry=Reactions_Dict[rxn]["stoichiometry"]
    Rxn_Cpds_Array=ReactionsHelper.parseStoich(old_stoichiometry)

    is_transport = ReactionsHelper.isTransport(Rxn_Cpds_Array)

    if(is_transport != Reactions_Dict[rxn]["is_transport"]):
        print("Updating: ",rxn,is_transport)
        Reactions_Dict[rxn]["is_transport"] = is_transport
        Update_Reactions+=1

if(Update_Reactions>0):
    print("Saving adjusted is_transport flag for "+str(Update_Reactions)+" reactions")
    ReactionsHelper.saveReactions(Reactions_Dict)
