#!/usr/bin/env python
from BiochemPy import Reactions, Compounds

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
CompoundsHelper.saveCompounds(Compounds_Dict)
Aliases_Dict = CompoundsHelper.loadMSAliases()
CompoundsHelper.saveAliases(Aliases_Dict)
Names_Dict = CompoundsHelper.loadNames()
CompoundsHelper.saveNames(Names_Dict)

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
for rxn in Reactions_Dict:
    ReactionsHelper.rebuildReaction(Reactions_Dict[rxn])
ReactionsHelper.saveReactions(Reactions_Dict)
Aliases_Dict = ReactionsHelper.loadMSAliases()
ReactionsHelper.saveAliases(Aliases_Dict)
Names_Dict = ReactionsHelper.loadNames()
ReactionsHelper.saveNames(Names_Dict)
ECs_Dict = ReactionsHelper.loadECs()
ReactionsHelper.saveECs(ECs_Dict)
