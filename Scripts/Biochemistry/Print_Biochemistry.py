#!/usr/bin/env python
import os, sys
temp=list();
header=1;

sys.path.append('../../Libs/Python')
from BiochemPy import Reactions, Compounds

ReactionsHelper = Reactions()
Reactions_Dict = ReactionsHelper.loadReactions()
ReactionsHelper.saveReactions(Reactions_Dict)

CompoundsHelper = Compounds()
Compounds_Dict = CompoundsHelper.loadCompounds()
CompoundsHelper.saveCompounds(Compounds_Dict)

