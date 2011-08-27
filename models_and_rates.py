#!/usr/bin/env python

from subprocess import Popen, PIPE
import os.path
import sys
import json

tree = os.path.abspath("Euteleost.tree")
alignment = os.path.abspath(sys.argv[1])
hyphy = os.path.abspath("HYPHY")


pipe = Popen([hyphy, 'models_and_rates.bf'], stdin=PIPE, stdout=PIPE)
towrite = "\n".join([alignment, tree])
for x in pipe.communicate(towrite):
    print x
pipe.wait()
