#!/usr/bin/env python

from subprocess import Popen, PIPE
import os.path
import sys

tree = os.path.abspath("Euteleost.tree")
alignment = os.path.abspath(sys.argv[1])
hyphy = os.path.abspath("HYPHY")

pipe = Popen([hyphy, "scriptNucModelCompare.bf"], stdin=PIPE)
towrite = "\n".join([alignment, tree]) + "\n"
pipe.communicate(towrite)
pipe.wait()
