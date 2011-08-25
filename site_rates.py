#!/usr/bin/env python

from subprocess import Popen, PIPE
import os.path
import sys
import json

tree = os.path.abspath("Euteleost.tree")
alignment = os.path.abspath(sys.argv[1])
hyphy = os.path.abspath("HYPHY")

with open("all_models.json", "r") as rfile:
    all_models = json.load(rfile)

model = all_models[os.path.split(alignment)[1].split(".")[0]]

def enc(st):
# surrounds a string with braces
    return '{' + str(st) + '}'

model_str = ""
for key in ["AC", "AT", "AG", "CG", "CT", "GT"]:
    if key == "AG": # reference rate
        value = 1
    else:
        value = model[key]
    model_str += enc(value)
model_str = enc(model_str)

pipe = Popen([hyphy, 'SiteRates_GTR_template.bf'], stdin=PIPE, stdout=PIPE)
#import pdb; pdb.set_trace()
towrite = "\n".join([alignment, tree, model_str])
print pipe.communicate(towrite)
#pipe.stdin.write(towrite)
#print pipe.stdout.read()
pipe.wait()
