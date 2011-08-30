# single file run
python models_and_rates.py test/test-data/chr1_918.nex test/test-data/Euteleost.tree 0 100 --output test/

# multiple file run (untested)
ls *.nex | xargs -n1 -P8 python models_and_rates.py
