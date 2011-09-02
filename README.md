## single file run

    python models_and_rates.py test/test-data/chr1_918.nex test/test-data/Euteleost.tree --output trash --hyphy hyphy2 --epochs=10-20,20-30 --times=1,2,3,4,5

## multiple file run (untested)

    ls *.nex | xargs -n1 -P8 python models_and_rates.py
