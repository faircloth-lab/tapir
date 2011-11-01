# test/ manifest:

* chr1_918.jukescantor.phydesign.rates - what we get from phydesign for JC
* chr1_918.subsmodel.phydesign.rates - what we get from the phydesign website for locus-specific subs model
* SiteRates_jukescantor_template.bf - phydesign hyphy script with minor modifications to add Jukes Cantor model

# generating data similar to phydesign

* parse the input tree and alignment (test-data/chr1_918.nex and test-data/Euteleost.tree) and modify tree using phydesign code

        perl phydesign.pl

* run hyphy against tree (here use hyphy1 since it yields results similar to phydesign)

        hyphy1 SiteRates_jukescantor_template.bf > test.jukescantor.hyphy_output

* convert those sites rates in phydesign fashion

        perl correction.pl

Note now that the differences are not so large and generally on the order of the results returned from phydesign.  The difference is largely due to the binary of hyphy that we use.