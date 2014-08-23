tupleMaker
==========

Creating ROOT ntuple from HEP data. Newest version is tupleMaker3.


**NOTE** : You have to compile it on 32bit clued0 machine, but you can run it on 64bit. Don`t know why :(


* setup for run & compilation:
'''bash
setup D0RunII p21.26.00 -O SRT_QUAL=maxopt
export LINK_SHARED=yes
'''
* the patch of pmcs for inputs created by tupleMaker3. It includes the saving tree instead of histograms.
  - highes head revision is 1.248 (Oh My God I Hate CVS Very Much)

* my fitting script for output of patched pmcs with tupleMaker3 is in directory `wfitter`
  - don't forget to add script to wzfitter/bin/BINARIES
