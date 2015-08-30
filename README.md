AMSimulation
=============================
This software package contains three subpackages useful to perform L1 track trigger simulation using the associative memory approach within the official CMSSW framework

This package builds on top of a stand alone version developed here:
See [wiki_jf](https://github.com/jiafulow/SLHCL1TrackTriggerSimulations/wiki)

The initial code did
- Configuration: for generating particle gun samples
- NTupleTools: for flattening EDM ROOT files
- AMSimulation: for simulation of associative-memory-based track finding

This code aims to:
- Use the same Configuration step and use the NtupleTools to flatten the EDM Stubs collection
- Create a collection of L1 TTTracks that can be input to L1 Physics objects

Please do not develop on the 'master' branch.
A Full set of instructions can be found at: [wiki](https://github.com/rpatelCERN/AMSimulation/wiki) 
This will contain the full set of information for generating samples, the location of samples already made, and the instructions to compile and make ntuples for studies of L1 Physics Objects
