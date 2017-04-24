README for mmdecoder
====================

mmdecdoer is a simple HMM designed to annotate the genetic origins of a mouse
straing genotyped with the GigaMUGA array from Neogen.

## Usage ##


## Notes ##

The mmdecoder HMM has 55 fully connected states, each corresponding to a
reference genotype. The transition probabilities are the same for each
state. The chance to leave a state is called X. The chance to stay in a
state should therefore be 1 - 55X, but since X is so small, I just use
1.0.

The emission probabilities for each state depend entirely on the
corresponding reference genome. If the genotype is 'A', the emissions
are expected to be 100% 'A'. However, there is a small chance of error,
which is given by parameter <Y>.
