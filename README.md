# Bell’s theorem stands and falls with absolute space

I propose a reformulation of the Bohm-Bell Gedanken Experiment that is in the
spirit of Leibniz. Instead of assuming Newton’s absolute space as a container
of the experimental setup, I presuppose that space is composed of the relations
encountered in a Bell-type experiment. In this approach, outcomes of spin component
measurements can be explained causally, and Bell’s conclusion becomes
invalid. As an offshoot, the maximum degree to which Bell’s inequality is violated
is shown to be indicative of the number of spatial dimensions.

This program accompanies an article with the same title.
The program produces the data that in the article is treated as coming from an ideal
Bell-Bohm experiment. Both Alice and Bob have many more than two settings to choose
between. Settings are disguised as 'keys'. Shared hidden variable explain the actual
outcomes. The statistical correlations agree with those predicted by Quantum Mechanics.

This can be done because it is not assumed that Alice and Bob share a system of reference,
contrary to what Bell does. To compensate, Alice and Bob must each have many more than
two choices.

The program takes as input:
   number of spatial dimensions
   number of keys on Alice's apparatus
   number of keys on Bob's apparatus
   number of trials

For example:

    ./belldata 3 30 30 30000000

The program produces data that violate the Clauser-Horne-Shimony-Holt inequality if
the number of dimensions is >2, up to the super-quantum limit of 4 for very high
dimenions. For 3 dimensions, the limit is the Tsirelson's bound, the same as the 
limit according to QM.
