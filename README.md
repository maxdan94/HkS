#README#

##Info##

bb_dks.c finds the heaviest subgraph of k nodes in a weighted graph (that is a subgraph containing exactly k nodes such that the sum of the weight on its edges is maximized). The result is exact. It scales to real-world networks containing several billions of edges and for k up to 10, 20 or more depending on the structure of the graph.

bb_dks_approx.c solves an approximated version of the problem. It finds a subgraph of k nodes or less such that the sum of the weight on its edges is at least 1/alpha times the one of the heaviest subgraph.

The algorithms are both based on branch and bound, the branching phase is based on deciding whether to include an edge in the subgraph or not, edges are examined in non-increasing order of weight.

wcore.c computes an ordering of the vertices according to the weighted kcore (or weighted degeneracy ordering) and compute the weight of the subgraphs induced by all prefixes. It can be used to find a heavy subgraph of any k, it does not have any fixed parameter approximation guaranty but is faster than the two algorithms based on branch and bound.

A stackexchange question on the topic: http://cstheory.stackexchange.com/questions/20221/find-the-densest-subgraph-of-size-k

A paper on the subject was published at: http://damnet.reading.ac.uk/

##To compile##

gcc bb_dks.c -o bb_dks -O3

gcc bb_dks_approx.c -o bb_dks_approx -O3

gcc wcore.c -o wcore -O3

##To execute##

./bb_dks k net.txt

./bb_dks_approx k alpha net.txt

- k is an integer
- alpha should be greater than one.
- net.txt should contain the weighted graph "nodeID1 nodeID2 weight". It is better if the nodeIDs are from 0 to n-1.

It will print:

- some info
- total_weight
- the k nodes

./wcore net.txt res.txt

Will print in res.txt "size nodeID weight" on each line, where weight is the sum of the weights of the subgraph induced on all previous nodes.

##Initial Contributors##

Maximilien Danisch  
Adapted from the java version of Manthos Letsios  
Technical consultants: Oana Balalau, Emmanuel Orsini and Mauro Sozio  
September 2016  
http://bit.ly/maxdan94  
maximilien.danisch@telecom-paristech.fr