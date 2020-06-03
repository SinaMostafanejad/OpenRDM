Important Notes     {#notes}
===============

[TOC]

Every feature in <b>OpenRDM</b> is designed to provide the user with a
convenient way to reach its goal. Not surprisingly, there is a cost to be
paid for such convenience. As such, before using <b>OpenRDM</b>, 
we invite the user to spend some time to take a glance at some important points
and technical details about the <b>OpenRDM</b> library to become familiar with
pros and cons of each existing feature.

@section openmp OpenMP

The OpenMP application programming interface (API) has armed us with an efficient tool
to keep the sequential and parallelized versions of the code within the same
source code file. Such unification strategy presents a choice to the user between the
parallel or sequential implementations of the code. Due to the presence of many
nested loops of various depths, all efforts have been directed towards minimizing the 
number of OpenMP pragma directives and enlarging the parallel regions while minimizing
the number of implied and explicit barriers as well as the usage of single and critical
constructions. This conservative approach will not only reduces the overhead pertinent
to the OpenMP constructs, it bestows readability and performance upon our implementation.
Nevertheless, many important factors are required to be considered before making any judgements
about the performance of the code. First, apart from the author of the code, compiler also
plays a crucial role to the optimization of the code. Different compilers with various levels
of optimizations create an extra hidden layer behind the complexity of the sequential
or parallel code.
