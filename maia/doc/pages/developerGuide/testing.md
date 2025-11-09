# Testing # {#dvTesting}

Different testing possibilities for @maia are described in the following.

## Regression testing (canary)
The regression test framework is based on the tool
[canary](https://git.rwth-aachen.de/aia/MAIA/canary). Therefore, predefined test
cases are performed and checked against existing reference results. 

The continuous integration pipeline of [GitLab](@ref dvGitlab) uses canary under
the hood for automatic testing.

During every serious development process of @maia, canary should be used to
ensure the functionality of the code as you make changes and
[notice dangers as soon as possible](https://en.wiktionary.org/wiki/canary_in_a_coal_mine).
You can find the user guide in the
[canary](https://git.rwth-aachen.de/aia/MAIA/canary) repository.

## Unit testing
Unit test are based on [doctest](https://github.com/doctest/doctest) and located
in @maia repository under `test/`.  
Refer to the doctest
[documentation](https://github.com/doctest/doctest/tree/master/doc/markdown) and
the already existing cases, e.g. `test/test_lb_latticedescriptor.cpp`, to
understand the (quite simple) usage of doctest within @maia.
