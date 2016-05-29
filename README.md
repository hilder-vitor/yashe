A proof-of-concept implementation of the YASHE leveled homomorphic cryptosystems
=============================================================================================================================

This is a modification of Lepoint's implementation, which can be found on the following repository: https://github.com/tlepoint/homomorphic-simon.

It is done in C++ using the FLINT library (http://www.flintlib.org/). The C++ wrappers are used, so our implementation requires the use of FLINT versions >= 2.4. We only tested for FLINT compiles with GMP and not MPIR.


WARNING
-------

This academic implementation is NOT to be used, not to be considered secured
nor pretty code. However we publish it under license CeCILL as a way to
support code-sharing and to allow the community to verify easily both the
correctness and the efficiency of this homomorphic evaluation.

