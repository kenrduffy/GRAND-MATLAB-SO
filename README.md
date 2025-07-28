# GRAND-MATLAB-SO
Supersedes GRAND-MATLAB

Guessing Random Additive Noise Decoding (GRAND) (TM)

Subject to license: "GRAND Codebase Non-Commercial Academic Research Use License 021722.pdf"

Non-parallelized MATLAB implementations of: GRAND (hard detection); basic ORBGRAND (soft detection); 1-line ORBGRAND (soft detection); SGRAND (soft detection).

Code produces blockwise soft-output in the form of a post-decoding estimate of the likelihood that the decoding is correct.

Code uses the even code property for GRAND, basic ORBGRAND and 1-line ORBGRAND.

Simulation is setup and run with GRAND_Code/driver_GRAND.m and GRAND_Code/driver_SO.m

Sample output is in RESULTS, and sample plots from those results can be made with MAKE_FIGS/driver_sample_figs.m and MAKE_FIGS/driver_SO_figs.m

Note that for an [n,k] code, where k information bits become n coded bits, GRAND algorithms accurately and efficiently decode codes where n-k is moderate. This MATLAB implementation is solely intended to be instructive and is not parallelised, even though highly-parallelised implementations are possible. As a result, obtaining the full performance of a code with n-k>20 may prove time-consuming with the present implementation. 

The following should be cited in association with results from this code.

GRAND K. R. Duffy, J. Li, and M. Medard, "Capacity-achieving guessing random additive noise decoding," IEEE Trans. Inf. Theory, vol. 65, no. 7, pp. 4023–4040, 2019.

SGRAND A. Solomon, K. R. Duffy and M. Medard, "Soft maximum likelihood decoding using GRAND", IEEE ICC, 2020.

ORBGRAND K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, pp. 4528-4542, 2022.

SOGRAND P. Yuan, M. Médard, K. Galligan & K. R. Duffy, "Soft-output (SO) GRAND and long, low rate codes to outperform 5 LDPC codes", IEEE Trans. Wireless Commun., 24(4), 3386-3399, 2025.

For further details on GRAND, see: https://www.granddecoder.mit.edu/
