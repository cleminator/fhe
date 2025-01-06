# Clemens' FHE library

This repo contains my toy implementation of fully homomorphic encryption schemes, specifically the CKKS scheme.
My goal is simple: I want to write a fully functional implementation of CKKS in Python, without relying on any external libraries.

Therefore, while I will be implementing some commonly used optimization techniques to understand them better, the focus is neither on performance nor on security.

**BEWARE: This is not in any way a secure implementation, it is merely a learning tool for me. For example, I am not even using a CSPRNG, so use with caution.**


## References

Here are some references that I am using during the development of this library.

### General CKKS

- Stanford Lecture Notes on crypto, number theory, etc. (https://crypto.stanford.edu/pbc/notes/)
- Original CKKS paper (https://eprint.iacr.org/2016/421.pdf)
- Openmined Blog Series on CKKS (https://blog.openmined.org/ckks-explained-part-1-simple-encoding-and-decoding/)
- LWE Toolkit (https://eprint.iacr.org/2013/293.pdf)
- Sampling from Gaussian distributions (https://citeseerx.ist.psu.edu/document?repid=rep1&type=pdf&doi=4a3b8263ba0f8f8976499e7780e64decf2608412)
- Notes on Lattices, HE, and CKKS (https://arxiv.org/pdf/2205.03511)

### RNS & NTT:
- RNS-Version of CKKS: https://pmc.ncbi.nlm.nih.gov/articles/PMC8048025/pdf/nihms-1043784.pdf
- Efficient NTT for RLWE: https://eprint.iacr.org/2016/504.pdf
- Improvements on RNS-CKKS: https://eprint.iacr.org/2020/1118.pdf
