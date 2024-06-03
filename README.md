# Decoding CSIDH: A Comprehensive Guide to Isogeny-Based Cryptography
I plan to provide an uncompiled version of my thesis once it is published.
For now, I have added the code from the appendix, so people can get started.

## Sage Code
In the file `code_appendix.sage` I have added the code to compute the example in Section 5.2 of my thesis.
Upon removing every occurrence of `.montgomery_model()` in the code, one will get the example in Section 5.3.0 of my thesis.

The code was written with the mindset of providing an idea of how the theory translates into practice and isn't optimised for performance.
One can freely tweek the parameters in the code, like the prime `p` that defines the finite field you are working over as well as the ideals that are applied.
