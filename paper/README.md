# Slot protocol figure

LaTeX version of slot-protocol-figure_3.pdf, written for the protocol that
threshold-ml-dsa implements (src/protocol/sign_slot.h, src/circuit/wrk_phase2.h).

The notation of the original figure is kept. Added on top of it:

- superscript BDOZ on authenticated shares and superscript WRK on
  garbled-circuit wire encodings;
- the explicit BDOZ-to-WRK conversion of the fixed inputs of C_post before
  the challenge;
- step 6 as implemented: free XOR of the retained WRK encodings of u, soldering
  of the XOR-ed wires onto the pre-garbled input wires (one authenticated
  block per wire per garbler, in place of the original identity gates),
  local evaluation.

Where the code departs from the figure (e_w dealt directly, delta_w = w + R_w),
the code's version is written.

Build with a TeX distribution providing pdflatex, latexmk, amsmath, amssymb,
bm, geometry, enumitem and mdframed:

~~~sh
make -C threshold-ml-dsa/paper
~~~

Output: slot-protocol-figure_3.pdf. `make clean` removes auxiliary files;
`make distclean` also removes the PDF.
