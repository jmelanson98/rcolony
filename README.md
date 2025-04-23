# rcolony
The RColony package provides simple tools for the implementation of the Colony pedigree analysis software.

# small updates!
Forked from jonesor/rcolony (thank you thank you thank you for a great wrapper package!)
Here I have made some small updates (April 22, 2025) to integrate outputs with the most recent version of COLONY2 (download from https://www.zsl.org/about-zsl/resources/software/colony).

Modifications were made to the function build.colony.input to create a properly formatted .DAT file for a command line (non-GUI) serial run on MacOS Sonoma 14.0. These modifications included:

- Adjustments to the writing of paternal and maternal sibship exclusion tables (changed write.table to writeLines to avoid addition of unnecessary spaces at the end of each row (which confused the COLONY software) --> note that similar modifications may be necessary for other inputs where the number of columns per row are inconsistent in the input file--I only made modifications where necessary for my workflow.
- Addition of input queries for:
    - line 8 (after dioecious/monoecious): 0/1 for inbreeding
    - line 11 (after mating systems): 0/1 for clone/duplicate inference
    - line 12 (after clone inference): 0/1 for sibship size scaling
- Modification of known paternity/maternity subsection to write an exclusion threshold even when the number of known parentages is 0 (necessary for COLONY interpretation)


Note that because COLONY2 expects Intel architecture, command line runs on ARM-based Apple Silicon processors should be run via Rosetta.

