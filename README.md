Missing Citations project
=================
missingcitations.py takes the PMID of a published journal article X, or a list of PMIDs (e.g. references of an article X that hasn't yet been published) and returns the n most cited articles in its citation graph of depth k.

Limitations
-----------------
This search works only for articles in pubmed whose citation lists are searchable by Entrez.  (i.e. articles for which the "Related Information" include "References for this PMC article")

Does not include books or any papers not indexed by pubmed.
