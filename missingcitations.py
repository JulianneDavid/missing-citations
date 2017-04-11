"""
Python 2.7 code for determining whether a journal article X is missing any relevant citations
Please provide:
- either the PMID of a published journal article X or a list of PMIDs (e.g. references of X article not yet published)
    These should be entered as a list of strings.
- The desired depth k of the citation graph
- The desired number n of most-cited papers to return
- Your email for accessing Entrez

Limitations of using Entrez:
- Does not return any references to books or to papers not indexed by pubmed
- Not all articles in pubmed include citation lists searchable by Entrez:
    This only works for articles whose Related Information includes "References for this PMC Article."

Improvements to be made:
- Currently does not create or store an actual citation graph:
    Create citation graph!  association matrix?  ...other??
"""

from Bio import Entrez
from operator import itemgetter

Entrez.email = "please_enter@youremail.com"

starting_pmids = ["19304878"]       # initial PMID(s)
k = 3                               # depth of citation graph
n = 20                              # number of "top cited" papers to return

current_depth = 0                   # a counter for the depth of the graph currently


# This function is primarily from https://www.biostars.org/p/89478/
def get_citations(pmid, errorlist):
    """
    Returns the PMIDs of the papers this paper cites
         -or-
    if it does not have an Entrez-accessible reference list, updates the error list with the paper's PMID.
    """
    try:
        cites_list = []
        handle = Entrez.efetch("pubmed", id=pmid, retmode="xml")
        pubmed_rec = Entrez.read(handle)
#
        for ref in pubmed_rec["PubmedArticle"][0]['MedlineCitation']['CommentsCorrectionsList']:
            if ref.attributes['RefType'] == 'Cites':
                cites_list.append(str(ref['PMID']))

        return cites_list

    except KeyError:
        errorlist.append(pmid)
        return None


# Get the first level of citations: a graph of depth 1.
no_refs = []                        # for holding PMIDs for which Entrez doesn't return references

if len(starting_pmids) > 1:         # For checking a list of references of an unindexed paper
    paperX_refs = starting_pmids
    current_layer = starting_pmids
else:                               # If the initial paper is indexed by Pubmed
    current_depth += 1
    paperX_refs = get_citations(starting_pmids, no_refs)
    current_layer = paperX_refs

# The graph is now depth 1.  Add the first layer PMIDs to the dictionary and increase their values by 1.
all_refs = {}
for pmid in paperX_refs:
    all_refs[pmid] = 1

# Repeat citation collection for the desired depth of the graph
checked_refs = []

while current_depth < k:
    current_depth += 1
    new_layer = []
    # print("This is step")
    # print(current_depth)
    # print("The current set of PMIDs to be checked is")
    # print current_layer
    # print("The list of already-checked PMIDs is")
    # print checked_refs
    for pmid in current_layer:
        if pmid not in checked_refs:
            id_refs = get_citations(pmid, no_refs)
            checked_refs.append(pmid)
            if id_refs is None:
                continue
            else:
                for ref in id_refs:
                    # new_layer.append(ref)
                    all_refs[ref] = all_refs.get(ref, 0) + 1
                    if ref not in (checked_refs + current_layer + new_layer):
                        new_layer.append(ref)
    current_layer = new_layer


top_cited_refs = sorted(all_refs.iteritems(), key=itemgetter(1), reverse=True)[:n]
missing_citations = []

for pmid in top_cited_refs:
    if pmid[0] not in paperX_refs:
        missing_citations.append(pmid)


print("The papers referenced by paper X are:")
print(paperX_refs)

print("The PMIDs of the top cited articles and their citation values are:")
print(top_cited_refs)

print("Of those, the following are not cited by paper X:")
print(missing_citations)

print("The following do not have references accessible by Entrez:")
print(no_refs)
