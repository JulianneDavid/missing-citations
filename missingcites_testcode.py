"""
Test code for missingcitations.py: runs on the mock citation set:

X -> A, B, C, D, E
-----
A -> C, D, F, G
B -> C, F, H
C -> ?? (references not indexed)
D -> H, E
E -> H
-----
F -> H, J
G -> H, E
H -> J

"""

from operator import itemgetter

starting_pmids = ["X"]              # initial PMID(s)
k = 2                               # depth of citation graph
n = 5                               # number of "top cited" papers to return

current_depth = 0                   # a counter for the depth of the graph currently

# This function is primarily from https://www.biostars.org/p/89478/
def get_citations(pmid, errorlist):
    """
    Returns the PMIDs of the papers this paper cites
         -or-
    if it does not have an Entrez-accessible reference list, updates the error list with the paper's PMID.
    """
    if pmid == ['X']:
        cites_list = ["A", "B", "C", "D", "E"]
    elif pmid == 'A':
        cites_list = ["C", "D", "F", "G"]
    elif pmid == 'B':
        cites_list = ["C", "F", "H"]
    elif pmid == 'D' or pmid == "G":
        cites_list = ["H", "E"]
    elif pmid == 'E':
        cites_list = ["H"]
    elif pmid == 'F':
        cites_list = ["H", "J"]
    elif pmid == 'H':
        cites_list = ["J"]
    else:
        errorlist.append(pmid)
        cites_list = None

    return cites_list


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

while current_depth <= k:
    current_depth += 1
    new_layer = []
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
        missing_citations.append(pmid[0])


print("The papers referenced by paper X are:")
print(paperX_refs)

print("The PMIDs of the top cited articles and their citation values are:")
print(top_cited_refs)

print("Of those, the following are not cited by paper X:")
print(missing_citations)

print("The following do not have references accessible by Entrez:")
print(no_refs)
