#!/usr/bin/env python2

"""
Python 2.7 code for determining whether a journal article X is missing
    any relevant citations.
Please provide:
- either the PMID of a published journal article X or a list of PMIDs
    (e.g. references of X article not yet published)
    These should be entered as a list of strings.
- The desired depth k of the citation graph
- The desired number n of most-cited papers to return
- Your email for accessing Entrez

Limitations of using Entrez:
- Only returns references to articles indexed by pubmed
    No books or papers not indexed by pubmed.
- Not all articles in pubmed have citation lists searchable by Entrez:
    This only works for articles whose Related Information includes
    "References for this PMC Article."

Improvements to be made:
- Add argparse for command-line input
- Add option for pubmed "related articles" instead of references only.
"""
import argparse
from operator import itemgetter
from Bio import Entrez


parser = argparse.ArgumentParser(description='Returns Missing References.')
parser.add_argument('email', help='Enter your email address to access Entrez.')
parser.add_argument('PMIDs', nargs='+',
                    help='List: PMID of paper X, or PMIDs of its references')
parser.add_argument('-depth', type=int, default=2, dest='k',
                    help='add the desired citation graph depth to search; '
                         'default is 2.')
parser.add_argument('-toprefs', type=int, default=20, dest='n',
                    help='add the number of top-cited references to return; '
                         'default is 20.')

args = parser.parse_args()
Entrez.email = args.email
starting_pmids = args.PMIDs
k = args.k
n = args.n
no_refs = []
checked_refs = []

# Entrez.email = "sample@email.com"
# # starting_pmids = ["19304878"]       # initial PMID(s)
# # starting_pmids = ["18629091"]
# starting_pmids = ['11483584', '18628940', '11726920']
# k = 2                               # depth of citation graph
# n = 20                              # number of "top cited" papers to return


# This function is primarily from https://www.biostars.org/p/89478/
def get_citations(pmid, errorlist):
    """ Attempt to return the PMIDs of the papers that paper pmid cites.

    Input pmid, the pmid of a paper for which to fetch info from Pubmed,
    and errorlist, a list of pmids of papers for which this doesn't work.

    If the paper has an pmid does not have an Entrez-accessible reference list,
    update the error list with the paper's PMID.
    """
    try:
        cites_list = []
        handle = Entrez.efetch("pubmed", id=pmid, retmode="xml")
        pubmed_rec = Entrez.read(handle)
        medline_info = pubmed_rec["PubmedArticle"][0]['MedlineCitation']

        for ref in medline_info['CommentsCorrectionsList']:
            if ref.attributes['RefType'] == 'Cites':
                cites_list.append(str(ref['PMID']))

        return cites_list

    except KeyError:
        errorlist.append(pmid)
        return None


# Get the first level of citations: a graph of depth 1.
if len(starting_pmids) > 1:
    # For checking a list of references of an unindexed paper
    paperX_refs = starting_pmids
    current_layer = starting_pmids
else:
    # If the initial paper is indexed by Pubmed
    paperX_refs = get_citations(starting_pmids, no_refs)
    current_layer = paperX_refs

all_refs = {}
for pmid in paperX_refs:
    all_refs[pmid] = 1
current_depth = 1

# Repeat citation collection for the desired depth of the graph
while current_depth <= k:
    current_depth += 1
    new_layer = []
    # print("This is step")
    # print(current_depth)
    # print("The current set of PMIDs to be checked is")
    # print current_layer
    # print("The list of already-checked PMIDs is")
    # print checked_refs
    for pmid in current_layer:
        id_refs = get_citations(pmid, no_refs)
        checked_refs.append(pmid)
        if id_refs is None:
            continue
        else:
            for ref in id_refs:
                all_refs[ref] = all_refs.get(ref, 0) + 1
                if ref not in (checked_refs + current_layer + new_layer):
                    new_layer.append(ref)
    current_layer = new_layer

sorted_refs = sorted(all_refs.iteritems(), key=itemgetter(1), reverse=True)
top_cited_refs = sorted_refs[:n]
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
