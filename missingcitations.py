#!/usr/bin/env python2

"""
Python 2.7 code for determining whether a journal article X is missing
    any relevant citations.
Required input:
- email address for accessing Entrez
    ** Question: does this need to be sent for each Entrez access, or just once
    in the code?

- either the PMID of a published journal article X or a list of PMIDs
    (i.e. references of X article not yet published)
    These should be entered as a list of strings.

Optional input:
- run mode:
    'graph' mode searches the citation graph of paper X for papers with the
    highest indegrees, and returns any of those top-cited papers that are not
    cited by paper X.

    'related' mode searches for papers in pubmed with high neighbor scores to
    the reference list of paper X.

    'intersections' and 'Jaccard' modes looks for similarities in the reference
    lists of paper X and other articles that cite the same papers, and returns
    papers included in similar lists but not in paper X's

        'intersections' looks only at the size of the list overlap.

        'Jaccard' calculates similarity by intersection/union to downweight
            papers with long reference lists.

- for 'graph' mode, the desired depth of the citation graph
- for 'graph' mode, the desired number of most-cited papers to return
- for 'related' mode, the number of related papers per paper X citation
    to return (total number will depend on the length of X's citation list)
- for 'intersections' and 'Jaccard' modes, the number of papers with similar
    citation graphs to return

Limitations of using Entrez:
- Only returns references to articles indexed by pubmed
    No books or papers not indexed by pubmed.
- Not all articles in pubmed have citation lists searchable by Entrez:
    This only works for articles whose Related Information includes
    "References for this PMC Article."

Improvements to be made:
- Add 'compare' mode??
"""
import argparse
from datetime import date
from operator import itemgetter
from Bio import Entrez


def get_citations(pmid, error_list):
    """ Attempt to return the PMIDs of the papers that paper pmid cites.

    Input pmid, the pmid of a paper for which to fetch info from Pubmed,
    and errorlist, a list of pmids of papers for which this doesn't work.

    If the paper has an pmid does not have an Entrez-accessible reference list,
    update the error list with the paper's PMID.

    This function is primarily from https://www.biostars.org/p/89478/
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
        error_list.append(pmid)
        return None


def get_related_papers(pmid, n, date):
    """ Returns n related articles to paper pmid with highest neighbor scores.

    Input one pmid, the number of most-related papers to return, and the date
    to check related papers' dates against.

    Notes on neighbor_score command:
        - This helpfully returns pmids already sorted by decreasing score.
        - There are several options that could be explored later: all results
        (which is what I'm using now); reviews; cited-in, which is not useful;
        "combined"; and others.
    """
    handle = Entrez.elink(dbfrom="pubmed", db="pubmed", cmd="neighbor_score",
                          id=pmid)
    related = Entrez.read(handle)[0]["LinkSetDb"][0]["Link"]
    if date:
        related_articles = []
        for link in related:
            id_date = get_date(link["Id"])
            if id_date < date:
                related_articles.append(link["Id"])
            if len(related_articles) == n:
                return related_articles
    else:
        for link in related[:n]:
            related_articles = link["Id"]
        return related_articles


def get_date(pmid):
    """ Get the publication date of paper pmid """
    handle = Entrez.efetch("pubmed", id=pmid, retmode="xml")
    pubmed_data = Entrez.read(handle)["PubmedArticle"][0]['PubmedData']

    for ref in pubmed_data["History"]:
        if ref.attributes['PubStatus'] == 'pubmed':
            year = (int(ref['Year']))
            month = (int(ref['Month']))
            day = (int(ref['Day']))

    pubmed_date = date(year, month, day)
    return pubmed_date


def get_citing_papers(pmid):
    """ Attempt to return the PMIDs of papers that cite paper pmid.

    This function is from the Biopython tutorial, section 9.15.3:
    http://biopython.org/DIST/docs/tutorial/Tutorial.html
    """
    results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc",
                                       LinkName="pubmed_pmc_refs", id=pmid))
    pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]

    results2 = Entrez.read(Entrez.elink(dbfrom="pmc", db="pubmed",
                                        LinkName="pmc_pubmed",
                                        id=",".join(pmc_ids)))
    pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
    return pubmed_ids


def get_intersection(list1, list2):
    """ Returns the intersection of list1 and list2 """
    count = 0
    intersection_list = []
    for item in set(list1).intersection(set(list2)):
        intersection_list.append(item)
    return intersection_list


parser = argparse.ArgumentParser(description='Return Missing References.')
parser.add_argument('--email', '-e', required=True,
                    help='email address (required to access Entrez)')
parser.add_argument('--pmids', '-p', nargs='+', required=True,
                    help="Paper X's PMID or "
                         "space-separated list of its references' PMIDs")
parser.add_argument('--mode', '-m', choices=['graph', 'related', 'compare',
                                             'intersection', 'Jaccard'],
                    help="Choose what mode to run:  "
                         "'graph' (default) searches Paper X's citation graph "
                         " for papers with high in-degrees.  "
                         "'related' searches pubmed for related articles.  "
                         "'intersection' compares the citation list with "
                         " those of papers whose citation lists have large"
                         " intersections with paper X's.  "
                         "'Jaccard' is the same as intersection, but measures"
                         " similarity of citation lists by intersection/union"
                         " instead of intersection only.  "
                         "'compare' gives the top results return by multiple"
                         " methods.")
parser.add_argument('--depth', '-d', type=int, default=2,
                    help='for "graph" mode: enter the maximum distance from '
                         'query paper to probe for missing cites; default 2.')
parser.add_argument('--toprefs', '-t', type=int, default=20,
                    help='for "graph" mode: enter the number of top-cited '
                         'references to return; default is 20.')
parser.add_argument('--numrelated', '-n', type=int, default=3,
                    help='for "related" mode: enter the number of related '
                         'papers to return per X reference; default is 3.')
parser.add_argument('--overlaps', '-o', type=int, default=5,
                    help='for "intersection" and "Jaccard" modes: number of '
                         'papers with similar reference lists to paper X; '
                         'default is 5.')

args = parser.parse_args()
Entrez.email = args.email
starting_pmids = args.pmids
no_refs = []
checked_refs = []
missing_citations = []

if args.mode:
    mode = args.mode
else:
    mode = 'graph'
    print "No mode was selected.  Default 'graph' mode is running."

# Entrez.email = "sample@email.com"
# # starting_pmids = ["19304878"]       # initial PMID(s)
# # starting_pmids = ["18629091"]
# starting_pmids = ['11483584', '18628940', '11726920']

# Get the PMIDs of the reference list of paper X
if len(starting_pmids) > 1:
    # For checking a list of references of an unindexed paper
    paperX_refs = starting_pmids
    # Assuming Paper X isn't published, no need to check dates
    date_check = False
else:
    # If the initial paper is indexed by Pubmed
    paperX_refs = get_citations(starting_pmids, no_refs)
    date_check = True

print("The papers referenced by paper X are:")
print(paperX_refs)

if mode == 'graph':
    n = args.toprefs
    k = args.depth
    # k = 2                               # depth of citation graph
    current_layer = paperX_refs

    all_refs = {}
    for pmid in paperX_refs:
        all_refs[pmid] = 1
    current_depth = 1

    # Repeat citation collection for the desired depth of the graph
    while current_depth <= k:
        current_depth += 1
        new_layer = []
        for pmid in current_layer:
            citation_list = get_citations(pmid, no_refs)
            checked_refs.append(pmid)
            if citation_list is None:
                continue
            else:
                for ref in citation_list:
                    all_refs[ref] = all_refs.get(ref, 0) + 1
                    if ref not in (checked_refs + current_layer + new_layer):
                        new_layer.append(ref)
        current_layer = new_layer

    sorted_refs = sorted(all_refs.iteritems(), key=itemgetter(1), reverse=True)
    top_cited_refs = sorted_refs[:n]

    for pmid in top_cited_refs:
        if pmid[0] not in paperX_refs:
            missing_citations.append(pmid[0])

    print("The PMIDs of the top cited articles and their citation values are:")
    print(top_cited_refs)

    print("Of those, the following are not cited by paper X:")
    print(missing_citations)

elif mode == 'related':
    n = args.numrelated
    top_refs = []

    if date_check:
        paperX_date = get_date(starting_pmids)
    else:
        paperX_date = False

    for pmid in paperX_refs:
        related_papers = get_related_papers(pmid, n, paperX_date)
        for ref in related_papers:
            if ref not in top_refs:
                top_refs.append(ref)
            if ref not in paperX_refs:
                missing_citations.append(ref)

    print("The PMIDs of the related articles are:")
    print(top_refs)

    print("Of those, the following are not cited by paper X:")
    print(missing_citations)

elif mode == 'intersection' or mode == 'Jaccard':
    n = args.overlaps
    co_citing_papers = []
    cite_dict = {}
    intersections = {}
    matches = []
    jaccards = {}
    missing_citations = []

    for pmid in paperX_refs:
        citing_papers = get_citing_papers(pmid)
        for paper in citing_papers:
            if paper not in co_citing_papers:
                co_citing_papers.append(paper)

    for paper in co_citing_papers:
        citation_list = get_citations(paper, no_refs)
        if citation_list is None:
            continue
        else:
            ref_overlap = get_intersection(citation_list, paperX_refs)
            intersections[paper] = len(ref_overlap)
            cite_dict[paper] = list(set(citation_list) - set(ref_overlap))
            if mode == jaccards:
                union = list(set().union(citation_list, ref_overlap))
                jaccards[paper] = len(ref_overlap) / len(union)

    if mode == intersection:
        items = intersections.iteritems()
    else:
        items = jaccards.iteritems()

    largest_intersections = sorted(items, key=itemgetter(1), reverse=True)[:n]
    for pmid in largest_intersections:
        matches.append(pmid[0])
        for ref in cite_dict[pmid[0]]:
            if ref not in missing_citations:
                missing_citations.append(ref)

    print "Top papers with similar citation lists to Paper X are:"
    print(matches)

    print("Papers cited by these that are not cited by paper X are:")
    print(missing_citationsj)
