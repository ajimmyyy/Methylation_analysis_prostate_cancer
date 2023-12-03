from pygosemsim import download
from pygosemsim import graph
from pygosemsim import similarity
from pygosemsim import annotation
from pygosemsim import term_set
import networkx as nx
import functools
import pandas as pd

if __name__ == "__main__":
    fn_go = 'Data/110820004/Data-GO_term/comorbidity_group_go_term.csv'

    # download.obo("go-basic")
    # download.gaf("goa_human")
    # df_go = pd.read_csv(fn_go)

    G = graph.from_resource("go-basic")
    similarity.precalc_lower_bounds(G)

    print(nx.ancestors(G, "GO:0004022"))
    print(nx.descendants(G, "GO:0005515"))

    print(similarity.wang(G, "GO:0004022", "GO:0005515"))
    print(similarity.resnik(G, "GO:0004022", "GO:0005515"))

    # trpv1 = annot["Q8NER1"]["annotation"].keys()
    # trpa1 = annot["O75762"]["annotation"].keys()
    # sf = functools.partial(term_set.sim_func, G, similarity.lin)
    # print(term_set.sim_bma(trpv1, trpa1, sf))