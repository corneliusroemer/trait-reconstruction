import click
from Bio import Phylo
import pandas as pd
from click import echo

def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name:
            if clade.name in names:
                raise ValueError("Duplicate key: %s" % clade.name)
            names[clade.name] = clade
    return names

def get_terminal_names(clade):
    terminal_list = []
    for tip in clade.get_terminals():
        terminal_list.append(tip.name)
    return terminal_list

def get_lineage_counts(terminal_names,meta):
    lineages = meta.loc[meta.taxon.isin(terminal_names),'lineage'].value_counts().to_dict()
    return lineages

@click.command()
@click.option('-t', '--tree', type=click.File('r'), default = 'data/tree.nwk')
@click.option('-m', '--metadata',type=click.File('r'), default = 'data/pango_default.csv')
@click.option('-o', '--outfile',type=click.File('w'), default = 'reconstructed.tsv')
def main(tree,metadata,outfile):
    meta = pd.read_csv(metadata)[['taxon','lineage']]
    tree = Phylo.read(tree, 'newick')
    lookup_dict = lookup_by_names(tree)
    internals = []
    for internal in tree.get_nonterminals():
        internals.append(internal.name)
    count = 0
    for internal in internals:
        if count % 20 == 0:
            echo(get_lineage_counts(get_terminal_names(lookup_dict[internal]),meta))
        count += 1
    # echo(lookup_dict[internals[10]])
    # reconstruct_lineage(lookup_dict(internals[0]),meta)




    

if __name__ == "__main__":
    main()