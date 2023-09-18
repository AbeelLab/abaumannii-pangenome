#!/usr/bin/env python3
from sys import argv
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy
if __name__ == '__main__':
    obofile = argv[1] # obo file obtained from GO website
    gomapfile = argv[2] # tab-separated file mapping ORFs to GO terms
    studyfile = argv[3] # study input (ORFs)
    if len(argv) > 4:
        outfile = argv[4]
    else:
        outfile = 'goestudy.out'

    godag = GODag(obofile, optional_attrs={'relationship'}) # load obo file
    with open(gomapfile, 'r') as f: # load mapping orf -> go IDs
        gomap = {}
        for line in f:
            orf,go = line.strip().split('t')
            go = go.split(' ')
            gomap[orf] = go
            pop = set(gomap.keys())
    with open(studyfile, 'r') as f: # load study orf
        study = set([line.strip() for line in f])
    study = study.intersection(pop)
    
    # Create the GOE study object
    g = GOEnrichmentStudy(pop, gomap, godag, propagate_counts=False, alpha=0.01, \
                          methods='fdr_bh')
    results = g.run_study(study) # run study
    g.wr_tsv(outfile, results) # save the study results to a tab-separated file
