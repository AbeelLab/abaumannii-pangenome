#!/usr/bin/env python3
import os, re, sys
import pandas as pd
import networkx as nx
if __name__ == '__main__':

	panaroodir = sys.argv[1]
	ptolemydir = sys.argv[2]

	# Load and process panaroo outputs
	panaroomap = os.path.join(panaroodir, 'gene_data.csv')
	panaroograph = os.path.join(panaroodir, 'final_graph.gml')
	centroid2Loc = pd.read_csv(panaroomap, index_col=2, header=0, \
        	                    names=['strain','location','clusterID','annot','protseq',\
                                'dnaseq','gene','desc']).drop(columns=['desc'])
	centroid2Loc.update(centroid2Loc.strain.apply(lambda x: x.replace('_reformatted','')))
	g = nx.read_gml(panaroograph, label='id')
	g = nx.relabel_nodes(g,int)
 	panaroodf = pd.DataFrame([v for k,v in g.nodes.items()], index=g.nodes())
	panaroodf.loc[:,'locustag'] = panaroodf.apply(lambda x: \
                                            centroid2Loc.loc[x.seqIDs].annot \
                                            if isinstance(x.seqIDs,str) \
                                            else centroid2Loc.loc[x.seqIDs].annot.values, axis=1)

	# Load and process ptolemy outputs
 	ptolemymap = os.path.join(ptolemydir, 'orf2id_mapping.txt')
 	safile = os.path.join(ptolemydir, 'syntenic_anchors.txt')

 	id2orf = pd.read_table(ptolemymap, delimiter='\t', header=None, index_col=2, \
			       names=['orf', 'strain', 'id'])
 	id2orf['locustag'] = id2orf.apply(lambda x: re.sub('^'+x.strain+'_', '', x.orf), axis=1)
 	id2orf['locustag'].update(id2orf.locustag.apply(lambda x: '_'.join(x.split('_')[1:-2]) \
							if len(x.split('_'))<8 else \
							'_'.join(x.split('_')[2:-2])))


 	# Match the panaroo outputs to ptolemys ORF indices
 	keeplocus = set(sum(panaroodf.locustag.apply(lambda x: [x] if isinstance(x,str) \
						     else list(x)),[]))
 	keeplocus = set(id2orf.locustag.values).intersection(keeplocus)
 	df = id2orf[id2orf.apply(lambda x: x.locustag in keeplocus, axis=1)].dropna()
 	df['newid'] = df.index # df with common ORFs, renumbered (0-based)
 	df.to_csv(ptolemymap, sep='\t', columns=['orf','strain','newid'], header=False, index=False)
 	id2orf = df.set_index('locustag')

 	# Remove ORFs discarded by panaroo from all ptolemy index files in the database
 	keepID = df.index
 	df = pd.read_table(os.path.join(ptolemydir,'id2fasta.txt'), delimiter='\t', header=None, \
			   index_col=0, names=['id','seq']).loc[keepIDs]
 	df.to_csv(os.path.join(ptolemydir,'id2fasta.txt'), sep='\t', header=False, index=True)
 	df = pd.read_table(os.path.join(ptolemydir,'global_z.txt'), delimiter='\t', header=None, \
			   index_col=0).apply(lambda x: set(map(int,x[1].split(','))), axis=1)
 	df = df.apply(lambda x: ','.join(map(str,x.intersection(keepID))))
 	df = df.drop(labels=df[df.apply(lambda x: not(x))].index)
 	df.to_csv(os.path.join(ptolemydir,'global_z.txt'), sep='\t', header=False, index=True)
 	df = pd.read_table(os.path.join(ptolemydir,'global_z_prime.txt'), delimiter='\t', header=None, \
			   index_col=0)
 	df.loc[keepID].dropna()[1].to_csv(os.path.join(ptolemydir,'global_z_prime.txt'), \
					  sep='\t', header=False, index=True)

 	# Writing out the syntenic anchor file
 	df = panaroodf.apply(lambda x: id2orf.loc[set(x.locustag).intersection(keeplocus),'newid'].values,axis=1)
 	anchors = df.apply(lambda x: ''.join(['\t'.join([str(mem), \
			','.join(map(str,set(x).difference([mem])))+'\n']) \
			for mem in x if len(x)>1]))
 	anchors[anchors.apply(lambda x: len(x)==0)] = np.nan
 	anchors.dropna().to_csv(safile, sep='\t', header=False, index=False, na_rep=None)
