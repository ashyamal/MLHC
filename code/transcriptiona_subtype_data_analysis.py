#Extracting/formating AD transcriptional subtype data into gene sets for input into MAGMA (data from Dr. Yosuke Tanigawa)

import pandas as pd
import re
import mygene

df = pd.read_excel(r'Tx.def.xlsm')
df = df.astype({'#Tx_idx': 'int', 'BETA': 'float'})
tx_dict = {}
for i in range(len(df)):
    k = df['#Tx_idx'].iloc[i]
    v = df['gene_symbol'].iloc[i]
    #print(k)
    if k in tx_dict.keys():
        tx_dict[k] = tx_dict[k] + [v]
    else:
        tx_dict[k] = [v]

mg = mygene.MyGeneInfo()
df = pd.DataFrame(pd.read_excel("Tx.def.xlsm"))
new_df = df[['Tx', 'gene_symbol']].copy()
new_df = new_df.astype({'gene_symbol': 'str'})
def gene_symbol_to_id(x):
    try:
        gene_id = mg.query('symbol:' + x, species='human')['hits'][0]['entrezgene']
        return str(gene_id)
    except:
        return '0000'
#new_df = new_df.groupby(['Tx'])['gene_symbol'].apply(lambda x: ' '.join(gene_symbol_to_id(x)))
new_df['gene_symbol']= new_df['gene_symbol'].apply(gene_symbol_to_id)
#new_df.to_csv('Tx_subtypes_grouped.txt', index = True, sep = ' ', header = False)
new_df.to_csv('Tx_gene_ids.txt', index = False, sep = ' ')

df = pd.DataFrame(pd.read_csv("Tx_gene_ids.csv", delim_whitespace = True))
df = df.astype({'gene_symbol': 'str'})
new_df = df.groupby(['Tx'])['gene_symbol'].apply(lambda x: ' '.join(x))
new_df.to_csv('Tx_subtypes_grouped.txt', index = True, sep = ' ', header = False)
