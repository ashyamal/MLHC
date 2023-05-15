#Extracting/formating SNP Loc data from Alzheimer's GWAS dataset downloaded from NHGRI-EBI GWAS Catalog

import pandas as pd
import re
import mygene


df = pd.DataFrame(pd.read_csv("efotraits_EFO_0006930-associations-2023-05-8.csv"))

def convert_to_float(s):
    a = int(s[0])
    b = int(s[7])
    return a * 10 ** (- b)

def remove_html_tags(text):
    clean = re.compile('<.*?>')
    return re.sub(clean, '', text)

new_df = df[['Variant and risk allele', 'Location', 'Location']].copy()
new_df.iloc[:, 0] = new_df.iloc[:, 0].apply(lambda s : s[:s.index('-')])
new_df.iloc[:, 1] = new_df.iloc[:, 1].apply(lambda s : s[: s.rfind(':')])
new_df.iloc[:, 2] = new_df.iloc[:, 2].apply(lambda s : s[s.rfind(':') + 1 :])
pvals = df['P-value']
new_df['P'] = pvals
new_df.set_axis(['SNP', 'CHR', 'BP', 'P'], axis = 1, inplace =  True)
new_df = new_df[new_df.BP != 'Mapping not available']
new_df['P']= new_df['P'].apply(lambda x: convert_to_float(x))



df1 = pd.DataFrame(pd.read_csv("gwas_catalog_v1.0-associations_e109_r2023-05-07.csv", usecols = ['SNPS', 'INITIAL SAMPLE SIZE']))
df2 = df1.iloc[:,[1, 0]]
def parse_sample_col(s):
    return sum([int(x) for x in re.findall(r'\b\d+\b', str(s).replace(',', ''))])
def split_snps(s):
    if '; ' in s:
        return s.split('; ')
    elif ', ' in s:
        return s.split(', ')
    elif ',' in s:
        return s.split(',')
    else:
        return s.split(' x ')
df2['INITIAL SAMPLE SIZE'] = df2['INITIAL SAMPLE SIZE'].apply(parse_sample_col)
df2['SNPS'] = df2['SNPS'].apply(split_snps)
df2 = df2.explode('SNPS')
df2 = df2.astype({'SNPS': 'str', 'INITIAL SAMPLE SIZE': 'int'})


new_df1 = df2.groupby(['SNPS']).sum()


n_vals = pd.Series([new_df1.loc[s, 'INITIAL SAMPLE SIZE'] if s in new_df1.index else 0 for s in new_df.SNP])
new_df = new_df.set_index(n_vals.index)
new_df['N'] = n_vals
new_df.to_csv('SNP_p_vals_with_N.txt', index = False, sep = ' ')


'''
Supplemental code: performing MAGMA on Alzheimer's GWAS data (trait = Alzheimer's disease) from the NHGRI catalog, 
isolating significant genes using the Bonferonni correction

num_genes = 2525
bf_correction = .05/num_genes


df_list = []
count = 0
with open('/Users/ashyamal/Downloads/eMAGMA/Hippocampus.genes.out') as f:
    for line in f:
        if count != 0:
            line = line.strip()
            columns = re.split('\s+', line, maxsplit=10)
            df_list.append(columns)
        count += 1
df = pd.DataFrame(df_list)
df.set_axis(['GENE', 'CHR', 'START', 'STOP', 'NSNPS', 'NPARAM', 'N', 'ZSTAT', 'P', 'PERMP', 'NPERM'], axis = 1, inplace =  True)
df = df.astype({'P': 'float'})
new_df = df[df['P'] <= bf_correction]
new_df.sort_values(by = ['P'], inplace = True)
print(new_df)
new_df.to_csv('/Users/ashyamal/Downloads/eMAGMA/Hippocampus_signif_genes.txt', index = False, sep = ' ')
'''
