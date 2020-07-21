# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 20:57:31 2020

@author: 1
"""
import squiggle
import numpy as np
import pandas as pd
pd.plotting.register_matplotlib_converters()
import matplotlib.pyplot as plt
# %matplotlib inline
import os
import seaborn as sns


from Bio import SeqIO 
for sequence in SeqIO.parse('MN908947.fna', "fasta"):
  print(sequence.seq)
  print(len(sequence),'nucliotides')

#Loading Complementary DNA Sequence into an alignable file
from Bio.SeqRecord import SeqRecord
DNAsequence = SeqIO.read('MN908947.fna', "fasta")
'''
Since input sequence is FASTA (DNA), and Coronavirus is RNA type of virus, we need to:

Transcribe DNA to RNA (ATTAAAGGTT… => AUUAAAGGUU…)
Translate RNA to Amino acid sequence (AUUAAAGGUU… => IKGLYLPR*Q…)
In the current scenario, the .fna file starts with ATTAAAGGTT, then we call transcribe() so T (thymine) is replaced with U (uracil), so we get the RNA sequence which starts with AUUAAAGGUU


'''
DNA = DNAsequence.seq#Convert DNA into mRNA Sequence
mRNA = DNA.transcribe() #Transcribe a DNA sequence into RNA.
print(mRNA)
print('Size : ',len(mRNA))

'''
The difference between the DNA and the mRNA is just that the bases T (for Thymine) 
are replaced with U (for Uracil).

Next, we need to translate the mRNA sequence to amino-acid sequence using translate() 
method, we get something like IKGLYLPR*Q ( is so-called STOP codon, effectively is a 
separator for proteins).

'''
Amino_Acid = mRNA.translate(table=1, cds=False)
print('Amino Acid', Amino_Acid)
print("Length of Protein:",len(Amino_Acid))
print("Length of Original mRNA:",len(mRNA))

#import the codon table
from Bio.Data import CodonTable
print(CodonTable.unambiguous_rna_by_name['Standard'])
'''
Let’s now identify all the Proteins (chains of amino acids), basically separating 
at the stop codon, marked by *. Then let’s remove any sequence less than 20 amino 
acids long, as this is the smallest known functional protein (if curious).

'''
#Identify all the Proteins (chains of amino acids)
Proteins = Amino_Acid.split('*') # * is translated stop codon
df = pd.DataFrame(Proteins)
df.describe()
print('Total proteins:', len(df))
def conv(item):
    return len(item)
def to_str(item):
    return str(item)
df['sequence_str'] = df[0].apply(to_str)
df['length'] = df[0].apply(conv)
df.rename(columns={0: "sequence"}, inplace=True)
df.head()# Take only longer than 20
functional_proteins = df.loc[df['length'] >= 20]
print('Total functional proteins:', len(functional_proteins))
functional_proteins.describe()

'''
Available Tools in ProtParam:

count_amino_acids: Simply count the number times an amino acid is repeated in 
the protein sequence.
get_amino_acids_percent: The same as only returns the number in the percentage 
of the entire sequence.
molecular_weight: Calculates the molecular weight of a protein.
aromaticity: Calculates the aromaticity value of a protein according to Lobry 
& Gautier (1994, Nucleic Acids Res., 22, 3174-3180).
flexibility: Implementation of the flexibility method of Vihinen et al. 
(1994, Proteins, 19, 141-149).
isoelectric_point: This method uses the module IsoelectricPoint to calculate 
the pI of a protein.
secondary_structure_fraction: This method returns a list of the fraction of 
amino acids that tend to be in helix, turn, or sheet.
Amino acids in Helix: V, I, Y, F, W, L.
Amino acids in Turn: N, P, G, S.
Amino acids in Sheet: E, M, A, L.
The list contains 3 values: [Helix, Turn, Sheet].
'''
from __future__ import division
poi_list = []
MW_list=[]
from Bio.SeqUtils import ProtParam
for record in Proteins[:]:
    print("\n")
    
    X = ProtParam.ProteinAnalysis(str(record))
    POI = X.count_amino_acids()
    poi_list.append(POI)
    MW = X.molecular_weight()
    MW_list.append(MW)
    POI=set(POI.values())
    if len(POI)==1 and list(POI)[0]==0:
      continue
    else:
      print("Protein of Interest = ", POI)
      print("Amino acids percent =    ",str(X.get_amino_acids_percent()))
      
      print("Molecular weight = ", MW_list)
      print("Aromaticity = ", X.aromaticity())
      print("Flexibility = ", X.flexibility())
      print("Isoelectric point = ", X.isoelectric_point())
      print("Secondary structure fraction = ",   X.secondary_structure_fraction())

#Plot the results
MoW = pd.DataFrame(data = MW_list,columns = ["Molecular Weights"] )#plot POI
poi_list = poi_list[548]
plt.figure(figsize=(10,6));
plt.bar(poi_list.keys(), list(poi_list.values()), align='center')      

'''
Now let us compare the similarity among COVID-19/COV2, MERS, and SARS.

Load the DNA sequence file (FASTA) each of SARS, MERS, and COVID-19.
'''
#Comparing Human Coronavirus RNA
from Bio import pairwise2
SARS = next(SeqIO.parse("sars.fasta",'fasta'))
MERS = next(SeqIO.parse("mers.fasta", "fasta"))
COV2 = next(SeqIO.parse("cov2.fasta", "fasta"))

#Before comparing the similarity let us visualize the DNA each of COV2, SARS, and MERS respectively
## Alignments using pairwise2 alghoritm
SARS_COV = pairwise2.align.globalxx(SARS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('SARS/COV Similarity (%):', SARS_COV / len(SARS.seq) * 100)
MERS_COV = pairwise2.align.globalxx(MERS.seq, COV2.seq, one_alignment_only=True, score_only=True)
print('MERS/COV Similarity (%):', MERS_COV / len(MERS.seq) * 100)
MERS_SARS = pairwise2.align.globalxx(MERS.seq, SARS.seq, one_alignment_only=True, score_only=True)
print('MERS/SARS Similarity (%):', MERS_SARS / len(SARS.seq) * 100)


# Plot the data
X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
Y = [SARS_COV/ len(SARS.seq) * 100, MERS_COV/ len(MERS.seq)*100, MERS_SARS/len(SARS.seq)*100]
plt.title('Sequence identity (%)')
plt.bar(X,Y)

'''
Conclusion
 
So we have seen how we can interpret, analyze the COVID-19 DNA sequence data, and tried to get as many insights regarding the proteins that made it up. In our result, we got the number of Leucine(L) and Valines(V) high in this protein which indicates a good number of Alpha-Helices.

Later we compared the COVID-19 DNA with MERS and SARS. And we saw that COVID-19 is closely related to SARS.



'''