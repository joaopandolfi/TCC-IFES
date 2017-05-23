import pandas as pd
import glob
import os

# genes that didn't integrate with PPI
out = {'CHRNA7', 'DAOA', 'DTNBP1', 'MUTED', 'NPAS3', 'OFCC1', 'PRODH', 'SLC18A1'}


file_allseeds  = 'SZ_KATO/seeds/seeds_schizophrenia.txt'
suffix_seeds   = '/seeds/seeds_schizophrenia.txt'
suffix_exclude = '/seeds/seeds_schizophrenia.txt_EXCLUDED.txt'

df_allseeds = pd.read_table(file_allseeds, header=None)
allseeds = set(df_allseeds[1])

LOO = sorted(glob.glob('LOO_SZ_KATO_*'))
CV  = sorted(glob.glob('CV_*SZ_KATO_*'))

#outppi = ((allseeds-used) - excluded)

for DIR in CV:
    df_used     = pd.read_table(DIR + suffix_seeds, header=None)
    df_excluded = pd.read_table(DIR + suffix_exclude, header=None)
    #
    used = set(df_used[1])
    excluded = set(df_excluded[1])
    #
    print DIR, '\t', excluded


print out
for DIR in LOO:
    df_usedseeds = pd.read_table(DIR+'/seeds/seeds_schizophrenia.txt', header=None)
    used = set(df_usedseeds[1])
    loo  = allseeds - used
    print
    if loo-out == set([]):
        print 'Gene que não integrou:'
        print DIR, loo
    else:
        print 'Gene integrou OK:'
        print DIR, '\t', len(loo-out), len(loo), len(used), len(allseeds)
        print loo - out



