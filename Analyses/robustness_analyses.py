# -*- coding: utf-8 -*-
import pandas as pd
from scipy.stats import spearmanr

import glob
from selection import *

#------------------------------------------------------------------------------
#In [112]: sum(abs(df_nodes['Y']) - abs(df_nodes['(D-C)/(C+D)']))
#Out[112]: -2.595978042447579e-10
#In [114]: sum(abs(df_nodes['X']) - abs(df_nodes['(C+D)']))
#Out[114]: 1.0001921478541842e-09
#------------------------------------------------------------------------------

def create_result_dataframe(n, IndexList, columns_order):
    df_result = pd.DataFrame({
        columns_order[0]: np.zeros(n, dtype='float'),
        columns_order[1]: np.zeros(n, dtype='int'),
        columns_order[2]: np.zeros(n, dtype='float'),
        columns_order[3]: np.zeros(n, dtype='float'),
        columns_order[4]: np.zeros(n, dtype='float'),
        columns_order[5]: np.zeros(n, dtype='float'),
        columns_order[6]: np.zeros(n, dtype='float')},
        index = IndexList, columns=columns_order)
    return df_result

#------------------------------------------------------------------------------

def compute_intersection_correlation(dirCV, partition_out,
                                    dfA, dfB, key, TOP, col_rank, columns_order,
                                    sort_order=[0,0], print_table=False):
    '''
    Compute correlation between DataFrames.

    '''
    Fields=[col_rank,key]
    topA = dfA.sort_values(by=Fields, ascending=sort_order)
    topB = dfB.sort_values(by=Fields, ascending=sort_order)
    #
    n = len(TOP)
    IndexList = [dirCV] * n
    df_result = create_result_dataframe(n, IndexList, columns_order)
    #
    for i,top in enumerate(TOP):
        A = set(topA[key][:top])
        B = set(topB[key][:top])
        unionAB = A | B
        interAB = A & B
        jaccardAB  = float(len(interAB)) / len(unionAB)
        proportionAB = float(len(interAB)) / top
        #
        tA = topA[:top]
        tB = topB[:top]
        keyA = tA[tA[key].isin(interAB)][key]
        keyB = tB[tB[key].isin(interAB)][key]
        #
        r = spearmanr(keyA,keyB)[0]
        #
        df_result.ix[i] = [
            partition_out,
            top,
            proportionAB,
            jaccardAB,
            r,
            r*proportionAB,
            r*jaccardAB]
        #
    if print_table == True:
        print df_result, '\n'
    return df_result


#------------------------------------------------------------------------------

def analyse_CV(CV, df_cmp, TOP, col_rank,
                    columns_order, key='GENE', print_table=False):
    #
    df_all_results = create_result_dataframe(0, [], columns_order)
    for dirCV in CV['in']:
        path_in  = dir_in  + dirCV
        path_out = dir_out + dirCV
        #
        #df_used     = pd.read_table(path_in + suffix_seeds, header=None)
        df_excluded = pd.read_table(path_in + suffix_exclude, header=None)
        #
        #used = set(df_used[1])
        excluded = set(df_excluded[1])
        #
        files = glob.glob(path_out + pattern_file_nodes)
        if len(files)==1:
            print('%s\t%s' % (dirCV, excluded))
            file_out = files[0]
            df_nodes = read_scores(file_out, new_header=['GENE','C','D'])
            #
            partition_out = int(dirCV.split('_')[1])
            #
            df_current = compute_intersection_correlation(dirCV, partition_out,
                            df_cmp, df_nodes, key, TOP, col_rank,
                            columns_order, print_table=print_table)
            #
            df_all_results = pd.concat([df_all_results, df_current])
        else:
            print('<<< ERROR >>>\t%s\t%s' % (path_in, excluded))
    #
    return df_all_results

#------------------------------------------------------------------------------


pd.set_option('display.height', 50) # height has been deprecated.
pd.set_option('display.max_rows', 25)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


#[int(c) if c.isdigit() else c for c in re.split('([0-9]+)', s)]

# genes that didn't integrate with PPI
outseeds = {'CHRNA7', 'DAOA', 'DTNBP1', 'MUTED', 'NPAS3', 'OFCC1', 'PRODH', 'SLC18A1'}

# Directories (input and output)
dir_in  = 'input/'
dir_out = 'output/'


#------------------------------------------------------------------------------

file_allseeds  = dir_in + 'SZ_KATO/seeds/seeds_schizophrenia.txt'
df_allseeds = pd.read_table(file_allseeds, header=None)
allseeds = set(df_allseeds[1])


suffix_seeds   = '/seeds/seeds_schizophrenia.txt'
suffix_exclude = '/seeds/seeds_schizophrenia.txt_EXCLUDED.txt'

#------------------------------------------------------------------------------
LOO = {}
LOO['in']  = sorted(glob.glob1(dir_in,  'LOO_SZ_KATO_*'))
LOO['out'] = sorted(glob.glob1(dir_out, 'LOO_SZ_KATO_*'))

CV = {}
CV['in']   = sorted(glob.glob1(dir_in,  'CV_*SZ_KATO_*'))
CV['out']  = sorted(glob.glob1(dir_out, 'CV_*SZ_KATO_*'))

#------------------------------------------------------------------------------

print 'Teste: diretorios input == outputs (subdiretorios):'
print set(CV['in']) == set(CV['out'])
print set(LOO['in']) == set(LOO['out'])
#------------------------------------------------------------------------------

pattern_file_nodes = '/out_CV_*-3_extracted_features-NODES.txt'
#TOP = [10,20,30,40,50,75,100,150,200,300,400,500]
TOP = [10,20,30,40,50,100,200]
key='GENE'

#------------------------------------------------------------------------------

file_cmp = 'output/COMPARATIONS/out_KATO_epsilon_05_20140901202044-3_extracted_features-NODES.txt'
df_cmp = read_scores(file_cmp, new_header=['GENE','C','D'])
file_out = 'output/CV_10_2_SZ_KATO_1/out_CV_10_2_SZ_KATO_1_05_20161124185819-3_extracted_features-NODES.txt'
df_nodes = read_scores(file_out, new_header=['GENE','C','D'])

#------------------------------------------------------------------------------

columns_order = [
                'OUT',
                'TOP',
                'OVERLAP',
                'JACCARD',
                'SPEARMAN',
                'SPEARMAN_OVERLAP',
                'SPEARMAN_JACCARD']

dfX = analyse_CV(CV, df_cmp, TOP, 'X',
               columns_order, key, print_table=True)

dfS = analyse_CV(CV, df_cmp, TOP, 'S',
               columns_order, key, print_table=True)

dfLooX = analyse_CV(LOO, df_cmp, TOP, 'X',
               columns_order, key, print_table=True)


print("######################################################################")
print(dfLooX)

plot_analyse(dfX, 'X')
plot_analyse(dfS, 'S')

######################################################
def plot_df_boxplot(df_top, var, top, s_overlap, s_out,
                    title, xlabel, ylabel):
    OUT = df_top.index.unique()
    dict_plot = {}
    for out in OUT:
        dict_plot[out] = df_top.ix[out][s_overlap].reset_index(drop=s_out)
    #
    df_plot = pd.DataFrame(dict_plot)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.axis([0.5, 4.5, 0.0, 1.0])
    df_plot.boxplot(grid=False)
    filename = 'fig_'+var+'_'+str(top)+'_'+str(out)+'.pdf'
    plt.savefig(filename, format='pdf')
    plt.close()



def plot_analyse(dfA, varA,
                 keep_columns=['TOP','OUT','OVERLAP'],
                 new_index=['TOP','OUT']):
    s_out     = keep_columns[1]
    s_overlap = keep_columns[2]
    #
    df_analyse = dfA[keep_columns].set_index(new_index)
    TOP = df_analyse.index.levels[0]
    #OUT = df_analyse.index.levels[1]
    #
    xlabel = u'(%) Sementes excluídas'
    ylabel = u'(%) Interseção com original'
    for top in TOP:
        title = u'(%d primeiros, escore %s)' % (top, varA)
        df_top = df_analyse.ix[top]
        #
        plot_df_boxplot(df_top, varA, top, s_overlap, s_out, title, xlabel, ylabel)


######################################################



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Centrality Measures
#------------------------------------------------------------------------------
dfCM = pd.read_table('GS_Centrality_Measures.txt')

#------------------------------------------------------------------------------


# FAZER LOO
