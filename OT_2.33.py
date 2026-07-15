# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 09:49:34 2020

@author: Laurent
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os, six
import dabest
import pingouin as pg
from pingouin import kruskal, pairwise_corr
#import datetime
#from matplotlib.backends.backend_pdf import PdfPages
#from scipy.optimize import curve_fit
import scikit_posthocs as sp
from scipy import stats
from scipy.stats.stats import pearsonr
#from scipy.stats import wilcoxon
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")
# --- PHP
import dataframe_image as dfi # dfi.export(df, 'dataframe.png')

# 20210315 todo: compare Omar/gautier cond by cond; sort by 'Speed [nm/sec]' for Omar
#   compare tbk1 vs tau or tdrop1 vs tau

plt.close(fig='all')
path0="/home/php/Bureau/Fabio-data/211120-Fabio/"
#path0="/home/laurent/DATA/Experiences Diverses/OT/Fabio/22November2021/"

#path0="/home/laurent/DATA/Experiences Diverses/OT/Fabio/28October2021/"

OoG="Gautier"; cd11a0rConA='cd11a'
#OoG="Omar"; cd11a0rConA='ConA'
OoG="Lucie"; cd11a0rConA='cd11a'
OoG="All"; cd11a0rConA='cd11a'

print('Intern :', OoG)
#ListIntern = ["Gautier", "Omar", "Lucie"]
ListIntern = ["Gautier", "Lucie"]

#ListIntern=['Omar']

lpathA = ['fitted_data_160921_free/', 'fitted_data_160921_E2est/']
lpathB = ['_160921_allFreeExceptEta', '_160921_E2eqE2_EST']
lconstraint = [False, True]
lpathA = ['']
lpathB = ['']
lconstraint = [False]

OneFile = True

if OneFile :
    dFint = pd.read_excel(path0+'fitted_data_all_curves_ALL.xlsx')
    dFint['Constraint'] = False
    if "Gautier" not in ListIntern: dFint = dFint[dFint['InternType']!=0]
    if "Omar" not in ListIntern: dFint = dFint[dFint['InternType']!=1]
    if "Lucie" not in ListIntern: dFint = dFint[dFint['InternType']!=2]
else:
    dFboth = []
    for (pathA, pathB, constraint) in zip(lpathA, lpathB, lconstraint):
        path1 = path0+pathA
        if OoG == 'All':
            dFl = []
            for name in ListIntern:
                path2='fitted_data_all_curves_ordered_'+name+pathB
                FFx = path1+path2+'.xlsx'
                dFx = pd.read_excel(FFx)
                dFx['Constraint'] = constraint
                print('load ', FFx)
                dFl.append(dFx)    
            dFint0 = pd.concat(dFl, ignore_index=True)
        else:
            path2='fitted_data_all_curves_ordered'+OoG+pathB
            FF = path1+path2+'.xlsx' #+'_180321.csv'
            dFint0 = pd.read_excel(FF)
            dFint0['Constraint'] = constraint
            print('load ', FF)
        dFboth.append(dFint0)
    dFint = pd.concat(dFboth, ignore_index=True)

dFint['best_slope_near_t0 [pN/sec]'] = dFint['best_slope_retract_pos [pN/sec]']

SaveGraph=True; OutFormat=".jpg";
#OutFold=path1+'Output20210916_'+OoG+'_free/'
#OutFold=path0+'Output20210916_'+OoG+'_test/'
#OutFold=path0+'Output20211025_'+OoG+'_test/'
OutFold=path0+'Output20211028_'+OoG+'_php-test-NEWREJECTFABIO-all/'
OutFold=path0+'Output20210509_'+OoG+'_test-all-CritStricts/'
OutFold=path0+'Output20210509_'+OoG+'_test-all-FabioTRUECleanSlope/'
#OutFold=path0+'Output20210509_'+OoG+'_test-all-FabioFALSE/'
#OutFold=path0+'Output20210509_'+OoG+'_test-all-FabioTRUE/'
OutFold=path0+'Output20211122_'+OoG+'_test-all-ForceMin6FabioTRUECleanSlope10_RatiosPhysE2cut/'

#OutFold=path1+'Output20210916_'+OoG+'/'
#OutFold=path1+'Output20210916_'+OoG+'_E2est/'

if not os.path.exists(OutFold): os.makedirs(OutFold)
fused_CD3=True  # fusionb of cd3 and cd3_ucht1
printcolumns=True
countpopulations=False
plotswarm=False
plothist=False
plotcorr=False
plottime=False  #survival
plotcorrAllcorr=False
test=False; test2=False; pstar=0.05
testtable=False
estimation=False
estimation2=False
compareinterns=False
plotClassicalGraph=False
newplots = False
allscatter=False
CompareInternParam=True
correlE=False
RuptChargeAD=False
FinalGraphs=True
fractionsPH = False
SaveGraph=True
label=''
#--------------------
# modifs php
# added 211109
REJECTFABIO=True  # TRUE: reject via dFint[ (dFint['REJECTED_CURVE_ON'] == 0 ) & (dFint['BIG_RESIDUALS_ON'] == 0 ) & (dFint['REJECTED_BYHAND_ON'] == 0 ) &( dF['BIG_RELERR_FORCEtBK1PLUS_ON']==0) ]
# FALSE : reject only by hand
#---
CleanSlopes = False
DeltaSlopePos = 10
DeltaSlopeNeg = 10

RATIOREJECT=True
#RATIOREJECTADONLY=True #FALSE : apply to AD and TU before creation of subsets ; DOES NOT WORK since division by 0 on line 398
ratio_threshold=1
#---
REMOVESMALLFORCES=True
force_threshold=6.# pN, from choice Fabio on adhesion/contact
#---
TIMERANGE=True
TIMERANGEADONLY=False
Maxtime = 1.0
#---
STACKEDSTATS=True

# comparaison rupt vs aggregation finale de donnees sur les cas sans lat
COMPAGGREG=True

#---
CLOSEFIGSATEND=True
if ListIntern==['Omar']:test2=False#pour eviter erreurs dues au manque de lat

# CREATION OF CUSTOM PALETTES FOR FINAL GRAPHS

valeurs=sns.color_palette() # palette normale
valeurs2=sns.color_palette("muted") # palette pour les +lat

molecules=['igg2a', 'cd45', 'cd3', 'cd11a open', 'cd11a closed']
test=len(molecules)*[1]

ordre_presentation=['igg2a', 'cd11a open', 'cd11a closed',  'cd3', 'cd45']
control='igg2a'

couleurs={}
couleurs2={}
j=0
for i in molecules:
    couleurs[i]=valeurs[j]
    couleurs2[i]=valeurs2[j]
    j=j+1

couleurs_sorted={}
couleurs2_sorted={}

for i in ordre_presentation:
    couleurs_sorted[i]=couleurs[i]
    couleurs2_sorted[i]=couleurs2[i]

fig, (ax1, ax2) = plt.subplots(1, 2)
sns.barplot(molecules, test, palette=valeurs, ax=ax1)
ax1.set_title("-LAT")
sns.barplot(molecules, test, palette=valeurs2, ax=ax2)
ax2.set_title("+LAT")

palette_final=dict(couleurs_sorted)
del(palette_final[control])
palette_nolat=list(palette_final.values())

palette2_final=dict(couleurs2_sorted)
del(palette2_final[control])
palette_lat=list(palette2_final.values())

palette_coupled=[]

j=0
while j < len(palette_nolat):
    palette_coupled.append(palette_nolat[j])
    palette_coupled.append(palette_lat[j])
    j=j+1
    
palette_all={}
for i in ordre_presentation:
    palette_all[i+'-Lat']=couleurs_sorted[i]
    palette_all[i+'+Lat']=couleurs2_sorted[i]
    
# -----------------------------------------------------------------------

# Definition for stacked stats and plots
def StatsDataFrame(dF, condition, type_event): # dF, 'ConditionFull', 'Morpho'

    counts=dF.groupby([condition, type_event]).agg(['count']).unstack(fill_value=0).stack() #hangling cases with 0 of this type eg.AD
    unique=dF[condition].unique(); categ=len(unique) ; names=tuple(unique)
    
    r=[i for i in range(categ)]
    adhesion=[]
    rupture=[]
    charge=[]
    condition=[]
    for i in unique:
        print(i)
        condition.append(i)
        adhesion.append(counts.loc[(i,'Ad')][('Index','count')])
        rupture.append(counts.loc[(i,'Rupt')][('Index','count')])
        charge.append(counts.loc[(i,'Char')][('Index','count')])
    raw={'condition':condition, 'adhesion':adhesion, 'rupture':rupture, 'charge':charge}
    print(raw)
    dcount=pd.DataFrame(raw)
    # let's go for stacking bars
    totals = [i+j+k for i,j,k in zip(dcount['rupture'], dcount['charge'], dcount['adhesion'])]
    print('totals:', totals)
    totalstubes=[i+j for i,j in zip(dcount['rupture'], dcount['charge'])]
    print('totals tubes:', totalstubes)
    # PEUT ETRE QQCH MYSTERIEUX LA...
    type1 = [i / j * 100 for i,j in zip(dcount['rupture'], totals)]
    type2 = [i / j * 100 for i,j in zip(dcount['charge'], totals)]
    type3 = [i / j * 100 for i,j in zip(dcount['adhesion'], totals)]
    
    fracchargrupt=[i / j for i,j in zip(dcount['charge'],totalstubes)]

                                              # plot
    barWidth = 0.5
    # Create green Bars
    figStats = plt.figure("Fractions 1", dpi=100)

    plt.bar(r, type1, color='#b5ffb9', edgecolor='white', width=barWidth, label='rupture')
    # Create orange Bars
    plt.bar(r, type3, bottom=type1, color='#f9bc86', edgecolor='white', width=barWidth, label='adhesion')
    # Create blue Bars
    plt.bar(r, type2, bottom=[i+j for i,j in zip(type1, type3)], color='#a3acff', edgecolor='white', width=barWidth, label='charge')
    # Custom x axis
    plt.xticks(r, names)
    plt.ylabel("Fraction")
    plt.xticks(rotation=90)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
    plt.tight_layout()
    plt.savefig(OutFold+'FractionsRecup.png')
    
    figStats2 = plt.figure("Fractions 2", dpi=100)
    plt.bar(r, type1, color='#b5ffb9', edgecolor='white', width=barWidth, label='rupture')
    # Create orange Bars
    plt.bar(r, type2, bottom=type1, color='#f9bc86', edgecolor='white', width=barWidth, label='charge')
    # Create blue Bars
    plt.bar(r, type3, bottom=[i+j for i,j in zip(type1, type2)], color='#a3acff', edgecolor='white', width=barWidth, label='adhesion')
    # Custom x axis
    plt.xticks(r, names)
    plt.ylabel("Fraction")
    plt.xticks(rotation=90)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
    plt.tight_layout()
    plt.savefig(OutFold+'FractionsType.png')   
    
    figStats3 = plt.figure("Fractions 3", dpi=100)

    plt.bar(r, fracchargrupt, color='grey', edgecolor='white', width=barWidth)
    # Create orange Bars
    
    # Custom x axis
    plt.xticks(r, names)
    plt.ylabel("Charge/(Charge+Rupture)")
    plt.xticks(rotation=90)
    plt.legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
    plt.tight_layout()
    plt.savefig(OutFold+'FractionsRatioCharOverRupt.png')   

print("#####################   COLUMNS    ###############################")

listrejected = ['b5c5g-2019.06.06-15.22.08.373.txt', 'b7c7d-2019.06.25-17.25.57.680.txt']
dayrejected =[] #['2021.06.25', '2021-06-25']  # Lucie bad affect of Lat

for ir, r in enumerate(listrejected):
    dFint['REJECTED_BYHAND_ON'] = np.where( (dFint['Filename']==r) , 1 , dFint['REJECTED_BYHAND_ON'])    
for ir, r in enumerate(dayrejected):
    dFint['REJECTED_BYHAND_ON'] = np.where( (dFint['myDate']==r) , 1 , dFint['REJECTED_BYHAND_ON'])    

dFint['discontType_str'] = np.where( (dFint['categoryType_str']=='adhesion') , 'Adhesion' , dFint['discontType_str'])

# addition PHP for rejecting the data following fabio criteria or not
if REJECTFABIO :
    dF = dFint[ (dFint['REJECTED_CURVE_ON'] == 0 ) & (dFint['BIG_RESIDUALS_ON'] == 0 ) & (dFint['REJECTED_BYHAND_ON'] == 0 ) & ( dFint['BIG_RELERR_FORCEtBK1PLUS_ON']==0) ]
else:
  #  dF=dFint[(dFint['REJECTED_BYHAND_ON'] == 0 ) ]
#    dF = dFint[ (dFint['REJECTED_CURVE_ON'] == 0 ) & (dFint['BIG_RESIDUALS_ON'] == 0 ) & (dFint['REJECTED_BYHAND_ON'] == 0 ) ]
    dF = dFint[ (dFint['REJECTED_BYHAND_ON'] == 0 ) ]
    
dF=dF[dF['Speed [nm/sec]']<3000]      # filtre of velocity
#dF=dF[dF['Time.break..s.']<0.3]      # filtre of timewait
if printcolumns: print(dF.columns.tolist())
#if fused_CD3: dF['myCondition'] = dF['myCondition'].str.replace(r'_ucht1', '')
if fused_CD3: dF['myCondition'] = dF['myCondition'].str.replace(r'CD3_ucht1', 'cd3')
dF['Blebbistatin'] = False

if CleanSlopes: 
    dF = dF[ ( dF['EX [pN/nm]']*dF['Speed [nm/sec]'] - dF['slope_retract_pos [pN/sec]'] > DeltaSlopePos ) & (  dF['EX [pN/nm]']*dF['Speed [nm/sec]']  - dF['slope_retract_neg [pN/sec]'] >  DeltaSlopeNeg )  ]

dF['ratio_E1E2SUMretr_E2retr recalc'] = (  ( dF['slope_retract_pos [pN/sec]'] / dF['slope_retract_neg [pN/sec]'] )
                * ( dF['EX [pN/nm]']*dF['Speed [nm/sec]'] - dF['slope_retract_neg [pN/sec]'] ) / ( dF['EX [pN/nm]']*dF['Speed [nm/sec]'] - dF['slope_retract_pos [pN/sec]'] ) )

dF['ratio_slope'] =  ( dF['slope_retract_pos [pN/sec]'] / dF['slope_retract_neg [pN/sec]'] )

dF['ratio_E1E2SUM_E2'] = ( dF["E1 [pN/nm]"] + dF["E2 [pN/nm]"]) / dF["E2 [pN/nm]"]

#dF = dF[ (dF['ratio_E1E2SUMretr_E2retr'] >= 1)  & (dF['ratio_slopeLinearPartRetrPos'] >= 1) ]        
        
        

# dF['discontTypeFused_str'] = dF['discontType_str'].str.replace(r'curvature', 'curv')
# dF['discontTypeFused_str'] = dF['discontTypeFused_str'].str.replace(r'chargeDisc', 'nocurv')
# dF['discontTypeFused_str'] = dF['discontTypeFused_str'].str.replace(r'ruptureDisc', 'nocurv')
# dF['ConditionCurv'] = np.where(dF['discontTypeFused_str']=='curv' , dF['myCondition']+'_curv', dF['myCondition']+'_nocurv')

#if OoG=='Lucie':   
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='aCD3') , 'cd3', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='Temoin_aCD3') , 'cd3', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='LatA_aCD3_2 M') , 'cd3', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='LatA_aCD3') , 'cd3', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='aCD45') , 'cd45', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='aCD11a') , 'cd11a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='Temoin_aCD11a') , 'cd11a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='LatA_aCD11a_2 M') , 'cd11a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='LatA_CD11a') , 'cd11a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='LatA_aCD11a') , 'cd11a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='aCD4') , 'cd4', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='IgG2a') , 'igg2a', dF['myCondition'])

dF['myCondition'] = np.where( (dF['InternType']==2) & (dF['myCondition']=='cd11a') , 'cd11a open', dF['myCondition'])

#if OoG=='Omar':
#dF = dF[(dF['myCondition']!='Nue') & (dF['myCondition']!='streptaNue')]   # bien définir les soussets par la suite pour éviter pollution des data
dF['Blebbistatin'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='BBCD3') , True , dF['Blebbistatin'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='streptaCD45') , 'cd45', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='streptaCD3') , 'cd3', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='streptaConA') , 'ConA', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='beadStreptaConA') , 'ConA', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='streptaIGg2A') , 'igg2a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='IGg2A') , 'igg2a', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='antiCD45') , 'cd45', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='streptaAntiCD45') , 'cd45', dF['myCondition'])
dF['myCondition'] = np.where( (dF['InternType']==1) & (dF['myCondition']=='BBCD3') , 'cd3', dF['myCondition'])

dF['myCondition'] = np.where( (dF['InternType']==0) & (dF['myCondition']=='cd11a') , 'cd11a closed', dF['myCondition'])

dF['Abs Force.contact..pN.']=-dF['Force.contact..pN.']

dF['InternName'] = ''
dF['InternName'] = np.where( (dF['InternType']==0) , 'Gautier', dF['InternName'] )
dF['InternName'] = np.where( (dF['InternType']==1) , 'Omar', dF['InternName'] )
dF['InternName'] = np.where( (dF['InternType']==2) , 'Lucie', dF['InternName'] )

dF["tDROP2-tDROP1 [sec]"] = dF["tDROP2 [sec]"]-dF["tDROP1 [sec]"]
dF['force_tBk2-tBk1 [pN]'] = dF['force_tBk2 [pN]']-dF['force_tBk1 [pN]']
dF['force_tBk+-tBk1 [pN]'] = dF['force_tBk1_plus [pN]']-dF['force_tBk1 [pN]']
dF['delta_force_tDROP1 sign [pN]'] = dF['force_tDROP1 [pN]'] - dF['force_tDROP1_plus [pN]']
dF["E1+E2 [pN/nm]"] = dF["E1 [pN/nm]"] + dF["E2 [pN/nm]"]
dF["E1byE2"] = dF["E1 [pN/nm]"] / dF["E2 [pN/nm]"]
dF['E2ESTNew [pN/nm]'] = dF['E1E2SUM_EST [pN/nm]'] - dF['E1 [pN/nm]']

dF['ConditionFull'] = np.where(dF['Latrunculine']==True , dF['myCondition']+'+Lat', dF['myCondition']+'-Lat')
dF['ConditionFullPoolLat'] = np.where(dF['Latrunculine']==True , 'All+Lat', dF['ConditionFull'])

# dF['TypeFinal'] = ''
# dF['TypeFinal'] = np.where( (dF['InternType']==0) , 'Rupt', dF['TypeFinal'] )
# dF['TypeFinal'] = np.where( (dF['InternType']==1) , 'Char', dF['TypeFinal'] )
# dF['TypeFinal'] = np.where( (dF['categoryType_str']=='adhesion') , 'AD', dF['TypeFinal'] )

dF['ConditionADTU'] = np.where(dF['categoryType_str']=='adhesion' , dF['myCondition']+'_AD', dF['myCondition']+'_TU')

dF['Tension [pN/nm]'] = dF['force_tBk2 [pN]']**2/(8*np.pi**2*2e-19)*1e-21
dF['Tension from jump [pN/nm]'] = dF['Jump.force..pN.']**2/(8*np.pi**2*2e-19)*1e-21

dF['ratio_BestSlopeNeart0ApprNeg'] =  dF['best_slope_near_t0 [pN/sec]']/dF['slope_approach [pN/sec]']

# Warning : ajoute pour eviter pb avec omar seul
if ListIntern==['Omar']:
    rejectedconditions = ['nues']
else:
    rejectedconditions = ['bCD3_aCD11a_avant', 'SCD3_aCD11a', 'cd4', 'nues', 'SCD45_aCD11a', 'bCD45_aCD11a_apres', 'bCD3_aCD11a_apres', 'bCD45_aCD11a_avant', 'igg2a']

for cond in rejectedconditions: dF = dF[dF['myCondition']!=cond]

###############################   INITIALISATION    ###################################################     

wpool=[dF]; wpoolname=["All"]

#wyM=["E1 [pN/nm]", 'E1_OverSum_EST [pN/nm]' , "E2 [pN/nm]", "etaN [pN/nm.sec]", "E1N [pN/nm]"];
wyMSchmitz=[1.5e-3, 1.5e-3, 0.2, 6e-3, 6e-3, 6e-3, 1.5e-3, 0.215]
wyMmin=[0.00001, 0.00001, 0.001, 0.0001, 0.00001, 0.00001]; wyMmax=[2., 2., 200., 10., 1., 0.1,]; wdyM=[0.05, 0.05, 0.2, 0.01, 0.01, 0.001]

wyMrelErr= ['relErr_E1', 'relErr_E1_OverSum_EST', 'relErr_E2','relErr_eta','relErr_eta','relErr_etaN','relErr_E1N', 'relErr_E1E2SUM'] #; maxerror=0.8
wyMrelErrNormSum= ['relErr_E1', 'relErr_E1_OnSum', 'relErr_E2','relErr_eta','relErr_eta','relErr_etaN','relErr_E1N', 'relErr_E1E2SUM'] #; maxerror=0.8
wyM=["E1 [pN/nm]", "E1_EST_OnSum [pN/nm]", "E2 [pN/nm]", "eta [pN/nm.sec]", "eta_EST [pN/nm.sec]", "etaN [pN/nm.sec]", "E1N [pN/nm]", "E1+E2 [pN/nm]"]
wyM_filt=["E1 [pN/nm]_filt", "E1_EST_OnSum [pN/nm]_filt", "E2 [pN/nm]_filt", "eta [pN/nm.sec]_filt", "eta_EST [pN/nm.sec]_filt", "etaN [pN/nm.sec]_filt", "E1N [pN/nm]_filt", "E1+E2 [pN/nm]_filt"]
wyMminBis=[1e-5, 1e-5, 1e-5, 0.001, 0.001, 0.001, 0.00001, 0.01]; wyMmaxBis=[2., 2., 2., 10, 10., 0.1, 0.01, 2.]; wdyMBis=[0.1, 0.1, 0.1, 0.5, 0.5, 0.01, 0.001, 0.1]
wmaxerror = [0.5, 0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.3]
wmaxerror = [10,10,10,10,10,10,10,10]
wyM2=['E1 [pN/nm]', 'E1_EST [pN/nm]', 'E1_EST_OnSum [pN/nm]', 'E1_EST_OnDeltaForce [pN/nm]', "E2 [pN/nm]", 'E2_EST [pN/nm]', "E1+E2 [pN/nm]"]
wyM2minBis=[1e-5, 1e-5, 1e-5, 1e-5, 0.01, 0.01, 0.01]; wyM2maxBis=[2., 2., 2., 2., 1, 1, 2.]; wdyM2Bis=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

for m, param in enumerate(wyM): dF[param+'_filt'] = np.where(dF[wyMrelErrNormSum[m]]<wmaxerror[m], dF[param], np.nan)

wyETA = [ "eta [pN/nm.sec]_filt", "eta_EST [pN/nm.sec]_filt", 'eta_EST_04 [pN/nm]', "eta_E1NZero_04_EST [pN/nm]", "etaN [pN/nm.sec]_filt" ]
wyETAmin=[0.001, 0.001, 0.001, 0.001, 0.001]; wyETAmax=[10., 10., 10., 10., 10.]; wdyETA=[1, 1, 1, 1, 0.01]


wyComb = ['best_slope_near_t0 [pN/sec]', 'slope_near_t0_THEO [pN/sec]', 'E1+E2 [pN/nm]', "E1 [pN/nm]", "E2 [pN/nm]"]
wyCombmin=[0. , 0. , 0.01, 1e-5, 1e-5]; wyCombmax=[ 750., 750., 2., 2., 2., ]; wdyComb=[100.,100., 0.1, 0.1, 0.1]

wyT1=["tBk1 [sec]", 'tDROP1 [sec]', 'tDROP2 [sec]']; wyT1Schmitz=[0., 0., 0.]
wyT1min=[0.01, 0.01, 0.01]; wyT1max=[10., 10., 10.]; wdyT1=[0.5, 0.5, 0.5]

wyT2=['tauB_THEO [sec]', 'tCHARGE [sec]']
wyT2min=[0.01, 0.01]; wyT2max=[10., 10.]; wdyT2=[0.5, 0.5]

wyF1=["force_tBk1 [pN]", 'force_tBk1_THEO [pN]', 'force_tBk2 [pN]', 'force_tBk2_THEO [pN]']
wyF1min=[0., 0., 0., 0.]; wyF1max=[100., 100., 100., 100.]; wdyF1=[30., 30., 30., 30.]

wyS=['slope_near_t0 [pN/sec]', 'slope_near_t0_THEO [pN/sec]', 'slope_tube [pN/sec]', 'slope_tube_THEO [pN/sec]']
wySmin=[0. , 0. , -10, -10.]; wySmax=[ 750., 750., 30., 30.]; wdyS=[100.,100., 5., 5.]

wyT = wyT1+wyT2; wyTmin = wyT1min+wyT2min; wyTmax = wyT1max+wyT2max; wdyT = wdyT1+wdyT2

ordrepres=['cd3', cd11a0rConA, 'cd45'] if fused_CD3 else ['cd3', 'cd3_ucht1', cd11a0rConA, 'cd45']

ordredisc=["curvature", "chargeDisc",  "ruptureDisc" ]

###############################   DEFINITION DES SOUS SETS    ###################################################     
#------------------------------
# REMOVING RESIDUAL SMALL FORCES (Fabio)
#------------------------------

if REMOVESMALLFORCES:
    dF=dF[dF['force_tBk1 [pN]']>force_threshold]

# DEFINING SUBSETS

AD = dF[dF['categoryType_str']=='adhesion']
TF = dF[dF['categoryType_str']=='finTube']
TI = dF[dF['categoryType_str']=='infinTube']
TU = dF[(dF['categoryType_str']=='finTube')+(dF['categoryType_str']=='infinTube')]

AD['E1N [pN/nm]'] = np.nan # np.where( (dF['InternType']==1) & (dF['myCondition']=='BBCD3') , True , dF['Blebbistatin'])

#------------------------------
# SORTING DATA ON RATIO (Fabio)
#------------------------------


if RATIOREJECT:
    # if RATIOREJECTADONLY:

    AD=AD[AD['ratio_E1E2SUMretr_E2retr']>ratio_threshold]
    # nouvelle condition imposee sur la recup de data d'adhesion (29/11/21)
    AD = AD[ (AD["E2_EST [pN/nm]"]-AD["E2 [pN/nm]"])/AD["E2_EST [pN/nm]"] < 0.1] # on ne prend que les E2=E2_est dans le pool
    # else:
    #     AD=AD[AD['ratio_E1E2SUMretr_E2retr']>ratio_threshold]
    #     TU=TU[TU['ratio_E1E2SUMretr_E2retr']>ratio_threshold]

#------------------------------
# REMOVING ZERO TIME
#------------------------------

if TIMERANGE:
    if TIMERANGEADONLY:
        AD=AD[ (AD['Time.break..s.']>=0) & (AD['Time.break..s.']<=Maxtime) ]   # 'deltaTimeWait [sec]'
    else :
        AD = AD[ (AD['Time.break..s.']>=0) & (AD['Time.break..s.']<=Maxtime) ]
        TU = TU[ (TU['Time.break..s.']>=0) & (TU['Time.break..s.']<=Maxtime) ]
# medTUbasaltension = TU[TU['Latrunculine']==True].median()['Tension [pN/nm]']
# TU['W0 [pN/nm]'] = TU['Tension [pN/nm]']-medTUbasaltension

RUPT = TU[TU['discontType_str']=='ruptureDisc']
CHAR = TU[TU['discontType_str']=='chargeDisc']

TUNL = TU[TU['Latrunculine']==False]
ADNL = AD[AD['Latrunculine']==False]

RUPTFree = RUPT[ (RUPT['Constraint']==False)]
CHARFree = CHAR[ (CHAR['Constraint']==False)]
#ADCons = AD[ (AD['Constraint']==True)]     # 20211025 in case of 2 separate tables according to  
ADCons = AD[ (AD['Constraint']==False)]    

RUPTFree['Morpho'] = 'Rupt'
CHARFree['Morpho'] = 'Char'
ADCons['Morpho'] = 'Ad'

RUPTFreeLat = RUPTFree[RUPTFree['Latrunculine']==True]
CHARFreeLat = CHARFree[CHARFree['Latrunculine']==True]

# ratio_E1E2SUMretr_E2retr ≥ ratio_E1E2SUMretr_E2retr_Tol & rehabAdhE1min_cond > E1_min_Tol & rehabAdhE1max < E1_max_Tol );

# ratio_E1E2SUMretr_E2retr_Tol = 1.1;

# E1_min_Tol = 5*1E-3;

# E1_max_Tol = 9*1E1;

u2nolat = RUPTFree[RUPTFree['Latrunculine']==False]["E2 [pN/nm]"]
E2_75 =  np.quantile(u2nolat, 0.75);# E1_max = np.quantile(u1nolat, 0.75)
print('E2_75 from RUPT=', E2_75)
print('RUPTFree["E2 [pN/nm]"]<=E2_75', RUPTFree[RUPTFree["E2 [pN/nm]"]<=E2_75].count()["E2 [pN/nm]"])
print('RUPTFree', len(RUPTFree))
print('RUPTFreeLat[RUPTFreeLat["E2 [pN/nm]"]<=E2_75', RUPTFreeLat[RUPTFreeLat["E2 [pN/nm]"]<=E2_75].count()["E2 [pN/nm]"])
print('RUPTFreeLat', len(RUPTFreeLat))

u3nolat = CHARFree[CHARFree['Latrunculine']==False]["E2 [pN/nm]"]
E2_75ch =  np.quantile(u3nolat, 0.75);# E1_max = np.quantile(u1nolat, 0.75)
print('E2_75 from CHAR=', E2_75ch)

# E1_EST_min = 1e-4; E1_min = 0.003
E1_min = 0.003; E1_max = 0.5 ; E2_max = 1.
#ADcut = ADCons[ (ADCons['E2_EST_OnRetraction [pN/nm]']==ADCons['E2_EST [pN/nm]']) & (ADCons['E1 [pN/nm]']>E1_min) & (ADCons['E1_EST [pN/nm]']>E1_EST_min) ]
#ADcut = ADCons[ (ADCons['E2_EST_OnRetraction [pN/nm]']==ADCons['E2_EST [pN/nm]']) & (ADCons['E1 [pN/nm]']>E1_min) ]
# ADcut = ADCons[ (ADCons['ratio_E1E2SUMretr_E2retr']>1) & (ADCons['E1 [pN/nm]']>E1_min) ]
#ADcut = ADCons[ (ADCons['E1 [pN/nm]']>E1_min) &  (ADCons['E2 [pN/nm]']>E1_min) &  (ADCons['E2 [pN/nm]']<0.1) ]

ADcut = ADCons[ (ADCons['E1 [pN/nm]']>E1_min) &  (ADCons['E1 [pN/nm]']<E1_max) &  (ADCons['E2 [pN/nm]']<E2_max)]
#RUPTcut = RUPTFree[ (RUPTFree['E2_EST_OnRetraction [pN/nm]']==RUPTFree['E2_EST [pN/nm]']) & (RUPTFree['E1 [pN/nm]']>E1_min) & (RUPTFree['E1_EST [pN/nm]']>E1_EST_min) ]

ADcutLat = ADcut[ADcut['Latrunculine']==True]
#print(len(ADcut))
print('ADcutLat[ADcutLat["E2 [pN/nm]"]<=E2_75]', ADcutLat[ADcutLat["E2 [pN/nm]"]<=E2_75].count()["E2 [pN/nm]"])
print('ADcutLat', len(ADcutLat))

ADcut['E1 [pN/nm]'] = np.where( (ADcut["E2 [pN/nm]"]>E2_75) & (ADcut['Latrunculine']==True) , np.nan, ADcut["E1 [pN/nm]"] )
ADcut['eta [pN/nm.sec]'] = np.where( (ADcut["E2 [pN/nm]"]>E2_75) & (ADcut['Latrunculine']==True) , np.nan, ADcut['eta [pN/nm.sec]'] )
ADcut['etaN [pN/nm.sec]'] = np.nan
ADcut['E1N [pN/nm]'] = np.nan
ADcut['E2 [pN/nm]'] = np.where( (ADcut["E2 [pN/nm]"]>E2_75) & (ADcut['Latrunculine']==True) , np.nan, ADcut["E2 [pN/nm]"] )
#print(len(ADcut))

#u1nolat = RUPTFree[RUPTFree['Latrunculine']==False]['eta [pN/nm.sec]']
# ETA_75 =  np.quantile(u1nolat, 0.75);# E1_max = np.quantile(u1nolat, 0.75)
# ETA_25 =  np.quantile(u1nolat, 0.25);# E1_max = np.quantile(u1nolat, 0.75)

ADtime = ADcut[ (ADcut['deltaTimeWait [sec]']>=0) & (ADcut['deltaTimeWait [sec]']<=0.5) ] 


RUPTAD = pd.concat([RUPTFree, ADcut], ignore_index=True)
RUPTADCHAR = pd.concat([RUPTFree, ADcut, CHARFree], ignore_index=True)
RUPTADCHARLat = RUPTADCHAR[RUPTADCHAR['Latrunculine']==True]
RUPTCHAR = pd.concat([RUPTFree, CHARFree], ignore_index=True)

RUPTADnolim = pd.concat([RUPTFree, ADCons], ignore_index=True)

RUPTADCHARnolim = pd.concat([RUPTFree, ADCons, CHARFree], ignore_index=True)


# descriptive stats to know how much we gain
print('=============================================')
print('Descriptive stats')
for i, j in zip(['RUPTFree', 'RUPTAD', 'RUPTCHAR', 'RUPTADCHAR'], [RUPTFree, RUPTAD, RUPTCHAR, RUPTADCHAR]):
        print('-------------------------------------')
        print(i)
        dh=j.groupby(['Morpho', 'ConditionFull']).describe()
        print(dh)
        if i=='RUPTADCHAR':
            dfi.export(dh, OutFold+'RuptADChar-Counts.png', max_cols=1)

# ---------------------------------------------------
# Warning : on a du ajouter ceci pour eviter un pbv pour omar seul
if ListIntern==['Omar']:
    ordrepres3=['igg2a', 'cd45', 'cd3', 'ConA']
else:
    ordrepres3 = ['cd11a open', 'cd11a closed', 'cd3', 'cd45']
# ---------------------------------------------------

    
wp = ["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
wmed = np.zeros((5, 4))
for ip, p in enumerate(wp):
    RUPTAD['Delta'+p] = np.nan
    for isset, sset in enumerate(ordrepres3):
        dFx = RUPTAD[ (RUPTAD['myCondition']==sset) & (RUPTAD['Latrunculine']==True) ]
        wmed[ip, isset] = dFx[p].median()
        RUPTAD['Delta'+p] = np.where( (RUPTAD['myCondition']==sset) & (RUPTAD['Latrunculine']==False) , RUPTAD[p] - wmed[ip, isset], RUPTAD["Delta"+p] )
        print(sset, p, wmed[ip, isset], RUPTAD[RUPTAD['myCondition']==sset]['Delta'+p].median())

wpA = ['E1 [pN/nm]', 'eta [pN/nm.sec]']
wpB = ['E1N [pN/nm]', 'etaN [pN/nm.sec]']
wmed2 = np.zeros((2, 4))

for ip, p in enumerate(wpA):
    RUPTAD['Jump'+p] = np.nan
    for isset, sset in enumerate(ordrepres3):
        dFx = RUPTAD[ (RUPTAD['myCondition']==sset) ]
        wmed2[ip, isset] = dFx[wpB[ip]].median()
        RUPTAD['Jump'+p] = np.where( (RUPTAD['myCondition']==sset) , RUPTAD[p] - wmed2[ip, isset], RUPTAD["Jump"+p] )
        print(sset, p, wmed[ip, isset], RUPTAD[RUPTAD['myCondition']==sset]['Jump'+p].median())        
 
#ADCons = AD[ (AD['Constraint']==True) & (AD['E1 [pN/nm]']>E1_min) & (AD['E1_EST [pN/nm]']>E1_EST_min)]    
# u1nolat = RUPTFree[RUPTFree['Latrunculine']==False]["E1 [pN/nm]"]
#u1lat = RUPTFree[RUPTFree['Latrunculine']==True]["E1 [pN/nm]"]
# u2lat = RUPTFree[RUPTFree['Latrunculine']==True]["E2 [pN/nm]"]
 
# E2_min = np.quantile(u2lat, 0.25); E1_min = np.quantile(u1lat, 0.25)
#print(len(RUPT["E2 [pN/nm]"]))


#ADcut = AD[ (AD['E2 [pN/nm]']>E2_min) & (AD['E2 [pN/nm]']<E2_max) & (AD['E1 [pN/nm]']>E1_min) & (AD['E1 [pN/nm]']<E1_max) & (AD['E1_EST [pN/nm]']>E1_EST_min)]

wpoolM5=[AD, TF, RUPT]; wpoolM5name=["Adhesions", "Finites Tubes", "Rupture Disc"]

#wpoolM=[AD, TU, CU, CN]; wpoolMname=["Adhesions", "Tubes", "Tubes Curvature", "Tubes Nocurv"]
#wpoolC=[CU, CN, RUPT]; wpoolCname=["Tubes Curvature", "Tubes Nocurv", "Tubes RuptureDisc"]

wpoolTU=[TU]; wpoolTUname=["Tubes"]

wpoolM2=[AD, TF]; wpoolM2name=["Adhesions", "Finites Tubes"]
wpoolM3=[AD, TU, dF]; wpoolM3name=["Adhesions", "Tubes", "All"]
wpoolM4=[AD, TF, TI, dF]; wpoolM4name=["Adhesions", "Fin. Tubes", "Inf. Tubes", "All"]
    
a = 0.1
RUPTGoodE2EST = RUPT[ (RUPT['E2_EST [pN/nm]'] > (1-a)*RUPT['E2 [pN/nm]']) &  (RUPT['E2_EST [pN/nm]'] < (1+a)*RUPT['E2 [pN/nm]']) & (RUPT['Constraint']==False) ]
wtwait = [(0,np.inf), (0, 0.5), (0.5, 1.0), (1.01, np.inf), (0.05, 1)]
for twait in wtwait:
    strwait = 'deltaTimeWait [sec]'
    RUPTt = RUPT[ (RUPT[strwait]>=twait[0]) & (RUPT[strwait]<=twait[1]) & (RUPT['Constraint']==False)]
    RUPTGoodE2ESTt = RUPTGoodE2EST[ (RUPTGoodE2EST[strwait]>=twait[0]) & (RUPTGoodE2EST[strwait]<=twait[1])]
    print(a, str(twait), len(RUPTGoodE2ESTt), len(RUPTt))
    if len(RUPTt)>0: print( len(RUPTGoodE2ESTt)/len(RUPTt) )  
        
##########################################   FUNCTIONS    ###################################################     

def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return fig, ax        

# FractionPlot(r, counts, 'NAd', 'AD', 'TU', durees,k, unique, names, t, angle)
def FractionPlot(r, counts, df1, df2, df3, durees,k, unique, names, t, angle):
# def FractionPlot(df1, df2, df3, df4, unique, names, t, angle):
    raw={}
    noadh=[] 
    adhesion=[]
    tube=[]
    # infinite=[]
    totals=[]
    #print(r, names)    listemod=tuple(tmp)
    for j in unique:
        if counts.index.isin([(j,df1,t)]).any():
            noadh.append(counts.loc[(j,df1,t)][('X','count')])
        else:
            noadh.append(0)
        if counts.index.isin([(j,df2,t)]).any():
            adhesion.append(counts.loc[(j,df2,t)][('X','count')])
        else:
            adhesion.append(0)
        if counts.index.isin([(j,df3,t)]).any():
            tube.append(counts.loc[(j,df3,t)][('X','count')])
        else:
            tube.append(0) 
        # if counts.index.isin([(j,df4,t)]).any():
        #     infinite.append(counts.loc[(j,df4,t)][('X','count')])
        # else:
        #     infinite.append(0)
    raw={'no':noadh, 'adh':adhesion, 'tub':tube}#, 'inftub':infinite}
    #print(raw)
    dcount=pd.DataFrame(raw)
    # sums
    totals = [i+j+k for i,j,k in zip(dcount['no'], dcount['adh'], dcount['tub'])]#, dcount['inftub'])]
    # ratios - COMMENT : np.divide(x,0) renvoit 0 !!! pout Trueeviter les cas ou il y a pas de mesures
    no = np.nan_to_num([np.divide(i, j) * 100 for i,j in zip(dcount['no'], totals)])
    # print(no)
    adh = np.nan_to_num([np.divide(i, j)  * 100 for i,j in zip(dcount['adh'], totals)])
    # print(adh)
    tub =np.nan_to_num([np.divide(i, j)  * 100 for i,j in zip(dcount['tub'], totals)])
    # print(tub)
    # inftub = np.nan_to_num([np.divide(i, j)  * 100 for i,j in zip(dcount['inftub'], totals)])
    # print(inftub)
    barWidth = 0.5
    # Create green Bars : tubes
    # ax[k].bar(r, inftub, color='darkgreen', edgecolor='white', width=barWidth, label='Inf. Tubes')
    # Create orange Bars : adhesion
    ax[k].bar(r, tub, color='lightgreen', edgecolor='white', width=barWidth, label='Tubes')
    # Create blue Bars : non adhesion
    ax[k].bar(r, adh, bottom=tub, color='orange', edgecolor='white', width=barWidth, label='Adhesion')
    
    ax[k].bar(r, no, bottom=[i+j for i,j in zip(tub, adh)], color='lightblue', edgecolor='white', width=barWidth, label='No')
    # Custom x axis
    ax[k].set_xlabel("t="+str(t)+"sec")
    ax[k].set_xticks(r)
    ax[k].set_xticklabels(names, rotation=angle)
    ax[k].set_ylim(0,100)
    if k==durees-1:
        ax[k].legend(loc='upper left', bbox_to_anchor=(1,1), ncol=1)
    if k==0:
        ax[k].set_ylabel("Fraction")
    
    plt.tight_layout()

def MultiSwarmPlots(wy, wymax, wymin, wpool, colhue, figname):
    fig = plt.figure(figname, figsize=(len(wpool)*2, len(wy)*2), dpi=100)
    for (iy,y) in enumerate(wy):
        for (ipool,pool) in enumerate(wpool):
            fig.add_subplot( len(wy), len(wpool),len(wpool)*iy+ipool+1)
            if colhue=="": ax = sns.swarmplot(x="myCondition", y=y, data=pool,
                                              dodge=True, order=ordrepres)
            else: ax = sns.swarmplot(x="myCondition", y=y, hue=colhue, data=pool,
                                     dodge=True, order=ordrepres)
            ax.set_ylim(wymin[iy],wymax[iy]); ax.set_ylabel(y)
            if iy!=len(wy)-1: ax.axes.get_xaxis().set_visible(False)
            if ipool!=0: ax.axes.get_yaxis().set_visible(False)
            if iy==0: ax.set_title(wpoolname[ipool]+' n='+str(len(y)), fontsize=6)
            if not "[sec]" in y:  ax.set_yscale('log')
            ax.legend(loc="upper right", title=colhue, title_fontsize=5, prop={'size': 5})
            
def MultiHistPlots(axall, wy, wymax, wymin, wyRef, wpool, wpoolname, namecond, condCD, condLat, figname):
    fig.suptitle(figname, fontsize=14)
    for (iy,y) in enumerate(wy):
        if not "[sec]" in y:  bins = np.logspace(np.log10(wymin[iy]), np.log10(wymax[iy]), num=40)
        else: bins = np.linspace(wymin[iy], wymax[iy], num=40)
   #     print('pool length=', len(wpool))
        for (ipool, upool) in enumerate(wpool):
            if condCD=='':
                pool=upool[upool['Latrunculine']==(condLat=='Lat')]
            else:
                pool=upool[(upool[namecond]==condCD)&(upool['Latrunculine']==(condLat=='Lat'))]
            if len(wpool)==1: ax = axall[iy]
            else: ax = axall[ipool, iy]
            sns.distplot(pool[y], kde=False, axlabel=y, bins=bins, 
                         label=condCD+condLat+':'+str(pool[y].count()), ax=ax)

            if iy==0: print(condCD, condLat, wpoolname[ipool], pool[y].count() )
#            plt.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
            ax.set_xlim(wymin[iy],wymax[iy])
            if not "[sec]" in y: ax.set_xscale('log')
            ax.set_ylim(0,max(pool[y].count()/4.,30)); ax.set_ylabel(wpoolname[ipool])
            if ipool!=len(wpool)-1: ax.axes.get_xaxis().set_visible(False)       
            if iy!=0: ax.axes.get_yaxis().set_visible(False)
            #ax2=ax.twinx()
            if wyRef!=None: ax.plot([wyRef[iy],wyRef[iy]], [0,30], linestyle='--',
                                    c='k',  label='Schmitz '+ "%.4f" % wyRef[iy])
            ax.plot([pool[y].median(), pool[y].median()], [0,30], 
                         alpha=0.5, label='Med '+condLat+ " %.5f" % pool[y].median())
            print(condLat, wpoolname[ipool], y, 'N=', pool[y].count(), 'Median=',
                  " %.5f" % pool[y].median())
   #         if iy==0:
            ax.legend(loc="upper right", title_fontsize=4, prop={'size': 5})
            
def MultiCorr(wy, wpool, wpoolname, condCD, condLat, figname0):
    ratio=0.5
    for (ipool, upool) in enumerate(wpool):
        #print('ipool=', ipool)
        figname=figname0+wpoolname[ipool]
        pool=upool[upool['myCondition'].str.contains(condCD)&
                   (upool['Latrunculine']==(condLat=='Lat'))]
#https://spam.inserm.fr/fmlurlsvc/?fewReq=:B:JVUzOD85My9/NDsnOS9gbTQ5ODM5OC96YG5naH18e2w0OD45amtsbz9sPDs4ODFqPzE9Oz8wMD86Om06PTg4O2s/OTg7bz0/ay99NDg+MT04OTwxPTwveGBtND8/TzF/WWBmOTs/OD4/JD8/TzF/WWB4OTs/OD4/L3tqeX00eWBse3tsJGFsZ3tgJ3l8bGphSXxnYH8kaGR8J297L2o0PDsvYW1lNDk=&url=https%3a%2f%2fseaborn.pydata.org%2fexamples%2fmany_pairwise_correlations.html 
        corr = pool[wy].corr(min_periods=5)  #      print(corr)
        mask = np.triu(np.ones_like(corr, dtype=np.bool))
        f, ax = plt.subplots(num=figname, figsize=(11*ratio,9*ratio))
        f.suptitle(figname, fontsize=14)
        cmap = sns.diverging_palette(220, 10, as_cmap=True)
        sns.heatmap(corr, mask=mask, cmap=cmap, vmin=-1., vmax=1., annot=True,
                    center=0, square=True, linewidths=.5,cbar_kws={"shrink": .5}) 
        ax.tick_params(axis='both', which='major', labelsize=5)
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

def MultiEstimate(wy, wymax, wymin, wdy, pool, xCond, idx, MeanMed, CI=0.95, figname='Default', logscale=True):
    fig , axall = plt.subplots(1, len(wy), num=figname, figsize=(len(wy)*3.5, 6.),
                               gridspec_kw={'wspace': 0.25}, dpi=100)
    fig.suptitle(figname, fontsize=14)
    genwy = ((iy,y) for (iy,y) in enumerate(wy))
    for (iy,y) in genwy:
        print(figname, y)
        ax = axall[iy]
        try:
            multi_groups = dabest.load(pool, idx=idx, x=xCond, y=y) #dabest.load(pool, idx=idx, x=xCond, y=y, ci=CI)
            if logscale and (not ("[pN]" in y or "[pN/sec]" in y or "[pN/nm]" in y or "[sec]" in y or "[nm/sec]" in y)) :  ax.set_yscale('log')
            if logscale :  ax.set_yscale('log')
            if MeanMed=='Mean': multi_groups.mean_diff.plot(ax=ax)
            elif MeanMed=='Median': multi_groups.median_diff.plot(fig_size=(3.5, 6.), ax=ax)
            ax.set_title(y, fontsize=10); ax.set_ylabel('Parameter value')
            if iy!=0: ax.set_ylabel(''); ax.contrast_axes.set_ylabel('')
            ax.xaxis.label.set_size(4); ax.tick_params(axis='x', which='major', labelsize=5)
            ax.contrast_axes.tick_params(axis='x', which='major', labelsize=7)
            ax.tick_params(axis="y", labelsize=8)
            ax.contrast_axes.tick_params(axis="y", labelsize=8)
            ax.set_ylim(wymin[iy],wymax[iy])
            ax.contrast_axes.set_ylim(-wdy[iy],wdy[iy])
     #       if not ("[sec]" in y or "[pN]" in y or "[pN/sec]" in y) :  ax.set_yscale('log')
        except (IndexError, KeyError, ValueError):
            print("Index or Key or Value Error")#; genwy.__next__()
 #       perm_test = dabest.PermutationTest(control, test, effect_size="median_diff", is_paired=False)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    return fig


def SimpleEstimate(ax, y, iy, wymax, wymin, wdy, pool, xCond, idx, palette, MeanMed, CI=0.95,figname='Default', logscale=True):
    # fig , axall = plt.subplots(1, len(wy), num=figname, figsize=(len(wy)*3.5, 6.),
    #                            gridspec_kw={'wspace': 0.25}, dpi=100)
    # fig.suptitle(figname, fontsize=14)
    print(figname, y)
    try:
        multi_groups = dabest.load(pool, idx=idx, x=xCond, y=y) #dabest.load(pool, idx=idx, x=xCond, y=y, ci=CI)
        if logscale and (not ("[pN]" in y or "[pN/sec]" in y or "[pN/nm]" in y or "[sec]" in y or "[nm/sec]" in y)) :  ax.set_yscale('log')
        if logscale :  ax.set_yscale('log')
        if MeanMed=='Mean': multi_groups.mean_diff.plot(custom_palette=palette, ax=ax)
        elif MeanMed=='Median': multi_groups.median_diff.plot(fig_size=(3.5, 6.), custom_palette=palette, ax=ax)
        ax.set_title(y, fontsize=10); ax.set_ylabel('Parameter value')
        if iy!=0: ax.set_ylabel(''); ax.contrast_axes.set_ylabel('')
        ax.xaxis.label.set_size(4); ax.tick_params(axis='x', which='major', labelsize=5)
        ax.contrast_axes.tick_params(axis='x', which='major', labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
        ax.contrast_axes.tick_params(axis="y", labelsize=8)
        ax.set_ylim(wymin[iy],wymax[iy])
        ax.contrast_axes.set_ylim(-wdy[iy],wdy[iy])
    except (IndexError, KeyError, ValueError):
        print("Index or Key or Value Error")#; genwy.__next__()
 #   if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    return fig

def DoFit(ax, f, wx, wy, xmax, initialguessEq):
    x2 = wx[wx<xmax]
    y2 = wy[wx<xmax]   #    y2log = dur2s_logy[dur2s<maxtft]
    while True:
        try:
            pEq, pcovEq=curve_fit(f, x2, y2, p0=initialguessEq)
            break
        except RuntimeError:
            print("No convergence"); pEq=np.zeros(len(initialguessEq))
            break
    y2_fit=f(x2, pEq[0], pEq[1])
    ax.plot(x2, y2_fit, alpha=0.5, label='p0='+"%.5f" % pEq[0]+'/ p1='+"%.5f" % pEq[1])
    return pEq[0], pEq[1]

k0=1.; kBT=4; d=1.  # kinetic parameters
def FitExp1(x, *p): return np.exp(-(x-p[1])/p[0])
def Qr_(t, *p): return np.exp(-p[0]*(np.exp(t*p[1])-1))
def Qf_(t, *p): return np.exp(-t*p[0])

# listx=['best_slope_near_t0 [pN/sec]', 'force_tBk2 [pN]']; listy=["tBk1 [sec]", "tBk2 [sec]"]
def MultiSurvival(axall, axallparam, axallft, icondCD, ncondCD, listx, listz, wpool, wpoolname, nbins, condCD, condLat, figname):
    if nbins==3: axallparam[1, icondCD].set_prop_cycle(color=['red', 'green', 'blue', 'red' ,'green', 'blue'])
    if nbins==1: axallparam[1, icondCD].set_prop_cycle(color=['red', 'red'])
    for (ipool, upool) in enumerate(wpool):        
         pool=upool[upool['myCondition'].str.contains(condCD)&
                    (upool['Latrunculine']==(condLat=='lat'))]
         ax = axall[ipool, icondCD]
         axparam = axallparam[ipool, icondCD]
         axft = axallft[ipool, icondCD]
         if ipool==0 or ipool==2: axparamR=axparam.twinx()
         x=listx[ipool]; z=listz[ipool]   # x="tBk1 [sec]"; z="Slope [pN/nm]"
         p=pool.sort_values(by=[z])
         p["bins"] = pd.qcut(p[z], nbins, labels=False, duplicates='drop')#   print(p["bins"])
         wz=np.zeros(nbins); wdz=np.zeros(nbins)
         wp0=np.zeros(nbins); wp1=np.zeros(nbins)
         binsForce_tBk2 = np.linspace(0, 60, num=30)
         print(wpoolname[ipool], condCD,'z=', z)
         for n in range(nbins):
             zn=p["bins"]==n; wz[n]=p[zn][z].median(); wdz[n]=p[zn][z].std()
             wx=np.sort(p[zn][x]); wy=1-np.arange(len(wx))/len(wx)
            # print(ipool, n, len(wx), len(wy))
             ax.scatter(wx, wy, alpha=0.5, label=wpoolname[ipool]+' N='+str(len(wx))+' '+"%.5f" % wz[n])
             if (ipool==0 or ipool==2) and len(wx)>1:
                 wp0[n], wp1[n] = DoFit(ax, Qr_, wx, wy, 0.5, [1., 0.2])             
             if ipool==1 and len(wx)>1 :
                 wp0[n], wp1[n] = DoFit(ax, Qf_, wx, wy, 7., [1., 0.2])                   
             print('wz=', wz); print('wdz=', wdz); print('wp0=', wp0); print('wp1=', wp1)
         if ipool==0 or ipool==2:
             wk0=wp0*wp1; wf0=wz/wp1  # caution use real speed
             for n in range(nbins):
                 wkf = wk0[n]*np.exp(binsForce_tBk2/wf0[n])
                 if ipool==0: axallparam[1, icondCD].plot(binsForce_tBk2, wkf, '-', label='Out_Adh'+str(n))
                 if ipool==2: axallparam[1, icondCD].plot(binsForce_tBk2, wkf, '--', label='Intra'+str(n))
             print('wk0=',wk0);print('wf0=',wf0)
             xt = pool["tBk1 [sec]"]; yf = pool['force_tBk1 [pN]']
             axft.set_ylabel(wpoolname[ipool]+'\n'+'force_tBk1 [pN]'); axft.set_xlabel('tBk1 [sec]')
             axft.axis([0., 1.5, 0., 120.])
         if ipool==1:
             wkf=wp0 #wkf=wp0*np.exp(-wz/wf0)
             print('wkf=',wkf)
             xt = pool["tBk2 [sec]"]; yf = pool['force_tBk2 [pN]']
             axft.set_ylabel(wpoolname[ipool]+'\n'+'force_tBk2 [pN]'); axft.set_xlabel('tBk2 [sec]')
             axft.axis([0., 10., 0., 120.])
         axft.scatter(xt, yf, alpha=0.5, label='wk0')
         
         ax.set_xlabel(x, fontsize=7); ax.set_ylabel(wpoolname[ipool]+'\n'+'Survival')
         if ipool==0: ax.set_title(condCD, fontsize=7)
         ax.set_yscale('log'); ax.axis([0., 0.5, 0.05, 1.])
         if ipool==0 or ipool==2: ax.axis([-0.02, .75, 0.01, 1.5])
         if ipool==1: ax.axis([-0.2, 10., 0.01, 1.5])
         ax.legend(loc="lower left", title='Number', title_fontsize=5, prop={'size': 5})
         if icondCD!=0: ax.axes.get_yaxis().set_visible(False); axft.axes.get_yaxis().set_visible(False)
         if ipool==0 or ipool==2:
             axparam.scatter(wz, wk0, alpha=0.5, label='wk0')
             axparamR.scatter(wz, wf0, marker='+', c='r', alpha=0.5, label='wf0')
             axparam.axis([0., 500, 0.0001, 100]); axparam.set_yscale('log')
             axparam.set_title(condCD, fontsize=7)
        #     if icondCD==ncondCD-1: axparamR.legend(loc="upper right")
             axparamR.axis([0., 500, 0., 80])
             axparam.set_ylabel(wpoolname[ipool]+'\n'+'Off-rate@ 0pN'); axparamR.set_ylabel('Bell force (pN)')
             if icondCD!=ncondCD-1: axparamR.axes.get_yaxis().set_visible(False)       
         if ipool==1:
             axparam.scatter(wz, wkf, alpha=0.5, label='wkf')
             axparam.axis([0., 60., 0.0001, 100]); axparam.set_yscale('log')
             axparam.set_ylabel(wpoolname[ipool]+'\n'+'Off-rate')         
         axparam.set_xlabel(z, fontsize=7)
         axparam.grid()
         if icondCD==ncondCD-1: axparam.legend(loc="center right"); axparamR.legend(loc="upper right")
    #     if icondCD!=0: axparam.axes.get_yaxis().set_visible(False)
         axallparam[1, icondCD].legend(loc="lower right", prop={'size': 7})
         
    return (wz, wdz, wp0, wp1, wk0, wf0, wkf)


##########################################   MAIN    ###################################################     

if countpopulations:
    listgb0=['myCondition','categoryType_str']
    listgb1=['ConditionFull','categoryType_str']    
    listgb2=['categoryType_str','discontType_str']
    listgb3=['InternName', 'myCondition']
    listgb5=['discontType_str','myCondition']
    listgb6=['discontType_str','Latrunculine']
    listgb7=['ConditionFull', 'discontType_str']
    multilist = [listgb0, listgb1, listgb2, listgb3, listgb5, listgb6, listgb7]
    dFsort = dF.sort_values(by=['InternName', 'ConditionFull', 'categoryType_str', 'discontType_str'])
    ADsort = dFsort[dFsort['categoryType_str']=='adhesion']
    TFsort = dFsort[dFsort['categoryType_str']=='finTube']
    TIsort = dFsort[dFsort['categoryType_str']=='infinTube']
    TUsort = dFsort[(dFsort['categoryType_str']=='finTube')+(dF['categoryType_str']=='infinTube')]
    wpoolM3sort=[ADsort, TUsort, dFsort]
    for (ipool,pool) in enumerate(wpoolM3sort):
        figname = 'MultiCount_'+wpoolM3name[ipool]+'_'+OoG
        fig, axall= plt.subplots(4,2, num=figname, figsize=(10,12), dpi=100); fig.suptitle(figname, fontsize=14) 
        for igb, ngb in enumerate(multilist):            
            print("########################   COUNT    #########################", wpoolM3name[ipool])
            result = pool.groupby(ngb).size()
            print(result)
            ax=axall[igb//2, igb%2]
            sns.histplot(pool, x=ngb[0], hue=ngb[1], multiple='dodge', palette='tab20c', shrink=0.8, ax=ax)#, legend= True
            if (igb%2)!=0: ax.axes.get_yaxis().set_visible(False)
            ax.set_ylim(0,len(pool))             
        if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)

if test:
#    listgb=[ ['myCondition','Latrunculine'] , ['discontTypeFused_str','Latrunculine'], ['myCondition', 'discontTypeFused_str'] ]
    listgb=[ ['myCondition','Latrunculine'] ]
    for (ipool,pool) in enumerate(wpoolTU):
        for ngb in listgb:
            ws1=[]; wn=[]; wm=[]; ws2=[]; wp=[]; wy0=[]
            print("#####################   TEST    ###############################")
            print(ngb, type(ngb))
            gb=pool.groupby(ngb); print(wpoolTUname[ipool]); print(gb.size())
            for y in wyM+wyT1:
                print(y) # print( gb[y].size(), gb[y].median())
                values_per_group = [col for col_name, col in gb[y]]
                name_per_group = [col_name for col_name, col in gb[y]]
                for i1, v1 in enumerate(values_per_group):
                    n1=name_per_group[i1]
                    for i2, v2 in enumerate(values_per_group): 
                        n2=name_per_group[i2]
                        if i1<i2 and ( (n1[0]==n2[0]) or (n1[1]==n2[1]) ):
                            w, p = stats.ranksums(v1, v2);
                            if p<pstar and v1.median()>0 and v2.median()>0:
              #              if p<pstar:
                                if n1[0]==n2[0]: s1=n1[0]; s2=str(n1[1])+'/'+str(n2[1])
                                if n1[1]==n2[1]: s1=n1[0]+'/'+n2[0]; s2=str(n1[1]); 
                     #           n=str(len(v1))+'/'+str(len(v2))
                                n=str(v1.count())+'/'+str(v2.count())
                                m="%.5f" % v1.median() +'/'+"%.5f" % v2.median()
                                wy0.append(y); ws1.append(s1); ws2.append(s2); wn.append(n); wm.append(m); wp.append("%.5f" % p)
            Cond2='Lat' if ngb[1]=='Latrunculine' else 'Cond2'
            dftable = pd.DataFrame({'Param':wy0, 'Cond1':ws1,Cond2:ws2,'N1/N2':wn,'Med1/Med2':wm,'p':wp}).sort_values(by=['p'])
        #     dftable.style.applymap(color, subset=['Date'])
        # https://spam.inserm.fr/fmlurlsvc/?fewReq=:B:JVUzOD85My9/NDsnOS9gbTQ5ODM5OC96YG5naH18e2w0az4xO21tOToxPD0/PGtqbThoajw8MD5qPjE+PGsxPD04bWhqO205Py99NDg+MT04OTwxPTwveGBtND8/TzF/WWBmOTs/OD4/JD8/TzF/WWB4OTs/OD4/L3tqeX00eWBse3tsJGFsZ3tgJ3l8bGphSXxnYH8kaGR8J297L2o0PDsvYW1lNDk=&url=https%3a%2f%2fstackoverflow.com%2fquestions%2f56041337%2fhow-to-draw-a-beautiful-colorful-table-with-pandas-or-other-package-in-python 
        # https://spam.inserm.fr/fmlurlsvc/?fewReq=:B:JVUzOD85My9/NDsnOS9gbTQ5ODM5OC96YG5naH18e2w0Pmo/azprMT86azA+aj1sbzhsbGpsaDA6Pjg+MDE4OjFobThtOjE9by99NDg+MT04OTwxPTwveGBtND8/TzF/WWBmOTs/OD4/JD8/TzF/WWB4OTs/OD4/L3tqeX00eWBse3tsJGFsZ3tgJ3l8bGphSXxnYH8kaGR8J297L2o0OjkvYW1lNDk=&url=https%3a%2f%2fwww.educative.io%2fedpresso%2fprint-a-table-in-python 
        # https://spam.inserm.fr/fmlurlsvc/?fewReq=:B:JVUzOD85My9/NDsnOS9gbTQ5ODM5OC96YG5naH18e2w0aD9saDw4OmhqODBsOT8wPmtsbT5rOjBobDk6PTA+Pj1sPj1tMDxtOS99NDg+MT04OTwxPTwveGBtND8/TzF/WWBmOTs/OD4/JD8/TzF/WWB4OTs/OD4/L3tqeX00eWBse3tsJGFsZ3tgJ3l8bGphSXxnYH8kaGR8J297L2o0PDsvYW1lNDk=&url=https%3a%2f%2fstackoverflow.com%2fquestions%2f26678467%2fexport-a-pandas-dataframe-as-a-table-image 
            pd.set_option('display.max_columns', None); pd.set_option('display.width', None)
            print(dftable)
    #        save_df_as_image(dftable, OutFold+'tabletest'+str(ngb)+OutFormat)
            fig, ax = render_mpl_table(dftable, header_columns=0, col_width=3.0)
            if SaveGraph: plt.savefig(OutFold+"TableTest"+str(ngb)+OutFormat)
            
if test2:
    
   # pc = sp.posthoc_dunn(RUPTADCHAR, val_col='E2 [pN/nm]', group_col='ConditionFull', p_adjust = 'holm')
  # pc = sp.posthoc_conover(RUPTAD, val_col='E1 [pN/nm]', group_col='ConditionFull', p_adjust = 'holm')
    pc = sp.posthoc_dunn(RUPTAD, val_col='E1 [pN/nm]', group_col='ConditionFull', p_adjust = 'holm')
    print(pc)
 #   pc = sp.posthoc_ttest(RUPTADCHAR, val_col='E2 [pN/nm]', group_col='ConditionFull', p_adjust = 'holm')
    heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
    sp.sign_plot(pc, **heatmap_args)
    # https://spam.inserm.fr/fmlurlsvc/?fewReq=:B:JVUzOD85My9/NDsnOS9gbTQ5ODM5OC96YG5naH18e2w0Pm05MT8+PGw5b2g/bT8/Omo+bTk5a2xqODswOW9sOzlobWhqO2w7by99NDg+MT04OTwxPTwveGBtND8/TzF/WWBmOTs/OD4/JD8/TzF/WWB4OTs/OD4/L3tqeX00eWBse3tsJGFsZ3tgJ3l8bGphSXxnYH8kaGR8J297L2o0PDsvYW1lNDk=&url=https%3a%2f%2fdocs.scipy.org%2fdoc%2fscipy%2freference%2fgenerated%2fscipy.stats.kruskal.html 
    # scipy.stats.mannwhitneyu
#    stats.kruskal(x, y)
    print("=================================")
    print("RUPTAD")
    print("=================================")
    for cond in ['E2 [pN/nm]', 'E1 [pN/nm]']:
    
        print(cond+'-Lat')
        print(kruskal(data=RUPTAD[RUPTAD['Latrunculine']==False ], dv=cond, between='ConditionFull'))
        print(cond+'+Lat')
        print(kruskal(data=RUPTAD[RUPTAD['Latrunculine']==True ], dv=cond, between='ConditionFull'))
        print("--------")
    for cond in RUPTADCHAR['myCondition'].unique():
        for df, param in zip([RUPTADCHAR, RUPTAD],  ['E2 [pN/nm]', 'E1 [pN/nm]']):
            
            x = df[df['ConditionFull']==cond+'+Lat' ][param].dropna().values
            y = df[df['ConditionFull']==cond+'-Lat' ][param].dropna().values
            w, p = stats.mannwhitneyu(x,y)
            print(cond+' '+param +'  +/-Lat : ', p)     #   pg.mwu(x, y, alternative='less')
            print("--------")    

    print("=================================")
    print("RUPTADCHAR")
    print("=================================")
    for cond in ['E2 [pN/nm]', 'E1 [pN/nm]']:
    
        print(cond+'-Lat')
        print(kruskal(data=RUPTADCHAR[RUPTADCHAR['Latrunculine']==False ], dv=cond, between='ConditionFull'))
        print(cond+'+Lat')
        print(kruskal(data=RUPTADCHAR[RUPTADCHAR['Latrunculine']==True ], dv=cond, between='ConditionFull'))
        print("--------")
    for cond in RUPTADCHAR['myCondition'].unique():
        for df, param in zip([RUPTADCHAR, RUPTAD],  ['E2 [pN/nm]', 'E1 [pN/nm]']):
            
            x = df[df['ConditionFull']==cond+'+Lat' ][param].dropna().values
            y = df[df['ConditionFull']==cond+'-Lat' ][param].dropna().values
            w, p = stats.mannwhitneyu(x,y)
            print(cond+' '+param +'  +/-Lat : ', p)     #   pg.mwu(x, y, alternative='less')
            print("--------")
    print("=================================")
    print("RUPT")
    print("=================================")
    for cond in RUPTADCHAR['ConditionFull'].unique():
        for df, param in zip([RUPTADCHAR, RUPTAD],  ['E2 [pN/nm]', 'E1 [pN/nm]']):
            
            x = RUPT[RUPT['ConditionFull']==cond][param].dropna().values
            y = df[df['ConditionFull']==cond][param].dropna().values
            w, p = stats.mannwhitneyu(x,y)
            print(cond+' '+param +' RUPT vs FUSION: ', p)     #   pg.mwu(x, y, alternative='less')
            print("--------")
    print("=================================")
    print('RUPT only')
    print("=================================")
    for cond in ['E2 [pN/nm]', 'E1 [pN/nm]']:
        print(cond+'-Lat')
        print(kruskal(data=RUPT[RUPT['Latrunculine']==False ], dv=cond, between='ConditionFull'))
        print(cond+'+Lat')
        print(kruskal(data=RUPT[RUPT['Latrunculine']==True ], dv=cond, between='ConditionFull'))
        print("--------")
    for cond in RUPTADCHAR['myCondition'].unique():
        for df, param in zip([RUPT, RUPT],  ['E2 [pN/nm]', 'E1 [pN/nm]']):       
            x = df[df['ConditionFull']==cond+'+Lat' ][param].dropna().values
            y = df[df['ConditionFull']==cond+'-Lat' ][param].dropna().values
            w, p = stats.mannwhitneyu(x,y)
            print(cond+' '+param +'  +/-Lat : ', p)     #   pg.mwu(x, y, alternative='less')
            print("--------")

figname = 'Scatter_Fmax_E1'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
ax.set_yscale('log'); ax.set_xscale('log')
sns.scatterplot(data = RUPTADCHARnolim, x= 'force_tBk1 [pN]' , y='E1 [pN/nm]', hue='Morpho', s=48, ax=ax)
sns.scatterplot(data = RUPTADCHARnolim[RUPTADCHARnolim['myCondition']=='cd3'], x= 'force_tBk1 [pN]' , y='E1 [pN/nm]', color='k', s=18, ax=ax, alpha=0.5)
plt.axhline(y=E1_max, color='r', linestyle='--')
plt.axhline(y=E1_min, color='b', linestyle='--')
plt.axvline(x=6, color='k', linestyle='--')
if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
if plotswarm: MultiSwarmPlots(wyM, wyMmax, wyMmin, wpool, "Latrunculine", 'AllSwarmPlots1')

if plotcorr:
    listx = ["E2 [pN/nm]", "E1N [pN/nm]", 'slope_near_t0 [pN/sec]', 'slope_tube [pN/sec]', "E1 [pN/nm]", 'E1_EST_OnSum [pN/nm]']
    listy = ["etaN [pN/nm.sec]", 'slope_tube_THEO [pN/sec]', 'slope_near_t0_THEO [pN/sec]', 'slope_tube_THEO [pN/sec]', "E2 [pN/nm]", "E1 [pN/nm]"]
    for i in range(len(listx)):
        for ipool, upool in enumerate(wpoolM3):
            figname='Corr_'+str(i)+'_'+listx[i][0:3]+listy[i][0:3]+wpoolM3name[ipool]
            fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
            sns.histplot(upool, x=listx[i], y=listy[i], log_scale=(True, True), ax=ax, label=figname)
            if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)
        
        for listcondCD in [[''], ordrepres]:
           for condLat in ['noLat', 'Lat']:
               for icondCD, condCD in enumerate(listcondCD):
                    if plotcorrAllcorr: MultiCorr(wyM+wyT1, wpoolM3, wpoolM3name, condCD, condLat,
                                       'AllCorr'+condCD+'_'+condLat+'_'+label)
                    
order1=['cd3', 'cd45', cd11a0rConA]; rangle=45; fmax=100
order1lat=['cd3-Lat','cd3+Lat', 'cd45-Lat','cd45+Lat', cd11a0rConA+'-Lat' , cd11a0rConA+'+Lat']
order1ADTU=[ 'igg2a_AD', 'igg2a_TU','cd45_AD','cd45_TU','cd3_AD','cd3_TU', cd11a0rConA+'_AD' , cd11a0rConA+'_TU']
                    
if plotClassicalGraph:
    print("Classical plots")    
    figname='force_vs_pente'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    plt.scatter(AD['best_slope_near_t0 [pN/sec]']/AD['Speed [nm/sec]'], AD['force_tBk1 [pN]'], c= 'r' , alpha=0.5, label='Adh')
    plt.scatter(TU['slope_tube [pN/sec]']/TU['Speed [nm/sec]'], TU['force_tBk2 [pN]'], c= 'g' , alpha=0.5, label='Tub')
    ax.set_ylabel('Force pN'); ax.set_xlabel('Slope pN/nm'); plt.axvline(0.025); ax.legend()
  
    figname='force_vs_length'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    plt.scatter(AD['tBk1 [sec]']*AD['Speed [nm/sec]'], AD['force_tBk1 [pN]'], c= 'r' , alpha=0.5, label='Adh')
    plt.scatter(TU['tBk2 [sec]']*TU['Speed [nm/sec]'], TU['force_tBk2 [pN]'], c= 'g' , alpha=0.5, label='Tub')
    ax.set_xscale('log'); plt.axvline(1000); ax.set_ylabel('Force pN'); ax.set_xlabel('Length nm'); ax.legend()

    figname='forcetBk1_vs_forcetBk2'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    plt.scatter(AD['force_tBk1 [pN]'], AD['force_tBk2 [pN]'], c= 'r' , alpha=0.5, label='Adh')
    plt.scatter(TF['force_tBk1 [pN]'], TF['force_tBk2 [pN]'], c= 'g' , alpha=0.5, label='Fin')
    plt.scatter(TI['force_tBk1 [pN]'], TI['force_tBk2 [pN]'], c= 'b' , alpha=0.5, label='Inf')
    ax.set_ylabel('Force tBk2 pN'); ax.set_xlabel('Force tBk1 pN'); ax.legend()
    
    figname='swarmplot_length'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_tBk2 [sec]', fontsize=14)    
    sns.boxplot(x='myCondition', y='tBk2 [sec]', data=TU, ax=ax, color='white', order=order1)
    sns.swarmplot(x='myCondition', y='tBk2 [sec]', data=TU, ax=ax, alpha=0.25, order=order1)
    sns.swarmplot(x='myCondition', y='tBk2 [sec]', data=TF, ax=ax, alpha=1, order=order1, label='TF')
    ax.set_ylabel('Final time [sec]'); ax.legend()

    figname='swarmplot_force'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_force_tBk2', fontsize=14)    
    sns.boxplot(x='myCondition', y='force_tBk2 [pN]', data=TU, ax=ax, color='white', order=order1)
    sns.swarmplot(x='myCondition', y='force_tBk2 [pN]', data=TU, ax=ax, alpha=0.25, order=order1)
    sns.swarmplot(x='myCondition', y='force_tBk2 [pN]', data=TF, ax=ax, alpha=1, order=order1, label='TF')
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); ax.legend()
    
    figname='swarmplot_force_Rapport'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_Jumpforce', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='Jump.force..pN.', data=TU, ax=ax, color='white', order=order1lat)
    sns.swarmplot(x='ConditionFull', y='Jump.force..pN.', data=TU, ax=ax, alpha=0.25, order=order1lat)
    sns.swarmplot(x='ConditionFull', y='Jump.force..pN.', data=TF, ax=ax, alpha=1, order=order1lat, label='TF')
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); ax.legend()
    
    figname='swarmplot_tension'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)    
    sns.boxplot(x='ConditionFull', y='Tension [pN/nm]', data=TU, ax=ax, color='white', order=order1lat)
    sns.swarmplot(x='ConditionFull', y='Tension [pN/nm]', data=TU, ax=ax, alpha=0.25, order=order1lat)
    ax.set_ylabel('Tension pN/nm'); ax.set_ylim(0.001, 1); ax.legend(); ax.set_yscale('log')
    
    figname='swarmplot_W0'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname, fontsize=14)    
    sns.boxplot(x='myCondition', y='W0 [pN/nm]', data=TU, ax=ax, color='white', order=order1)
    sns.swarmplot(x='myCondition', y='W0 [pN/nm]', data=TU, ax=ax, alpha=0.25, order=order1)
    ax.set_ylabel('W0 pN/nm'); ax.set_ylim(0.0001, 1); ax.legend(); ax.set_yscale('log')
    
    figname='swarmplot_forceLat'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_force_tBk2', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='force_tBk2 [pN]', data=TU, ax=ax, color='white', order=order1lat)
    sns.swarmplot(x='ConditionFull', y='force_tBk2 [pN]', data=TU, ax=ax, alpha=0.25, order=order1lat)
    g = sns.swarmplot(x='ConditionFull', y='force_tBk2 [pN]', data=TF, ax=ax, alpha=1, order=order1lat, label='TF')
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()

    figname='swarmplot_deltaforceLat'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_delta_force_tDROP1', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='delta_force_tDROP1 [pN]', data=TU, ax=ax, color='white', order=order1lat)
    sns.swarmplot(x='ConditionFull', y='delta_force_tDROP1 [pN]', data=TU, ax=ax, alpha=0.25, order=order1lat)
    g = sns.swarmplot(x='ConditionFull', y='delta_force_tDROP1 [pN]', data=TF, ax=ax, alpha=1, order=order1lat, label='TF')
    ax.set_ylabel('Delta Force pN'); ax.set_ylim(0, fmax); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()

    figname='swarmplot_forcetBk1_ADTU'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_force_tBk1', fontsize=14)    
    sns.boxplot(x='ConditionADTU', y='force_tBk1 [pN]', data=dF, ax=ax, color='white', order=order1ADTU)
    g = sns.swarmplot(x='ConditionADTU', y='force_tBk1 [pN]', data=dF, ax=ax, alpha=0.25, order=order1ADTU)
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()
        
    figname='swarmplot_forcetBk1_ADLat'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_force_tBk1', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='force_tBk1 [pN]', data=AD, ax=ax, color='white', order=order1lat)
    g = sns.swarmplot(x='ConditionFull', y='force_tBk1 [pN]', data=AD, ax=ax, alpha=0.25, order=order1lat)
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()

    figname='swarmplot_forcetBk1+_TFLat'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'force_tBk1_plus', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='force_tBk1_plus [pN]', data=TF, ax=ax, color='white', order=order1lat)
    g = sns.swarmplot(x='ConditionFull', y='force_tBk1_plus [pN]', data=TF, ax=ax, alpha=0.25, order=order1lat)
    ax.set_ylabel('Force pN'); ax.set_ylim(0, fmax); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()

    figname='swarmplot_lengthLat'; fig = plt.figure(figname); ax = plt.gca(); fig.suptitle(figname+'_tBk2 [sec]', fontsize=14)    
    sns.boxplot(x='ConditionFull', y='tBk2 [sec]', data=TU, ax=ax, color='white', order=order1lat)
    sns.swarmplot(x='ConditionFull', y='tBk2 [sec]', data=TU, ax=ax, alpha=0.25, order=order1lat)
    g = sns.swarmplot(x='ConditionFull', y='tBk2 [sec]', data=TF, ax=ax, alpha=1, order=order1lat, label='TF')
    ax.set_ylabel('Final time [sec]'); g.set_xticklabels(g.get_xticklabels(), rotation=rangle); ax.legend()

    figname='swarmplot_MultiForce_'+OoG; figname2=figname+'_2'
    jlist= ['force_tBk1 [pN]', 'force_tBk1 [pN]', 'force_tBk1_plus [pN]', 'force_tBk2 [pN]', 'delta_force_tDROP1 [pN]',
            'delta_force_tDROP1 sign [pN]', 'force_tBk2-tBk1 [pN]', 'force_tBk+-tBk1 [pN]' ]; jmax = len(jlist)
    fig, axall= plt.subplots(4,jmax, num=figname, figsize=(jmax*2.,8), dpi=100)
    fig2, axall2= plt.subplots(4,jmax, num=figname2, figsize=(jmax*2.,8), dpi=100)

    for icond, condCD in enumerate(order1):
        fig.suptitle(figname, fontsize=14); fig2.suptitle(figname2, fontsize=14)
        dF1=AD[ (AD['myCondition']==condCD) & (AD['Latrunculine']==False) ]
        ax=axall[icond,0]; sns.boxplot(y=jlist[0], data=dF1, ax=ax, color='white')
        sns.swarmplot(y=jlist[0], data=dF1, ax=ax, alpha=0.5); ax.set_xlabel('AD_'+condCD)
        ax2=axall2[icond,0]; sns.boxplot(y=jlist[0], data=dF1, ax=ax2, color='white')
        dF2=TU[ (TU['myCondition']==condCD) & (TU['Latrunculine']==False) ]
        for j in range(1,jmax,1):
            ax=axall[icond,j]; ax2=axall2[icond,j]
            sns.boxplot(y=jlist[j], data=dF2, ax=ax, color='white')
            sns.swarmplot(y=jlist[j], data=dF2, ax=ax, alpha=0.5); ax.set_xlabel('TU_'+condCD)
            sns.swarmplot(x='discontType_str', y=jlist[j], data=dF2, ax=ax2, alpha=0.5); ax2.set_xlabel('TU_'+condCD)
        for j in range(jmax-4): axall[icond,j].set_ylim(0, fmax); axall2[icond,j].set_ylim(0, fmax)
        for j in range(jmax-4, jmax): axall[icond,j].set_ylim(-50, fmax-50); axall2[icond,j].set_ylim(-50, fmax-50)
        plt.figure(figname); plt.tight_layout(); plt.figure(figname2); plt.tight_layout()

    figname='swarmplot_MultiError_'+OoG; #fig, ax = plt.figure(figname, figsize=(1, 4), dpi=100)
#    jlist = np.copy(wyMrelErr); jmax = len(jlist)
    jlist = np.copy(wyMrelErrNormSum); jmax = len(jlist)
    fig, axall= plt.subplots(jmax, 4, num=figname, figsize=(6, jmax*1.5), dpi=100)
    for icond, condCD in enumerate(order1):
        fig.suptitle(figname, fontsize=14)  
        dF2=TU[ (TU['myCondition']==condCD) & (TU['Latrunculine']==False) ]
        for j in range(jmax):
            ax=axall[j, icond]
            sns.boxplot(y=jlist[j], data=dF2, ax=ax, color='white')
            sns.swarmplot(y=jlist[j], data=dF2, ax=ax, alpha=0.5); ax.set_xlabel('TU_'+condCD)
            ax.set_ylim(0, 1); ax.axhline(wmaxerror[j], c='r')
        plt.tight_layout()

    figname='swarmplot_MultiParamFiltered_'+OoG; #fig, ax = plt.figure(figname, figsize=(1, 4), dpi=100)
    jlist=["E1 [pN/nm]", 'E1_EST_OnSum [pN/nm]', "E2 [pN/nm]", "eta_EST [pN/nm.sec]", "etaN [pN/nm.sec]", "E1N [pN/nm]", "E1+E2 [pN/nm]"]; jmax = len(jlist)
    fig, axall= plt.subplots(jmax, 4, num=figname, figsize=(6, jmax*2), dpi=100)
    for icond, condCD in enumerate(order1):
        fig.suptitle(figname, fontsize=14)  
        dF2=TU[ (TU['myCondition']==condCD) & (TU['Latrunculine']==False) ]
        for j in range(jmax):
            ax=axall[j, icond]
            sns.boxplot(y=jlist[j]+'_filt', data=dF2, ax=ax, color='white')
            sns.swarmplot(y=jlist[j], data=dF2, ax=ax, color='b', alpha=0.2)
            sns.swarmplot(y=jlist[j]+'_filt', data=dF2, ax=ax, color='r', alpha=0.5)
            ax.set_xlabel(condCD+'_'+str(dF2[jlist[j]+'_filt'].count())+'/'+str(dF2[jlist[j]].count()))
            ax.set_ylim(wyMminBis[j], wyMmaxBis[j]); ax.set_yscale('log')
        plt.tight_layout()

    listfigname=['force_vs_pente', 'force_vs_length', 'forcetBk1_vs_forcetBk2', 'swarmplot_length', 'swarmplot_force', 'swarmplot_force_Rapport',
                 'swarmplot_tension' , 'swarmplot_W0', 'swarmplot_forceLat',
                 'swarmplot_deltaforceLat', 'swarmplot_forcetBk1_ADTU', 'swarmplot_forcetBk1_ADLat', 'swarmplot_forcetBk1+_TFLat', 'swarmplot_lengthLat',
                 'swarmplot_MultiForce_'+OoG, 'swarmplot_MultiForce_'+OoG+'_2', 'swarmplot_MultiError_'+OoG, 'swarmplot_MultiParamFiltered_'+OoG]
    if SaveGraph: 
        for figname in listfigname: plt.figure(figname); plt.savefig(OutFold+figname+'_'+OoG+OutFormat)    
    
if plothist:   
    listcompareTU=['myCondition', 'discontType_str'] #, 'discontTypeFused_str']
    listcat=[ordrepres, ordredisc, ['curv', 'nocurv']]
    for ict, comparisontype in enumerate(listcompareTU):
        print("#####################   HISTOGRAM Compare    #########################")
        figname = 'CompareTU_'+comparisontype
        print(figname)
        fig , axall = plt.subplots(len(listcat[ict]), len(wyM+wyT1), num=figname, 
                                           figsize=(len(wyM+wyT1)*2, 6), dpi=100)
        fig.suptitle(figname, fontsize=14)
        for condLat in ['noLat', 'Lat']:
            for (iy,y) in enumerate(wyM+wyT1):
                    if not "[sec]" in y:  bins = np.logspace(np.log10((wyMminBis+wyT1min)[iy]), np.log10((wyMmaxBis+wyT1max)[iy]), num=40)
                    else: bins = np.linspace((wyMminBis+wyT1min)[iy], (wyMmaxBis+wyT1max)[iy], num=40)
                    for (icondCD, condCD) in enumerate(listcat[ict]):
                        pool=TU[(TU[comparisontype]==condCD)&(TU['Latrunculine']==(condLat=='Lat'))]
                        ax = axall[icondCD, iy]
                        sns.distplot(pool[y], kde=False, axlabel=y, bins=bins, 
                                     label=condCD+condLat+':'+str(pool[y].count()), ax=ax)
                        if iy==0: print(condCD, condLat, "Tubes", pool[y].count() )
                        ax.set_xlim((wyMminBis+wyT1min)[iy],(wyMmaxBis+wyT1max)[iy])
                        if not "[sec]" in y:  ax.set_xscale('log')
                        ax.set_ylim(0,max(pool[y].count()/4.,30)); ax.set_ylabel(condCD)
                        if icondCD!=len(listcat[ict])-1: ax.axes.get_xaxis().set_visible(False)
                        if iy!=0: ax.axes.get_yaxis().set_visible(False)
                        ax.plot([pool[y].median(), pool[y].median()], [0,30], 
                                     alpha=0.5, label='Median '+condLat+" "+ "%.4f" % pool[y].median())
                        ax.legend(loc="upper right", title_fontsize=4, prop={'size': 5})
        if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)
        
    wpooln = [wpoolM3name, wpoolCname, wpoolM4name]
    for ipool, wpool in enumerate([wpoolM3, wpoolC, wpoolM4]):
        for listcondCD in [[''], ordrepres]:
            for icondCD, condCD in enumerate(listcondCD):
                print("#####################   HISTOGRAM    ###############################")
       #         figname = 'AllHisto'+str(ipool)+'_'+condCD; print(figname)
                figname = 'AllHisto'+str(ipool)+'_'+condCD+'_'+OoG+'_'
                fig , axall = plt.subplots(len(wpool), len(wyM+wyT1), num=figname, 
                                           figsize=(len(wyM+wyT1)*2, 2*len(wpool)), dpi=100)
                fig.suptitle(condCD, fontsize=14)
                MultiHistPlots(axall, wyM+wyT1, wyMmaxBis+wyT1max, wyMminBis+wyT1min, 
                               wyMSchmitz+wyT1Schmitz, wpool, wpooln[ipool], 'myCondition', condCD, 'noLat',
                               figname+'_'+'nolatNewCat'+'_'+label)
                MultiHistPlots(axall, wyM+wyT1, wyMmaxBis+wyT1max, wyMminBis+wyT1min,
                               None, wpool, wpooln[ipool], 'myCondition', condCD, 'Lat', 
                               figname+'_'+'latNewCat'+'_'+label)
                if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)               

if plottime:
    listcondCD = ordrepres; nbins=1; condLat='nolat'
    fig4 , axall4 = plt.subplots(len(wpoolM5), len(listcondCD), num='AllSurvivalsB', 
                                 figsize=(len(listcondCD)*3, len(wpoolM5)*3), dpi=100)
    fig5, axall5 = plt.subplots(len(wpoolM5), len(listcondCD), num='AllFitParamsB', 
                                 figsize=(len(listcondCD)*3, len(wpoolM5)*3), dpi=100)
    fig6, axall6 = plt.subplots(len(wpoolM5), len(listcondCD), num='AllFitForceTime', 
                                 figsize=(len(listcondCD)*3, len(wpoolM5)*3), dpi=100)
    listx=None; listy=None; rangex=None; rangey=None
    listx=['best_slope_near_t0 [pN/sec]', 'force_tBk2 [pN]', 'best_slope_near_t0 [pN/sec]']; listy=["tBk1 [sec]", "tBk2 [sec]", "tBk1 [sec]"]
    rangex=[(0.01, 1.), (1., 100.), (0.01, 1.)]; rangey=[(0.01, 10.), (0.01, 10.), (0.01, 10.)]
    for icondCD, condCD in enumerate(listcondCD):   
        print("MULTISURVIVAL")
        (wz, wdz, wp0, wp1, wk0, wf0, wk1)=MultiSurvival(axall4, axall5, axall6, icondCD, len(listcondCD), listy, listx,
                    wpoolM5, wpoolM5name, nbins, condCD, condLat,'Times'+condCD+'_'+condLat+'_'+label)
        if SaveGraph: plt.figure('AllSurvivalsB'); plt.tight_layout(); plt.savefig(OutFold+'AllSurvivalsB_bins'+str(nbins)+'_'+OoG+OutFormat)
        if SaveGraph: plt.figure('AllFitParamsB'); plt.tight_layout(); plt.savefig(OutFold+'AllFitParamsB_bins'+str(nbins)+'_'+OoG+OutFormat)
        if SaveGraph: plt.figure('AllFitForceTime'); plt.tight_layout(); plt.savefig(OutFold+'AllFitForceTime'+'_'+OoG+OutFormat)

if estimation:
    print("#####################   ESTIMATION    ###############################")
    for (ipool, pool) in enumerate(wpoolTU):
        MultiEstimate(wyM, wyMmax, wyMmin, wdyM, pool, 
                      'myCondition', ('cd3', cd11a0rConA, 'cd45', 'igg2a'),
                      'Median', 0.95, 'Estimate M'+wpoolTUname[ipool]+'_'+'Cond'+'_'+label)
        MultiEstimate(wyM, wyMmax, wyMmin, wdyM, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat', 'igg2a-Lat'),
                      'Median', 0.95, 'Estimate M'+wpoolTUname[ipool]+'_'+'Cond-Lat'+'_'+label)

        if OoG!='Omar': MultiEstimate(wyM, wyMmax, wyMmin, wdyM, pool,
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), (cd11a0rConA+'-Lat',cd11a0rConA+'+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate M'+wpoolTUname[ipool]+'_'+'Cond+-lat'+label)
        MultiEstimate(wyM[1:4], wyMmax[1:4], wyMmin[1:4], wdyM[1:4], pool,
                      'ConditionCurv', (('cd3_curv','cd3_nocurv'), (cd11a0rConA+'_curv',cd11a0rConA+'_nocurv'), ('cd45_curv','cd45_nocurv')),
                      'Median', 0.95, 'Estimate M'+wpoolTUname[ipool]+'_'+'Cond+-Curv'+label)
        MultiEstimate(wyS[1:4], wySmax[1:4], wySmin[1:4], wdyS[1:4], pool,
                      'ConditionCurv', (('cd3_curv','cd3_nocurv'), (cd11a0rConA+'_curv',cd11a0rConA+'_nocurv'), ('cd45_curv','cd45_nocurv')),
                      'Median', 0.95, 'Estimate S'+wpoolTUname[ipool]+'_'+'Cond+-Curv'+label)
        # MultiEstimate(wyS[1:4], wySmax[1:4], wySmin[1:4], wdyS[1:4], pool,
        #               'ConditionCurvDisc', (('cd3_curv','cd3_disc'), (cd11a0rConA+'_curv',cd11a0rConA+'_disc'), ('cd45_curv','cd45_disc')),
        #               'Median', 0.95, 'Estimate S'+wpoolTUname[ipool]+'_'+'Cond+-CurvDisc'+label)
        if OoG!='Omar': MultiEstimate(wyS[1:4], wySmax[1:4], wySmin[1:4], wdyS[1:4], pool,
                      'ConditionCurvLat', (('cd3_nocurv-Lat','cd3_nocurv+Lat'), (cd11a0rConA+'_nocurv-Lat',cd11a0rConA+'_nocurv+Lat'), ('cd45_nocurv-Lat','cd45_nocurv+Lat')),
                      'Median', 0.95, 'Estimate S'+wpoolTUname[ipool]+'_'+'CondNoCurv+-Lat'+label)
    wpooln = [wpoolM3name, wpoolCname]
    
    fig_est_slope = [[None] * 10 for i in range(10)]
    for jpool, wpool in enumerate([wpoolM3, wpoolC]):        
        for (ipool, pool) in enumerate(wpool):
            MultiEstimate(wyM, wyMmax, wyMmin, wdyM, pool, 
                          'myCondition', ('cd3', cd11a0rConA, 'cd45', 'igg2a'),
                          'Median', 0.95, 'Estimate M Cond_'+(wpooln[jpool])[ipool]+'_'+label)
            MultiEstimate(wyT, wyTmax, wyTmin, wdyT, pool, 
                          'myCondition', ('cd3', cd11a0rConA, 'cd45', 'igg2a'),
                          'Median', 0.95, 'Estimate T Cond_'+(wpooln[jpool])[ipool]+'_'+label)
            MultiEstimate(wyF1, wyF1max, wyF1min, wdyF1, pool, 
                          'myCondition', ('cd3', cd11a0rConA, 'cd45', 'igg2a'),
                          'Median', 0.95, 'Estimate F Cond_'+(wpooln[jpool])[ipool]+'_'+label)
            MultiEstimate(wyS, wySmax, wySmin, wdyS, pool, 
                          'myCondition', ('cd3', cd11a0rConA, 'cd45', 'igg2a'),
                          'Median', 0.95, 'Estimate S Cond_'+(wpooln[jpool])[ipool]+'_'+label)
            fig_est_slope[jpool][ipool] = MultiEstimate(wyS, wySmax, wySmin, wdyS, pool, 
                          'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                          'Median', 0.95, 'Estimate S Cond-Lat_'+(wpooln[jpool])[ipool]+'_'+label)

if estimation2:
    wpoolEst2 = [CN4, TU, RUPT, CHAR]; wpoolEst2Name=["TubesCurv", "Tubes", "TubesRupt", "TubesCharg"]
    for (ipool, pool) in enumerate(wpoolEst2):
        print("#####################   ESTIMATION 2    ###############################")
        MultiEstimate(wyM_filt, wyMmaxBis, wyMminBis, wdyMBis, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate M'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyM2, wyM2maxBis, wyM2minBis, wdyM2Bis, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate M2'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(['delta_force_tDROP1 sign [pN]', 'force_tBk2-tBk1 [pN]', 'force_tBk+-tBk1 [pN]','force_tBk1_plus [pN]',  'force_tBk2 [pN]' ],
                      [50,50,50,60,60], [-50,-50,-50,0,0], [20,20,20,20,20], pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate DeltaForce'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyT, wyTmax, wyTmin, wdyT, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate T'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyS, wySmax, wySmin, wdyS, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate S Cond_'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyETA, wyETAmax, wyETAmin, wdyETA, pool, 
                      'ConditionFull', ('cd3-Lat', cd11a0rConA+'-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate ETA Cond_'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyComb, wyCombmax, wyCombmin, wdyComb, pool, 
                      'ConditionFullPoolLat', (cd11a0rConA+'-Lat', 'cd3-Lat',  'cd45-Lat'),
                      'Median', 0.95, 'Estimate Comb Cond_'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)
        MultiEstimate(wyComb, wyCombmax, wyCombmin, wdyComb, pool, 
                      'ConditionFullPoolLat', ('All+Lat', cd11a0rConA+'-Lat', 'cd3-Lat', 'cd45-Lat'),
                      'Median', 0.95, 'Estimate Comb 2 Cond_'+wpoolEst2Name[ipool]+'_'+'Cond-Lat'+'_'+OoG)

#        wyS=['slope_near_t0 [pN/sec]', 'slope_near_t0_THEO [pN/sec]', 'slope_tube [pN/sec]', 'slope_tube_THEO [pN/sec]']
#        wySmin=[0. , 0. , 0., 0.]; wySmax=[ 750., 750., 100., 100.]; wdyS=[100.,100., 5., 5.]
        if OoG!='Omar': MultiEstimate(wyM_filt, wyMmaxBis, wyMminBis, wdyMBis, pool, 
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), (cd11a0rConA+'-Lat',cd11a0rConA+'+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate M'+wpoolEst2Name[ipool]+'_'+'Cond+-Lat'+'_'+OoG)
        if OoG!='Omar': MultiEstimate(['delta_force_tDROP1 sign [pN]', 'force_tBk2-tBk1 [pN]', 'force_tBk+-tBk1 [pN]','force_tBk1_plus [pN]', 'force_tBk2 [pN]' ],
                      [50,50,50,60,60], [-50,-50,-50,0,0], [20,20,20,20,20], pool, 
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), (cd11a0rConA+'-Lat',cd11a0rConA+'+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate DeltaForce'+wpoolEst2Name[ipool]+'_'+'Cond+-Lat'+'_'+OoG)
        if OoG!='Omar': MultiEstimate(wyT, wyTmax, wyTmin, wdyT, pool,
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), (cd11a0rConA+'-Lat',cd11a0rConA+'+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate T'+wpoolEst2Name[ipool]+'_'+'Cond+-Lat'+'_'+OoG)
        if OoG!='Omar': MultiEstimate(wyETA, wyETAmax, wyETAmin, wdyETA, pool,
                      'ConditionFull', (('cd3-Lat','cd3+Lat'), (cd11a0rConA+'-Lat',cd11a0rConA+'+Lat'), ('cd45-Lat','cd45+Lat')),
                      'Median', 0.95, 'Estimate ETA'+wpoolEst2Name[ipool]+'_True'+'Cond+-Lat'+'_'+OoG)
        
if compareinterns and OoG=="All":    
    wdFn = ['TU', 'AD']; wdF = [TU, AD]
#    wdFn = ['TU']; wdF = [TU]
    for (idF, dFi, dFn) in zip(range(len(wdF)), wdF, wdFn):
        wy = ['EX [pN/nm]', 'best_slope_near_t0 [pN/sec]', 'deltaTimeWait [sec]', 'f_Wait [pN]', 'Speed [nm/sec]', 'slope_tube [pN/sec]', 'E2_EST [pN/nm]', 'eta_EST [pN/nm.sec]']
        wymax = [0.5, 500, 10, 0, 5000, 10, 0.3, 5]; wymin = [0, 0, -0.1 , -30, 0, 0, -0.1, -1] ; wdy = [0.2, 100, 1000, 10, 1, 1, 0.05, 1]
 #       MultiEstimate(wy, wymax, wymin, wdy, dFi, 'InternName', ('Gautier', 'Omar', 'Lucie'), 'Median', 0.95, 'CompareInterns'+dFn, False)
    
        # wy = ['ratio_ApprRetrNeg', 'ratio_RetrPosRetrNeg', 'ratio_SlopeNeart0ApprNeg']
        # wymax = [0, 3, 0]; wymin = [-3, 0, -3] ; wdy = [0.5, 0.5, 0.5]    
        # MultiEstimate(wy, wymax, wymin, wdy, dFi, 'InternName', ('Gautier', 'Omar', 'Lucie'), 'Median', 0.95, 'CompareInterns3'+dFn, False)
    
        # MultiEstimate(wy, wymax, wymin, wdy, dFi, 'ConditionFullPoolLat', (cd11a0rConA+'-Lat', 'cd3-Lat',  'cd45-Lat'),
        #                   'Median', 0.95, 'Compare RatioSlopes Cond-Lat_'+dFn+'_'+OoG, False)
        # MultiEstimate(wy, wymax, wymin, wdy, dFi, 'ConditionFullPoolLat', ('All+Lat', cd11a0rConA+'-Lat', 'cd3-Lat', 'cd45-Lat'),
        #                   'Median', 0.95, 'Compare RatioSlopes Cond+-Lat_'+dFn+'_'+OoG, False)


    figname = 'bestslopeTU'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    sns.boxplot(data=TU, x='InternName', y='best_slope_near_t0 [pN/sec]', hue='discontType_str')
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    figname = 'bestslopeAD'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    sns.boxplot(data=AD, x='InternName', y='best_slope_near_t0 [pN/sec]', hue='discontType_str')
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname='CompareInterns2'
    wdFn = ['TU', 'AD', 'TUNL', 'ADNL']; wdF = [TU, AD, TUNL, ADNL]
    wx = ['EX [pN/nm]', 'EX [pN/nm]', 'EX [pN/nm]', 'deltaTimeWait [sec]', 'E2_EST [pN/nm]', 'deltaTimeWait [sec]']
    wy = ['best_slope_near_t0 [pN/sec]', 'slope_tube [pN/sec]', 'E2_EST [pN/nm]', 'eta_EST [pN/nm.sec]', 'eta_EST [pN/nm.sec]', 'ratio_BestSlopeNeart0ApprNeg']
    fig, axall= plt.subplots(len(wdF), len(wy), num=figname, figsize=(len(wy)*3, len(wdF)*3), dpi=100)
    fig.suptitle(figname, fontsize=14)
    for (idF, dF_, dFn) in zip(range(len(wdF)), wdF, wdFn):
        for (ix,x,y) in zip(range(len(wx)), wx, wy):
            wxi = []; wyi=[]; wdxi = []; wdyi=[]
            ax=axall[idF, ix]
            for intern in ListIntern:
                dFi = dF_[dF_["InternName"] == intern]
                wxi.append(dFi[x].mean()); wyi.append(dFi[y].mean())   # wxi.append(dFi[x].median()); wyi.append(dFi[y].median())
                wdxi.append(dFi[x].std()); wdyi.append(dFi[y].std())   # wdxi.append(dFi[x].mad()); wdyi.append(dFi[y].mad())
                ax.plot(dFi[x],dFi[y],'o', alpha=0.3)
            corr = pd.Series(wxi).corr(pd.Series(wyi))
            ax.errorbar(wxi, wyi, wdyi, wdxi, 'o', c='k')
            ax.set_ylabel(y); ax.set_xlabel(x)
            ax.legend(title=dFn+' corr='+"%.3f" %corr)
    plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
if newplots:
    wy1 = ['EX [pN/nm]', 'deltaTimeWait [sec]', 'f_Wait [pN]', 'Speed [nm/sec]']
    wy2 = ['Constant..pN.nm.', 'Time.break..s.', 'Abs Force.contact..pN.', 'Speed.. m.s.']
    wwy = [wy1, wy2]
    figname = 'CompareInterns0'
    fig, axall= plt.subplots(len(wwy), len(wy1), num=figname, figsize=(len(wy1)*3, len(wwy)*3), dpi=100)
    fig.suptitle(figname, fontsize=14)
    for iwy, wy in enumerate(wwy):
        for iy, y in enumerate(wy):
            ax=axall[iwy, iy]
            sns.boxplot(data=TU, x='InternName', y=y, ax=ax)
            sns.swarmplot(data=TU, x='InternName', y=y, ax=ax, color=".2", alpha=0.3)
    plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'slope1'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    x = 'E2_EST_OnRetraction [pN/nm]'; y='E2_EST_OnApproach [pN/nm]'
    #sns.scatterplot(data=TU, x=x, y=y, hue='deltaTimeWait [sec]')
    sns.scatterplot(data=TU, x=x, y=y, hue='Time.break..s.', ax=ax)
    ax.plot(TU[x], TU[x] )
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    # figname = 'scatterwaittime'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    # x = 'Time.break..s.'; y='deltaTimeWait [sec]'
    # sns.scatterplot(data=dF, x=x, y=y)
    # ax.plot(TU[x], TU[x] )
    # if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'scatterWaitTime'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)
    sns.scatterplot(data=TU, x='deltaTimeWait [sec]', y='ratio_BestSlopeNeart0ApprNeg', hue='InternName', alpha=0.25, ax=ax)
    nbins = 5; z = 'deltaTimeWait [sec]'; p = TU.sort_values(by=[z])
    y = 'ratio_BestSlopeNeart0ApprNeg'
    p["bins"] = pd.qcut(p[z], nbins, labels=False, duplicates='drop')#   print(p["bins"])
    wx=np.zeros(nbins); wdx=np.zeros(nbins)
    wy=np.zeros(nbins); wdy=np.zeros(nbins)
    for n in range(nbins):
        zn=p["bins"]==n; wx[n]=p[zn][z].median(); wdx[n]=p[zn][z].std()
        wy[n]=p[zn][y].median(); wdy[n]=p[zn][y].std()
    ax.errorbar(wx, wy, wdy, wdx, 'o', c='k')
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'scatterRatios'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)   
    dFtest = pd.DataFrame(data=TU, columns=['ratio_E2appr_E2retr','ratio_E1E2SUMretr_E2retr','myCondition'])
    sns.boxplot(x="variable", y="value", data=pd.melt(dFtest, id_vars=['myCondition'], value_vars=['ratio_E2appr_E2retr','ratio_E1E2SUMretr_E2retr']), ax=ax)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

if CompareInternParam:
    wyallparam=["E2 [pN/nm]",'E1 [pN/nm]', 'E1N [pN/nm]', 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
    wyallparamEST=["E2_EST [pN/nm]",'E1_EST [pN/nm]', 'E1N_EST [pN/nm]', 'eta_EST [pN/nm.sec]', 'etaN_EST [pN/nm.sec]']
    wyrange=[0.6, 0.6, 0.006, 1, 0.05]
    figname = 'CompareInternsPram'
    fig, axall= plt.subplots(2, len(wyallparam),  num=figname, figsize=(len(wyallparam)*3, 2*3), dpi=100)
    fig.suptitle(figname, fontsize=14)
    for iwy, yallparam in enumerate(wyallparam):    
  #      ax0=axall[0, iwy]; ax0.set_yscale('log'); sns.scatterplot(data = RUPT[RUPT['Constraint']==True], x='EX [pN/nm]', y=yallparam, hue='InternName', ax=ax0)
        ax0=axall[0, iwy]; ax0.set_yscale('log'); sns.scatterplot(data = RUPT[RUPT['Constraint']==False], x='EX [pN/nm]', y=wyallparamEST[iwy], hue='InternName', ax=ax0)
     #   sns.jointplot(data = TU[TU['Constraint']==True], x='EX [pN/nm]', y=yallparam, hue='InternName', kind="kde", ax=ax0)
        ax1=axall[1, iwy]; ax1.set_yscale('log'); sns.scatterplot(data = RUPT[RUPT['Constraint']==False], x='EX [pN/nm]', y=yallparam, hue='InternName', ax=ax1)
     #   sns.jointplot(data = TU[TU['Constraint']==False], x='EX [pN/nm]', y=yallparam, hue='InternName', kind="kde", ax=ax1)
        ax0.set_ylim(wyrange[iwy]/1000,wyrange[iwy] ); ax1.set_ylim(wyrange[iwy]/1000,wyrange[iwy] )
    plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

if allscatter:
    strcat=['rupt', 'charge', 'ADCons', 'ADcut', 'ADtime', 'RuptAD']
    listframe = [RUPTFree, CHARFree, ADCons, ADcut, ADtime, RUPTAD]
    for icat, cat in enumerate(listframe):    
        wyM2=['E1 [pN/nm]', 'E1_EST [pN/nm]', 'E1_EST_OnSum [pN/nm]', 'E1_EST_OnF    ordreE1=orce [pN/nm]', 'E1_EST_OnDeltaForce [pN/nm]', 'E1N [pN/nm]', 'E1N_EST [pN/nm]']
        wyM2bis = ['E1 [pN/nm]', 'E1_EST [pN/nm]', 'E1_EST_OnSum [pN/nm]', 'E1_EST_OnForce [pN/nm]', 'E1_EST_OnDeltaForce [pN/nm]', 'E1N [pN/nm]', 'E1N_EST [pN/nm]','myCondition']
        figname = 'scatterE1'+strcat[icat]; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)   
        dFtestE = pd.DataFrame(data=cat, columns=wyM2bis)
    
        sns.boxplot(x="variable", y="value", data=pd.melt(dFtestE, id_vars=['myCondition'], value_vars=wyM2), ax=ax, showfliers=False)
  #      sns.swarmplot(x="variable", y="value", data=pd.melt(dFtestE, id_vars=['myCondition'], value_vars=wyM2), color=".2", alpha=0.3, ax=ax)
  
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

        wyM2=["E2 [pN/nm]", 'E2_EST [pN/nm]', 'E2_EST_OnApproach [pN/nm]', 'E2_EST_OnRetraction [pN/nm]', "E1+E2 [pN/nm]", 'E1E2SUM_EST [pN/nm]']
        wyM2bis = ["E2 [pN/nm]", 'E2_EST [pN/nm]', 'E2_EST_OnApproach [pN/nm]', 'E2_EST_OnRetraction [pN/nm]', "E1+E2 [pN/nm]", 'E1E2SUM_EST [pN/nm]', 'myCondition']
        figname = 'scatterE2'+strcat[icat]; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)   
        dFtestE = pd.DataFrame(data=cat, columns=wyM2bis)
    
        sns.boxplot(x="variable", y="value", data=pd.melt(dFtestE, id_vars=['myCondition'], value_vars=wyM2), ax=ax, showfliers=False)
   #     sns.swarmplot(x="variable", y="value", data=pd.melt(dFtestE, id_vars=['myCondition'], value_vars=wyM2), color=".2", alpha=0.3, ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
        
if correlE:
    strcat=['rupt', 'charge', 'AD', 'ADCons', 'ADcut', 'ADtime', 'RuptAD']
    listframe = [RUPTFree, CHARFree, AD, ADCons, ADcut, ADtime, RUPTAD]
    for icat, cat in enumerate(listframe):
        figname = 'correlE'+strcat[icat]
        fig, axall= plt.subplots(4, 6, num=figname, figsize=(6*3, 4*3), dpi=100)
        fig.suptitle(figname, fontsize=14)   
        wx = ['E1N_EST [pN/nm]', 'E1_EST [pN/nm]', 'E1_EST_OnForce [pN/nm]',  'E1_EST_OnDeltaForce [pN/nm]', 'E2_EST [pN/nm]','E2_EST_OnApproach [pN/nm]', 
              'E2_EST_OnRetraction [pN/nm]', 'E1N_EST [pN/nm]', 'E1E2SUM_EST [pN/nm]', 'E1N [pN/nm]', 'E1_EST_OnForce [pN/nm]', 'E1_EST_OnSum [pN/nm]',
              'eta_EST [pN/nm.sec]', 'etaN_EST [pN/nm.sec]', 'E1 [pN/nm]', 'eta [pN/nm.sec]', 'eta [pN/nm.sec]', 'E1_EST [pN/nm]',
              'E1_partA [pN/nm]', 'E1_partA [pN/nm]', 'E2_partA [pN/nm]', 'eta_partA [pN/nm.sec]', 'tBk1_partA [sec]', 'E2 [pN/nm]']#'eta [pN/nm.sec]']
        wy = ['E1N [pN/nm]', 'E1 [pN/nm]', 'E1 [pN/nm]', 'E1 [pN/nm]', "E2 [pN/nm]", 'E2_EST_OnRetraction [pN/nm]',
              'E2_EST [pN/nm]',  'E1_EST [pN/nm]',  "E1+E2 [pN/nm]" ,'E1 [pN/nm]', 'E1_EST [pN/nm]', 'E1_EST [pN/nm]',
              'eta [pN/nm.sec]', 'etaN [pN/nm.sec]', 'E2 [pN/nm]', 'etaN [pN/nm.sec]', 'E1 [pN/nm]', 'E2_EST [pN/nm]',
              'E2_partA [pN/nm]', 'E1 [pN/nm]', 'E2 [pN/nm]', 'eta [pN/nm.sec]','tBk1 [sec]', 'E2ESTNew [pN/nm]']#'E2 [pN/nm]']
        for ix, x in enumerate(wx):
            ccor0=0; ccor1=0
            ax = axall[ix//6, ix%6]
            ax.set_yscale('log'); ax.set_xscale('log')
            catC = cat[cat['Constraint']==False]
            ux = catC[wx[ix]][ ( (catC[wx[ix]])>0) & ((catC[wy[ix]])>0) ]
            uy = catC[wy[ix]][ ( (catC[wx[ix]])>0) & ((catC[wy[ix]])>0) ]
            if len(ux)>0: ccor0 = pearsonr(np.log(ux), np.log(uy))[0]
            print(icat, strcat[icat], ix, wx[ix], wy[ix], 'n=', len(ux), 'ccor0 = %.3f' % ccor0)
            catC = cat[cat['Constraint']==True]
            ux = catC[wx[ix]][ ( (catC[wx[ix]])>0) & ((catC[wy[ix]])>0) ]
            uy = catC[wy[ix]][ ( (catC[wx[ix]])>0) & ((catC[wy[ix]])>0) ]
            if len(ux)>0: ccor1 = pearsonr(np.log(ux), np.log(uy))[0]
            print(icat, strcat[icat], ix, wx[ix], wy[ix], 'n=', len(ux), 'ccor1 = %.3f' % ccor1)
            hue = 'Time.break..s.'  #  'Constraint'
            if strcat[icat] in ['AD', 'ADCons', 'ADcut', 'ADtime']:
                sns.scatterplot(data=cat, x=wx[ix], y=wy[ix], hue=hue, alpha=0.5, ax=ax, label='%.3f' % ccor1)                
            elif strcat[icat] in ['rupt', 'charge']:
                sns.scatterplot(data=cat, x=wx[ix], y=wy[ix], hue=hue, alpha=0.5, ax=ax, label='%.3f' % ccor0)
            else: sns.scatterplot(data=cat, x=wx[ix], y=wy[ix], hue=hue, alpha=0.5, ax=ax)
            ax.plot(TU[x], TU[x], '-', c='k' );     plt.tight_layout()
            if SaveGraph: plt.savefig(OutFold+figname+OutFormat)   
    
        # figname = 'correlEbis'+strcat[icat]; fig = plt.figure(figname, figsize=(3,3), dpi=100); ax = plt.gca()
        # ax.set_yscale('log'); ax.set_xscale('log')
        # sns.scatterplot(data=cat, x='E1_EST [pN/nm]', y='E1N_EST [pN/nm]', ax=ax, label='E1N')
        # sns.scatterplot(data=cat, x='E1_EST [pN/nm]', y='E1_EST_OnForce [pN/nm]', ax=ax, label='E1onForce')
        # ax.plot(TU[x], TU[x] , '-', c='k' ); ax.legend();     plt.tight_layout()
        # if SaveGraph: plt.savefig(OutFold+figname+OutFormat)   

if RuptChargeAD:
    listfig = ['E2RuptvsChargeTU', 'E2RuptvsChargeAD', 'E1RuptvsChargeTU', 'E1RuptvsChargeAD', 'E2_EST_OnApproachRuptvsChargeTU', 'E2_EST_OnRetractionRuptvsChargeTU']
    listdata = [TU, AD, TU, AD, TU, TU]
    listy = ["E2 [pN/nm]", "E2 [pN/nm]", "E1 [pN/nm]", "E1 [pN/nm]",'E2_EST_OnApproach [pN/nm]', 'E2_EST_OnRetraction [pN/nm]']
    
    for ifig, figname in enumerate(listfig):
        fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)   
        dFtest = pd.DataFrame(data=listdata[ifig], columns=[listy[ifig],'discontType_str', 'Constraint']); ax.set_yscale('log')
        sns.boxplot(x="discontType_str", y="value", data=pd.melt(dFtest, id_vars=['discontType_str', 'Constraint'], value_vars=[listy[ifig]]), hue='Constraint', ax=ax)
        sns.swarmplot(x="discontType_str", y="value", data=pd.melt(dFtest, id_vars=['discontType_str', 'Constraint'], value_vars=[listy[ifig]]), color='0.2', alpha=0.5, hue='Constraint', ax=ax, dodge=True)
        ax.set_ylim(1e-5,1)
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

if FinalGraphs:
    
    # Warning : ajout pour eviter pb omar tout seul
    if ListIntern==['Omar']:
        ordrepres2=['igg2a', 'ConA']
        ordrepres3=['igg2a', 'cd45', 'cd3', 'ConA']
    else:
        ordrepres2 = [ 'cd11a open-Lat', 'cd11a open+Lat', 'cd11a closed-Lat', 'cd11a closed+Lat', 'cd3-Lat', 'cd3+Lat', 'cd45-Lat', 'cd45+Lat']
        ordrepres3 = ['cd11a open', 'cd11a closed', 'cd3', 'cd45']

#    dFL = RUPTAD[ (RUPTAD['InternName']=='Lucie') & (RUPTAD['myDate'] in ['2021.07.06', '2021.06.28', '2021.06.25']) ]
    for intern in ListIntern:
        dFL = RUPTAD[ (RUPTAD['InternName']==intern)  ]
        figname = intern+'DateE1'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.boxplot(x='myDate', y="E1 [pN/nm]", data=dFL, hue='Latrunculine', ax=ax)
        sns.stripplot(data= dFL, x='myDate', y="E1 [pN/nm]", hue='Latrunculine', color='0.2', alpha=0.5, ax=ax, dodge=True)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
        figname = intern+'DateE2'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.boxplot(x='myDate', y="E2 [pN/nm]", data=dFL, hue='Latrunculine', ax=ax)
        sns.stripplot(data= dFL, x='myDate', y="E2 [pN/nm]", hue='Latrunculine', color='0.2', alpha=0.5, ax=ax, dodge=True)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    wdFE1 = [ADcut, RUPTFree, RUPTAD, CHARFree, RUPTCHAR, RUPTADCHAR]
    wdFE1name = ['ADcut', 'Rupt','RuptAD' ,'Char', 'RuptChar', 'RuptAdChar']
    l = "Latrunculine"; cf = "ConditionFull"; my="myCondition"
    # wname = ['E1byE2 ', 'E1+E2 ', 'E1byE2 Ab ', 'E1 Ab ', 'E2 Ab ', 'E2_EST_OnApproach Ab', 'E1+E2 Ab ', 'E1byE2 Ab ']
    # wy = ["E1byE2", "E1+E2 [pN/nm]", "E1byE2", "E1 [pN/nm]", "E2 [pN/nm]", "E2_EST_OnApproach [pN/nm]", "E1+E2 [pN/nm]", "E1byE2"]
    # l = "Latrunculine"; cf = "ConditionFull"
    # wx = [l, l, cf, cf, cf, cf, cf, cf]
    # wsetlim = [0, 0, 0, 1, 1, 1, 0]
    wname = ['E1 Ab ', 'E2 Ab ', 'E1N Ab', 'eta Ab', 'etaN Ab']
    wy = ["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
    wx = [l, l, l, l, l]
    wsetlim = [1, 1, 0, 0, 0]
    
    wname2 = ['E1 Ab Intern', 'E2 Ab Intern', 'E1N Ab Intern', 'eta Ab Intern', 'etaN Ab Intern']
    wy2 = ["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
    wx2 = [cf, cf, cf, cf, cf]
    wsetlim2 = [1, 1, 0, 0, 0]
    
    wname3 = ['E1 Ab All', 'E2 Ab All', 'E1N Ab All', 'eta Ab All', 'etaN Ab All']
    wy3 = ["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
    wx3 = [my, my, my, my, my]
    wsetlim3 = [1, 1, 0, 0, 0]
    wsetlim3b = [0, 0, 0, 0, 0]

    wname4 = ['Jump E1 Ab All', 'Jump eta Ab All']
    wy4 = ["JumpE1 [pN/nm]", 'Jumpeta [pN/nm.sec]']
    wx4 = [cf, cf]
    wsetlim4 = [0, 0]
    E1minimal = 5e-4; E1maximal = 1
    
    orderfull =  ((cd11a0rConA+' open-Lat',cd11a0rConA+' open+Lat'), (cd11a0rConA+' closed-Lat',cd11a0rConA+' closed+Lat'),('cd3-Lat','cd3+Lat'), ('cd45-Lat','cd45+Lat'))
    orderfull2 =  ( cd11a0rConA+' open-Lat', cd11a0rConA+' closed-Lat','cd3-Lat', 'cd45-Lat' )
    orderfullfinal =  ( 'cd3-Lat', cd11a0rConA+' open-Lat', cd11a0rConA+' closed-Lat', 'cd45-Lat' )
    orderfullfinalarticle =  (  cd11a0rConA+' open-Lat', cd11a0rConA+' closed-Lat', 'cd3-Lat', 'cd45-Lat' )
    
    figname = 'RUPTADnolim E1'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.set_yscale('log') 
    sns.boxplot(data= RUPTADnolim, x="myCondition", y="E1 [pN/nm]", hue='Latrunculine', order=ordrepres3, ax=ax)
    sns.swarmplot(data= RUPTADnolim, x="myCondition", y="E1 [pN/nm]", hue='Latrunculine', order=ordrepres3, dodge=True, color='0.2', alpha=0.5, ax=ax)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    figname = 'RUPTADCHARnolim E1'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.set_yscale('log');
    distrib = RUPTADCHARnolim#[RUPTADCHARnolim['Latrunculine']==False]
 #   sns.boxplot(data= distrib , x="myCondition", y="E1 [pN/nm]", hue=, order=ordrepres3, ax=ax)
    sns.swarmplot(data= distrib, x="ConditionFull", y="E1 [pN/nm]", hue='Morpho', order=ordrepres2, dodge=True, alpha=0.5, ax=ax)
    plt.axhline(y=E1_max, color='r', linestyle='--')
    plt.axhline(y=E1_min, color='b', linestyle='--')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)   

    figname = 'Ratio Swarm'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    ax.set_yscale('log')
    sns.swarmplot(data= RUPTADCHAR, x="Morpho", y="ratio_E1E2SUM_E2", s=4, hue='Latrunculine', dodge=True, alpha=0.5, ax=ax)
    plt.axhline(y=1, color='r', linestyle='--')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)   

    figname = 'RatioSlope Swarm'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()
    sns.swarmplot(data= RUPTADCHAR, x="Morpho", y="ratio_slope", s=4, hue='Latrunculine', dodge=True, alpha=0.5, ax=ax)
    plt.axhline(y=1, color='r', linestyle='--')
    ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)   

    figname = 'Ratio_Ratio'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
    sns.scatterplot(x="ratio_E1E2SUM_E2", y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Morpho', ax=ax, alpha=0.3)
    sns.scatterplot(x="ratio_E1E2SUM_E2", y='ratio_E1E2SUM_E2', data=RUPTADCHAR, ax=ax, alpha=0.3, color='k')
    ax.errorbar(RUPTADCHARLat["ratio_E1E2SUM_E2"], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr'], color='0.2', fmt='.', alpha=0.25)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'RatioSlope_E1'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
    ax.set_xscale('log');# ax.set_yscale('log')
     #   ax.set_xlim(1e-3,1);
 #   ax.set_ylim(-5,10)
    sns.scatterplot(x="E1 [pN/nm]", y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Morpho', ax=ax, alpha=0.3)
    ax.axhline(y=1)
    ax.errorbar(RUPTADCHARLat["E1 [pN/nm]"], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr'], color='0.2', fmt='.', alpha=0.25)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'RatioSlope_Denominateur'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
    ax.set_xscale('log'); ax.set_yscale('log')
     #   ax.set_xlim(1e-3,1);
 #   ax.set_ylim(-5,10)
    sns.scatterplot(x='ratio_E1E2SUMretr_E2retr recalc', y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Morpho', ax=ax, alpha=0.3)
 #   ax.axhline(y=1)
    ax.errorbar(RUPTADCHARLat['ratio_E1E2SUMretr_E2retr recalc'], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr'], color='0.2', fmt='.', alpha=0.25)
    ax.plot(RUPTADCHARLat['ratio_E1E2SUMretr_E2retr recalc'], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr recalc'], color='0.2', alpha=0.25)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    figname = 'RatioSlope_vsSlopeRetractNeg'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
  #  ax.set_xscale('log');# ax.set_yscale('log')
     #   ax.set_xlim(1e-3,1);
 #   ax.set_ylim(-5,10)
    sns.scatterplot(x="slope_retract_neg [pN/sec]", y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Morpho', ax=ax, alpha=0.3)
    ax.axhline(y=1)
    ax.errorbar(RUPTADCHARLat["slope_retract_neg [pN/sec]"], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr'], color='0.2', fmt='.', alpha=0.25)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'RatioSlope_E2'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
    ax.set_xscale('log');# ax.set_yscale('log')
     #   ax.set_xlim(1e-3,1);
 #   ax.set_ylim(-5,10)
    sns.scatterplot(x="E2 [pN/nm]", y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Morpho', ax=ax, alpha=0.3)
    ax.axhline(y=1)
    ax.errorbar(RUPTADCHARLat["E2 [pN/nm]"], RUPTADCHARLat['ratio_E1E2SUMretr_E2retr'], color='0.2', fmt='.', alpha=0.25)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

    figname = 'RatioSlope_E1Violinplot'
    fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
  #  ax.set_xscale('log');# ax.set_yscale('log')
     #   ax.set_xlim(1e-3,1); ax.set_ylim(1e-3,1)
    sns.violinplot(x="Morpho", y='ratio_E1E2SUMretr_E2retr', data=RUPTADCHAR, hue='Latrunculine', ax=ax, inner="quartile")
    ax.set_ylim(-5,10)
    ax.axhline(y=1)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    figname = 'ParamFinal'
    fignameLat = 'ParamFinal_Lat'
    fig , axall = plt.subplots(1, len(wy), num=figname, figsize=(len(wy)*3.5, 6.), gridspec_kw={'wspace': 0.25}, dpi=100)
    figLat , axallLat = plt.subplots(1, len(wy), num=fignameLat, figsize=(len(wy)*3.5, 6.), gridspec_kw={'wspace': 0.25}, dpi=100)
    #fig.suptitle(figname, fontsize=14)
    #figLat.suptitle(fignameLat, fontsize=14)
 #   wnamefinal = ['E1 Ab ', 'E2 Ab ', 'E1N Ab', 'eta Ab', 'etaN Ab']
 #   wy = ["E1 [pN/nm]", "E2 [pN/nm]", "E1N [pN/nm]", 'eta [pN/nm.sec]', 'etaN [pN/nm.sec]']
    wdFy = [RUPTAD, RUPTADCHAR, RUPTCHAR, RUPTCHAR, RUPTCHAR]
    wymax = [1, 1, 1e-2, 10, 1e-1]
    wymin = [1e-3, 1e-3, 1e-6, 1e-3, 1e-3]
    wdy = [0.1, 0.05, 0.001, 1, 0.005]
 #   wdFEyname = ['ADcut', 'Rupt','RuptAD' ,'Char', 'RuptChar', 'RuptAdChar']
    for (iy, y) in enumerate(wy):
        SimpleEstimate(axall[iy], wy[iy], iy, wymax, wymin, wdy, wdFy[iy], 'ConditionFull', orderfullfinalarticle,
                           palette_all, 'Median', 0.95,  figname, logscale=True)
        SimpleEstimate(axallLat[iy], wy[iy], iy, wymax, wymin, wdy, wdFy[iy], 'ConditionFull', orderfull,
                           palette_all, 'Median', 0.95, fignameLat, logscale=True)
    plt.tight_layout(figname)
    plt.tight_layout(fignameLat)
    if SaveGraph: plt.figure(figname); plt.savefig(OutFold+figname+OutFormat)
    if SaveGraph: plt.figure(fignameLat); plt.savefig(OutFold+fignameLat+OutFormat)


    
    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPT, 'ConditionFull', orderfull,
                           'Median', 0.95, 'Rupt All Param' + 'Estimate Lat', logscale=True)
     
    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTAD, 'ConditionFull', orderfull,
                           'Median', 0.95, 'RuptAD E1' + 'Estimate Lat', logscale=True)

    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTADCHAR, 'ConditionFull', orderfull,
                           'Median', 0.95, 'RuptADCHAR E2 eta_' + 'Estimate Lat', logscale=True)
    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTCHAR, 'ConditionFull', orderfull,
                           'Median', 0.95, 'RuptCHAR E1N etaN eta_' + 'Estimate Lat', logscale=True)

    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPT, 'ConditionFull', orderfull2,
                           'Median', 0.95, 'Rupt All Param' + 'Estimate', logscale=True)
     
    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTAD, 'ConditionFull', orderfull2,
                           'Median', 0.95, 'RuptAD E1' + 'Estimate', logscale=True)

    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTADCHAR, 'ConditionFull', orderfull2,
                           'Median', 0.95, 'RuptADCHAR E2 eta_' + 'Estimate', logscale=True)
    MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.05,0.,0.001], RUPTCHAR, 'ConditionFull', orderfull2,
                           'Median', 0.95, 'RuptCHAR E1N etaN eta_' + 'Estimate', logscale=True) 
    
    # zone de test pour concatener les differentes choses
    if COMPAGGREG:
    
        RUPTFree['test']=RUPTFree['ConditionFull']+'-RUPT'
        RUPTAD['test']=RUPTAD['ConditionFull']+'-RUPTAD'
        RUPTCHAR['test']=RUPTCHAR['ConditionFull']+'-RUPTCHAR'
        RUPTADCHAR['test']=RUPTADCHAR['ConditionFull']+'-RUPTADCHAR'
        
        E1test=pd.concat([RUPTFree, RUPTAD], ignore_index=True)
        E1cond=E1test['test'].unique()
        E1order=(('cd11a open-Lat-RUPT','cd11a open-Lat-RUPTAD'),('cd11a closed-Lat-RUPT','cd11a closed-Lat-RUPTAD'),('cd3-Lat-RUPT', 'cd3-Lat-RUPTAD'), ('cd45-Lat-RUPT','cd45-Lat-RUPTAD'))
        dfdabest=E1test.pivot(columns='test', values='E1 [pN/nm]')
        multi_2group = dabest.load(dfdabest, idx=E1order)
        multi_2group.median_diff.plot(swarm_label='E1 (pN/nm)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'E1_aggregation'+OutFormat)
        
        figname = 'E1_aggregation_Intern'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= E1test, x='test', y='E1 [pN/nm]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
        E2test=pd.concat([RUPTFree, RUPTADCHAR], ignore_index=True)
        E2cond=E2test['test'].unique()
        E2order=(('cd11a open-Lat-RUPT','cd11a open-Lat-RUPTADCHAR'),('cd11a closed-Lat-RUPT','cd11a closed-Lat-RUPTADCHAR'),('cd3-Lat-RUPT', 'cd3-Lat-RUPTADCHAR'), ('cd45-Lat-RUPT','cd45-Lat-RUPTADCHAR'))
        dfdabest=E2test.pivot(columns='test', values='E2 [pN/nm]') 
        multi_2group = dabest.load(dfdabest, idx=E2order)
        multi_2group.median_diff.plot(swarm_label='E2 (pN/nm)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'E2_aggregation'+OutFormat)

        E2testVar=pd.concat([RUPTFree, RUPTAD, RUPTADCHAR], ignore_index=True)
        E2cond=E2testVar['test'].unique()
        E2orderVar=(('cd11a open-Lat-RUPT','cd11a open-Lat-RUPTAD','cd11a open-Lat-RUPTADCHAR'),('cd11a closed-Lat-RUPT','cd11a closed-Lat-RUPTAD','cd11a closed-Lat-RUPTADCHAR'),('cd3-Lat-RUPT', 'cd3-Lat-RUPTAD', 'cd3-Lat-RUPTADCHAR'), ('cd45-Lat-RUPT','cd45-Lat-RUPTAD','cd45-Lat-RUPTADCHAR'))
        dfdabest=E2testVar.pivot(columns='test', values='E2 [pN/nm]') 
        multi_3group = dabest.load(dfdabest, idx=E2orderVar)
        multi_3group.median_diff.plot(swarm_label='E2 (pN/nm)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'E2_aggregationVar'+OutFormat)
        
        figname = 'E2_aggregation_Intern'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= E2test, x='test', y='E2 [pN/nm]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

        figname = 'E2_aggregation_InternVar'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= E2testVar, x='test', y='E2 [pN/nm]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
        
        OTHERtest=pd.concat([RUPTFree, RUPTCHAR], ignore_index=True)
        OTHERcond=OTHERtest['test'].unique()
        OTHERorder=(('cd11a open-Lat-RUPT','cd11a open-Lat-RUPTCHAR'),('cd11a closed-Lat-RUPT','cd11a closed-Lat-RUPTCHAR'),('cd3-Lat-RUPT', 'cd3-Lat-RUPTCHAR'), ('cd45-Lat-RUPT','cd45-Lat-RUPTCHAR'))
        dfdabest=OTHERtest.pivot(columns='test', values='E1N [pN/nm]') 
        multi_2group = dabest.load(dfdabest, idx=OTHERorder)
        multi_2group.median_diff.plot(swarm_label='E1N (pN/nm)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'E1N_aggregation'+OutFormat)
        
        figname = 'E1N_aggregation_Intern'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= OTHERtest, x='test', y='E1N [pN/nm]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

        dfdabest=OTHERtest.pivot(columns='test', values='eta [pN/nm.sec]') 
        multi_2group = dabest.load(dfdabest, idx=OTHERorder)
        multi_2group.median_diff.plot(swarm_label='eta (pN/nm.sec)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'eta_aggregation'+OutFormat)
        
        figname = 'eta_aggregation_Intern'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= OTHERtest, x='test', y='eta [pN/nm.sec]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

        dfdabest=OTHERtest.pivot(columns='test', values='etaN [pN/nm.sec]') 
        multi_2group = dabest.load(dfdabest, idx=OTHERorder)
        multi_2group.median_diff.plot(swarm_label='etaN (pN/nm.sec)', custom_palette="Paired")#, color_col='InternName')
        if SaveGraph: plt.savefig(OutFold+'etaN_aggregation'+OutFormat)
        
        figname = 'etaN_aggregation_Intern'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
        fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
        sns.swarmplot(data= OTHERtest, x='test', y='etaN [pN/nm.sec]',hue='InternName', ax=ax)
        ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    
    for dFE1name, dFE1 in zip(wdFE1name, wdFE1):
        for name, x, y, setlim in zip(wname, wx, wy, wsetlim):
            figname = dFE1name + name; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca()#; f, dodge=True)
            fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
            sns.boxplot(data= dFE1, x=x, y=y, ax=ax)
            sns.swarmplot(data= dFE1, x=x, y=y, color='0.2', alpha=0.5, ax=ax)
            ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
            if setlim == 1: ax.set_ylim(E1minimal,E1maximal) 
            if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
        
        for name, x, y, setlim in zip(wname2, wx2, wy2, wsetlim2):    
            figname = dFE1name + name; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
            sns.boxplot(data= dFE1, x=x, y=y, ax=ax, order=ordrepres2, hue='InternName')
            sns.swarmplot(data= dFE1, x=x, y=y, color='0.2', alpha=0.5, ax=ax, order=ordrepres2, hue='InternName', dodge=True)
            ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
            if setlim == 1: ax.set_ylim(E1minimal,E1maximal)
            if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
            
      #  if dFE1name == 'RuptAD':
        if estimation2:
            MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.1,0.,0.001], dFE1,'ConditionFull', (cd11a0rConA+' open-Lat', cd11a0rConA+' closed-Lat', 'cd3-Lat', 'cd45-Lat'),
                          'Median', 0.95, dFE1name + 'Estimate', logscale=True)
            # MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.1,0.,0.001], dFE1,'ConditionFull', ('igg2a-Lat', cd11a0rConA+' open-Lat', cd11a0rConA+' closed-Lat', 'cd3-Lat', 'cd45-Lat'),
            #               'Median', 0.95, dFE1name + 'Estimate IgG2', logscale=True)
            MultiEstimate(wy3, [1,1,np.inf], [1e-3, 1e-3,-np.inf ], [0.1,0.,0.001], dFE1,'ConditionFull', orderfull,
                          'Median', 0.95, dFE1name + 'Estimate Lat', logscale=True)
        if dFE1name == 'RuptAD':
            for name, x, y, setlim in zip(wname3, wx3, wy3, wsetlim3b):    
                 figname = dFE1name + "Delta"+ name; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
                 sns.boxplot(data= dFE1, x=x, y="Delta"+y, ax=ax, order=ordrepres3)
                 sns.swarmplot(data= dFE1, x=x, y="Delta"+y, color='0.2', alpha=0.5, ax=ax, order=ordrepres3)
                 ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
            #     if setlim == 1: ax.set_ylim(1e-3,1)
                 if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
                 
            for name, x, y, setlim in zip(wname4, wx4, wy4, wsetlim4):    
                 figname = dFE1name + name; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
                 sns.boxplot(data= dFE1, x=x, y=y, ax=ax, order=ordrepres2)
                 sns.swarmplot(data= dFE1, x=x, y=y, color='0.2', alpha=0.5, ax=ax, order=ordrepres2)
                 ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
            #     if setlim == 1: ax.set_ylim(1e-3,1)
                 if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
            
        for name, x, y, setlim in zip(wname3, wx3, wy3, wsetlim3):    
            figname = dFE1name + name; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
            sns.boxplot(data= dFE1, x=x, y=y, ax=ax, order=ordrepres3, hue='Latrunculine')
            sns.swarmplot(data= dFE1, x=x, y=y, color='0.2', alpha=0.5, ax=ax, order=ordrepres3, hue='Latrunculine', dodge=True)
            ax.set_xticklabels(ax.get_xticklabels(),rotation=90); plt.tight_layout()
            if setlim == 1: ax.set_ylim(E1minimal,E1maximal)
            if SaveGraph: plt.savefig(OutFold+figname+OutFormat)


        
wtwait = [(0,np.inf), (0, 1), (1.01, np.inf), (0.05, 1)]
wtwait = [(0,np.inf)]

for twait in wtwait:
    strwait = 'deltaTimeWait [sec]'
    ADt = ADCons[ (ADCons[strwait]>=twait[0]) & (ADCons[strwait]<=twait[1])]
    RUPTt = RUPTFree[ (RUPTFree[strwait]>=twait[0]) & (RUPTFree[strwait]<=twait[1])]
    CHARt = CHARFree[ (CHARFree[strwait]>=twait[0]) & (CHARFree[strwait]<=twait[1])]
    
    # ADCons = ADt[ (ADt['Constraint']==True) & (ADt['E1 [pN/nm]']>E1_min) & (ADt['E1_EST [pN/nm]']>E1_EST_min)]    
    # RUPTFree = RUPTt[RUPTt['Constraint']==False]   
    # CHARFree = CHARt[CHARt['Constraint']==False]
    E2 = ADt.copy()
    E2 = E2.append(RUPTt, ignore_index=True)
    E2 = E2.append(CHARt, ignore_index=True)
    E1 = ADt.copy()
    E1 = E1.append(RUPTt, ignore_index=True)
    E1N = RUPTt.copy()
    E1N = E1N.append(CHARt, ignore_index=True)
    ETA = E2.copy()
    ETAN = E1N.copy()
    figname = 'E2 twaitrange='+str(twait); fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
    sns.boxplot(data= E2, x="ConditionFull", y="E2 [pN/nm]", ax=ax, order=ordrepres2)
    sns.swarmplot(data= E2, x="ConditionFull", y="E2 [pN/nm]", color='0.2', alpha=0.5, ax=ax, order=ordrepres2)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    figname = 'E1 ; E1min='+str(E1_min)+' twaitrange='+str(twait); fig = plt.figure(figname, figsize=(6,6), dpi=100)
    ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
    sns.boxplot(data= E1, x="ConditionFull", y="E1 [pN/nm]", ax=ax, order=ordrepres2)
    sns.swarmplot(data= E1, x="ConditionFull", y="E1 [pN/nm]", color='0.2', alpha=0.5, ax=ax, order=ordrepres2)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    figname = 'E2 Lat twaitrange='+str(twait); fig = plt.figure(figname, figsize=(6,6), dpi=100)
    ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
    sns.boxplot(data= RUPTFree, x='Latrunculine', y="E2 [pN/nm]", ax=ax)
    sns.swarmplot(data= E1, x="Latrunculine", y="E2 [pN/nm]", color='0.2', alpha=0.5, ax=ax)
    if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    # figname = 'E1+E2 ; E1min='+str(E1_min)+' twaitrange='+str(twait); fig = plt.figure(figname, figsize=(6,6), dpi=100)
    # ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
    # sns.boxplot(data= RUPTFree, x='Latrunculine', y="E1+E2 [pN/nm]", ax=ax)
    # sns.swarmplot(data= E1, x="Latrunculine", y="E1+E2 [pN/nm]", color='0.2', alpha=0.5, ax=ax)
    # if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
    # figname = 'E1E2SUM_EST; E1min='+str(E1_min)+' twaitrange='+str(twait); fig = plt.figure(figname, figsize=(6,6), dpi=100)
    # ax = plt.gca(); fig.suptitle(figname, fontsize=14); ax.set_yscale('log') 
    # sns.boxplot(data= RUPTFree, x='Latrunculine', y='E1E2SUM_EST [pN/nm]', ax=ax)
    # sns.swarmplot(data= E1, x="Latrunculine", y='E1E2SUM_EST [pN/nm]', color='0.2', alpha=0.5, ax=ax)
    # if SaveGraph: plt.savefig(OutFold+figname+OutFormat)    
    
    
    if twait == (0,np.inf):
        print(len(RUPTFree))
        figname = 'scatterRuptE1E2EST'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim(1e-3,1); ax.set_ylim(1e-3,1)
        sns.scatterplot(x='E1_EST_OnDeltaForce [pN/nm]', y='E2_EST_OnApproach [pN/nm]', data=RUPTFree, hue='Latrunculine', ax=ax, alpha=0.3)
        ax.plot(RUPTFree['E1_EST_OnDeltaForce [pN/nm]'], RUPTFree['E1_EST_OnDeltaForce [pN/nm]'], color='0.2')
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)
        figname = 'scatterRuptE1E2'; fig = plt.figure(figname, figsize=(6,6), dpi=100); ax = plt.gca(); fig.suptitle(figname, fontsize=14)  
        ax.set_xscale('log'); ax.set_yscale('log')
        ax.set_xlim(1e-3,1); ax.set_ylim(1e-3,1)
        sns.scatterplot(x='E1 [pN/nm]', y='E2 [pN/nm]', data=RUPTFree, hue='Latrunculine', ax=ax, alpha=0.3)
        ax.plot(RUPTFree['E1 [pN/nm]'], RUPTFree['E1 [pN/nm]'], color='0.2')
        if SaveGraph: plt.savefig(OutFold+figname+OutFormat)

if fractionsPH:
    #======================================================================
    # new code hereafter
    #======================================================================
    
    TUset=TU[TU['E1 [pN/nm]']>10^(-5)]# to remove saturating values
    ADset=AD[AD['E1 [pN/nm]']>10^(-5)]
    
    ADset['discontType_str']='adhesion' #fusion two types for adhesion
    
    TUNL = TUset[TUset['Latrunculine']==False]
    ADNL = ADset[ADset['Latrunculine']==False]
    TUL = TUset[TUset['Latrunculine']==True]
    ADL = ADset[ADset['Latrunculine']==True]
    
    seuil=0.02 #setting threshold on E1
    
    print('-------------------------------------')
    print('Threshold on E1 is set to ', seuil, 'pN/nm')
    print('-------------------------------------')
    print('Fractions above', seuil, ' (NO lat)')
    print('---')
    allmolecule=['cd3', 'cd45', 'cd11a']
    for molecule in allmolecule:
        TUx=TUNL[TUNL['myCondition']==molecule]
        ADx=ADNL[ADNL['myCondition']==molecule]
        stock=[TUx, ADx]
        stockid=['TU'+molecule, 'AD'+molecule]
        for idFraction, dFraction in enumerate(stock):
            value=dFraction.groupby('discontType_str')['E1 [pN/nm]'].apply(lambda c: (c>seuil).sum()/len(c))
            print(stockid[idFraction])
            print(value)
            print('---')
            
    print('-------------------------------------')
    print('Fractions above', seuil, '(ONLY lat)')
    print('---')
    allmolecule=['cd3', 'cd45', 'cd11a']
    for molecule in allmolecule:
        TUx=TUL[TUL['myCondition']==molecule]
        ADx=ADL[ADL['myCondition']==molecule]
        stock=[TUx, ADx]
        stockid=['TU'+molecule, 'AD'+molecule]
        for idFraction, dFraction in enumerate(stock):
            value=dFraction.groupby('discontType_str')['E1 [pN/nm]'].apply(lambda c: (c>seuil).sum()/len(c))
            print(stockid[idFraction])
            print(value)
            print('---')

if STACKEDSTATS:
    
    StatsDataFrame(RUPTADCHAR, 'ConditionFull', 'Morpho')
    
    

if CLOSEFIGSATEND: plt.close(fig='all')
