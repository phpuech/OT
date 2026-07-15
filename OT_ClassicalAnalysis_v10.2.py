# -*- coding: utf-8 -*-
"""
Éditeur de Spyder

Ceci est un script temporaire.

PHP -> adapte de la version 9 de lucie, elle meme adaptee de gautier like plots

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import dabest
import scipy.stats as ss
import statsmodels.api as sa
import scikit_posthocs as sp #  conda install -c conda-forge scikit-posthocs 
import warnings
warnings.filterwarnings("ignore")
import time, os

# user modifies this -------------------------------------------------------------------------
version=10.2
# where is my data
inputpath="/home/php/Bureau/Fabio-data/211108-Fabio/ClassicalAnalysis/"
files = [f for f in os.listdir(inputpath) if f.endswith('.xslx')]
files.sort()

SPECIF=True #only WHO considered pour les graphes de fraction if TRUE ; if FALSE : all interns
WHO='Gautier' #'Gautier', 'Lucie'

# pooling all dates if true
pooled=True 
# for plotting dabest , do we separate the different times or not
DabestSeparatedvsTime=False # False = all in one
# defining adhesion or not
CUTOFFFORCE=True
cutoffforce=5#pN #to get all values, set it to 0pN

showplots=False

# --------------------------------
# what is my data
conditionsref=['igg2a-Lat', 'cd11a open-Lat', 'cd11a closed-Lat', 'cd3-Lat', 'cd45-Lat']
rejectedconditions = ['bCD3_aCD11a_avant', 'SCD3_aCD11a', 'cd4', 'aCD4', 'nues', 'SCD45_aCD11a', 'bCD45_aCD11a_apres', 'bCD3_aCD11a_apres', 'bCD45_aCD11a_avant']

control="igg2a"+'-Lat'
reference='cd11a open'+'-Lat'

conditionsnoref=[ 'cd11a open-Lat', 'cd11a closed-Lat', 'cd3-Lat', 'cd45-Lat']
conditionsnorefLat=[ 'cd11a open-Lat', 'cd11a open+Lat','cd11a closed-Lat','cd11a closed+Lat', 'cd3-Lat','cd3+Lat', 'cd45-Lat', 'cd45+Lat']
orderLAT =  (('cd11a open-Lat','cd11a open+Lat'), ('cd11a closed-Lat','cd11a closed+Lat'),('cd3-Lat','cd3+Lat'), ('cd45-Lat','cd45+Lat'))
orderAB =  ( 'cd11a open-Lat','cd11a closed-Lat','cd3-Lat', 'cd45-Lat' )
orderABctrl=( 'igg2a-lat', 'cd11a open-Lat','cd11a closed-Lat','cd3-Lat', 'cd45-Lat' )
# --------------------------------

#for plotting
maxslope=0.5#pN/nm (only for plotting symetric, removes few points)
cutoffdist=1000#nm
cutoffslope=0.025#pN/nm
angle=90 # for rotating x axis label - set it to zero for parallel

# for php sum files
columns =["X","Filename","Type.of.event","Type.of.event..corrected.","Date",\
          "Hour","Condition","Couple","Bead","Cell","Direction","Sens",\
          "Constant..pN.nm.","Force.contact..pN.","Time.break..s.","Distance..µm.",\
          "Speed..µm.s.","Trace.s.fit.convergence","Slope..pN.nm.","X.error.",\
          "SD.baseline.retrace..pN.","Min.Force..pN.","Position.Min.Force..nm.",\
          "Jump.force..pN.","Jump.end..nm.","Retrace.s.fit.convergence",\
          "Retrace.Fitting.frame...pts.","Pente..pN.nm.","Tube.cassé..",\
          "Aire.au.min..pN.nm.","Aire.au.jump..pN.nm.", \
          "Latrunculine"]

# ------------------------------------------------------------------------
# function definitions

def CreatePopulateLat(df):
    df['Latrunculine']=False
    liste=['lat', 'Lat']
    for i in liste:
        df['Latrunculine'] = np.where(df['Condition'].str.contains(i), True, df['Latrunculine'])
    return df

def RenameCondition(df):
    print('Before', df['Condition'].unique())
    df['Condition'] = df['Condition'].str.lower()
    df['Condition'] = df['Condition'].str.replace(r'acd', 'cd')
    liste=['témoin', 'ucht1', 'lata', '2µM']
    for i in liste:
        df['Condition'] = df['Condition'].str.replace('_'+i, '')
        df['Condition'] = df['Condition'].str.replace(i+'_', '')
    print('After', df['Condition'].unique())
    return df

def AttributeIntegrinType(df):
    df.loc[(df['Condition']=='cd11a')&(df['Intern']=='Gautier') , 'Condition']='cd11a closed'
    df.loc[(df['Condition']=='cd11a')&(df['Intern']=='Lucie') , 'Condition']='cd11a open'
    return df

def CreateFullCondition(df):
    df['ConditionFull'] = np.where(df['Latrunculine']==True , df['Condition']+'+Lat', df['Condition']+'-Lat')
    return df

def Convert(list):
    return tuple(list)

def SetReference(liste, item):
    tmp=list(liste)
    tmp.remove(item)
    tmp.insert(0, item)
    listemod=tuple(tmp)
    return listemod
  
def DabestPlot(title, df, values, conditions, reference, pooled):
    subdf=df[(df['ConditionFull'].isin(conditions))]
    name=title+'-pooled'
    by=SetReference(conditions, reference)
    dfdabest=subdf.pivot(columns='ConditionFull', values=values) 
    shared_control2=dabest.load(dfdabest,idx=by)
    shared_control2.median_diff.plot(swarm_label=name,
                                         contrast_label="Median difference")
    plt.tight_layout()

def CoupledDabestPlot(title, df, values, conditions, by, pooled):
    subdf=df[(df['ConditionFull'].isin(conditions))]
    name=title+'-pooled'
    dfdabest=subdf.pivot(columns='ConditionFull', values=values) 
    shared_control2=dabest.load(dfdabest,idx=by)
    shared_control2.median_diff.plot(swarm_label=name,
                                         contrast_label="Median difference")
    plt.tight_layout()

def ChartBoxPlot(df, values,  hue, conditions, angle):
    chart=sns.catplot(x="ConditionFull", y=values,
             hue=hue,  data=df,
             kind='box', height=5, aspect=.7,  order=conditions)
    chart.set_xticklabels(rotation=angle, horizontalalignment='right')
    plt.tight_layout()
    return chart

def BoxDotsPlot(df1, values, order):
    sns.boxplot(x="ConditionFull", y=values,
                data=df1, order=order)
    sns.stripplot(x="ConditionFull",  y=values,
                data=df1, alpha=0.5,  dodge=True,  order=order)
    plt.tight_layout()
    
def BoxDotsPlotHue(df1, values, order, hue):
    sns.boxplot(x="ConditionFull", y=values,
                data=df1, order=order, hue=hue)
    sns.stripplot(x="Condition",  y=values,
                data=df1, alpha=0.5,  dodge=True,  order=order, hue=hue)
    plt.tight_layout()

def BoxDotsPlotTwoColors(df1, df2, values, order):
    sns.boxplot(x="ConditionFull", y=values,
                data=df1, order=order)
    sns.stripplot(x="ConditionFull",  y=values,
                data=df1, alpha=0.5,  dodge=True, order=order)
    sns.stripplot(x="ConditionFull",  y=values, #for two colors for infinite tubes
                data=df2, color='black', alpha=0.5,  dodge=True, order=order)
    plt.tight_layout()

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
            noadh.append(counts.loc[(j,df1,t)][('Filename','count')])
        else:
            noadh.append(0)
        if counts.index.isin([(j,df2,t)]).any():  
            adhesion.append(counts.loc[(j,df2,t)][('Filename','count')])
        else:
            adhesion.append(0)
        if counts.index.isin([(j,df3,t)]).any():
            tube.append(counts.loc[(j,df3,t)][('Filename','count')])
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
    
def FractionPlotTubes(r, counts, df1, df2,durees,k,  unique, names, t, angle):
# def FractionPlot(df1, df2, df3, df4, unique, names, t, angle):
    raw={}
    # noadh=[]
    # adhesion=[]
    tubes=[]
    infinites=[]
    totals=[]
    #print(r, names)    listemod=tuple(tmp)
    for j in unique:
        if counts.index.isin([(j,df1,t)]).any():
            tubes.append(counts.loc[(j,df1,t)][('Filename','count')])
        else:
            tubes.append(0)
        if counts.index.isin([(j,df2,t)]).any():
            infinites.append(counts.loc[(j,df2,t)][('Filename','count')])
        else:
            infinites.append(0)
        # if counts.index.isin([(j,df4,t)]).any():
        #     infinite.append(counts.loc[(j,df4,t)][('X','count')])
        # else:
        #     infinite.append(0)
    raw={'tub':tubes, 'tubinf':infinites}#, 'inftub':infinite}stock_nodup
    #print(raw)
    dcount=pd.DataFrame(raw);  
    # sums
    totals = [i+j for i,j in zip(dcount['tub'], dcount['tubinf'])]#, dcount['inftub'])]
    # ratios - COMMENT : np.divide(x,0) renvoit 0 !!! pout Trueeviter les cas ou il y a pas de mesures
    tube = np.nan_to_num([np.divide(i, j) * 100 for i,j in zip(dcount['tub'], totals)])
    # print(no)
    infinitetube = np.nan_to_num([np.divide(i, j)  * 100 for i,j in zip(dcount['tubinf'], totals)])
    # print(adh)
    #print(tube)
    # inftub = np.nan_to_num([np.divide(i, j)  * 100 for i,j in zip(dcount['inftub'], totals)])
    # print(inftub)
    barWidth = 0.5
    # Create green Bars : tubes
    # ax[k].bar(r, inftub, color='darkgreen', edgecolor='white', width=barWidth, label='Inf. Tubes')
    # Create orange Bars : adhesion
    ax[k].bar(r, tube, color='lightgreen', edgecolor='white', width=barWidth, label='Tubes')
    # Create blue Bars : non adhesion
    ax[k].bar(r, infinitetube, bottom=tube, color='darkgreen', edgecolor='white', width=barWidth, label='Tubes inf.')

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

# -------------------------------------------------------------------------
# to know when we perform the analysis
now = time.strftime("%c")
today =time.strftime("%Y")+time.strftime("%m")+time.strftime("%d")
heure = time.strftime("%H")+time.strftime("%M")+time.strftime("%S")
maintenant = today + "-" + heure
print(now)
print( "---------------------------------------------------")
# close the residual plots if anyplt.ylabel("Force (pN)")
plt.close(fig='all')

# -------------------------------------------------------------------------
# Merging the source files from force_classifier
# find them all
files = [f for f in os.listdir(inputpath) if f.endswith('.xlsx')]
# # prepare output as excel
outputfile='MERGED-RENAMED-STUDENTS.xlsx'
outputpath = inputpath + maintenant+ '/'
datasave = outputpath+maintenant +'-'+outputfile
# prepare saving
if not os.path.exists(outputpath):
     os.makedirs(outputpath)
# init the saving datafile    
stock=pd.DataFrame(columns=columns)
#load and merge plus adding names
for fichier in files :
    localfichier=inputpath+fichier
    work=pd.read_excel(localfichier, index_col=0, names=columns)
    if 'Gautier'  in fichier: 
        print('Gautier : '+str(len(work)))
        renamedwork=RenameCondition(work)
        renamedwork['Intern']='Gautier'
        renamedwork=AttributeIntegrinType(renamedwork)
    if 'Lucie' in fichier: 
        print('Lucie : '+str(len(work)))
        work=CreatePopulateLat(work)  
        renamedwork=RenameCondition(work)
        renamedwork['Intern']='Lucie'
        renamedwork=AttributeIntegrinType(renamedwork)
    print(fichier, 'modified')
    stock=stock.append(renamedwork, ignore_index=True)
    
# check for duplicates on filename which is unique since contains b, c, date, time
#duplicate detection 
dFload=stock.drop_duplicates(subset='Filename').copy()
# # saving new data without duplicates
dFload.to_excel(outputpath+maintenant+'-nodup-renamed-'+outputfile,  index=False) 

# print to screen basic descriptives
print('===========')
print('Original total #FC = ', len(stock)) 
print('Duplicates ?', len(stock) != len(dFload))
print('After duplicates removal, #FC = ', len(dFload))
print('===========')                

#---------------------------------------------------------------------------------------
# removing unwanted conditions
print('Keeping only wanted conditions')
print('Before', dFload['Condition'].unique())
for cond in rejectedconditions: dFload = dFload[dFload['Condition']!=cond]#; print('Removing '+cond)
print('After', dFload['Condition'].unique())
print('===========')

# Data with forces below seuilforce are addressed as NAd
if CUTOFFFORCE : dFload.loc[dFload["Jump.force..pN."]<cutoffforce, 'Type.of.event..corrected.']='NAd' ; print('WARNING, forces < '+str(cutoffforce)+ ' pN are turned to NAd');
#select only relevant conditions

# lazy continuity with previous code, creating complete Lat names
dF=CreateFullCondition(dFload)
  
# define subsets of data of interest
NAd=dF[dF['Type.of.event..corrected.']=='NAd'].copy()
AD = dF[dF['Type.of.event..corrected.']=='AD'].copy()
TU = dF[dF['Type.of.event..corrected.']=='TU'].copy()
TI = dF[(dF['Type.of.event..corrected.']=='TU') & (dF["Tube.cassé.."]=='no')].copy()
TF = dF[(dF['Type.of.event..corrected.']=='TU') & (dF["Tube.cassé.."]=='yes')].copy()

# create an easy df for adhesion and tubes only data, above threshold
Events=AD.copy()
Events=pd.concat([Events, TU])

#-------------------------------------------------------------------------
# recompose the classical scatter plots  for the entire set, no separation of condition    
fig1=plt.figure("F vs slope", dpi=100)
plt.scatter(dF["Pente..pN.nm."],dF["Jump.force..pN."], color='black', label='', alpha=0.25, s=48)
plt.scatter(TU["Pente..pN.nm."],TU["Jump.force..pN."], color='green', label='', alpha=0.75, s=24)
plt.scatter(AD["Pente..pN.nm."],AD["Jump.force..pN."], color='red', label='', alpha=0.75, s=24)
plt.axhline(y=cutoffforce, alpha=0.5, color='red', ls="--", label=str(cutoffforce)+' pN')
plt.axvline(x=cutoffslope, alpha=0.5, color='green', ls="--", label=str(cutoffslope)+' pN/nm')
plt.axvline(x=-cutoffslope, alpha=0.5, color='green', ls="--", label='')
plt.xlabel('Pente (pN/nm)')
plt.ylabel('Force (pN)')
plt.xlim(-maxslope, maxslope)
plt.tight_layout
plt.savefig(outputpath+'scatterFpentes.png')

fig2=plt.figure("F vs position", dpi=100)
cutoffforce=5#pN
cutoffdist=1000#nm
cutoffslope=0.025#pN/nm
plt.scatter(dF["Jump.end..nm."],dF["Jump.force..pN."], color='black', label='', alpha=0.25, s=48)
plt.scatter(TU["Jump.end..nm."],TU["Jump.force..pN."], color='green', label='', alpha=0.75, s=24)
plt.scatter(AD["Jump.end..nm."],AD["Jump.force..pN."], color='red', label='', alpha=0.75, s=24)
plt.axhline(y=cutoffforce, alpha=0.5, color='red', ls="--", label=str(cutoffforce)+' pN')
plt.axvline(x=cutoffdist, alpha=0.5, color='green', ls="--", label=str(cutoffdist)+' nm')
plt.ylabel('Force (pN)')
plt.xlabel('Distance (nm)')
plt.tight_layout   
plt.savefig(outputpath+'scatterFdists.png')
# -------------------------------------------------------
# print to screen basic descriptives
print("Total number of curves : ", len(dF))
print("Adhesion :", len(AD))
print("Tubes total :", len(TU))
print(" - Tubes finis :", len(TF))
print(" - Tubes infinis :", len(TI))
print("Non adhesion: ", len(NAd))
print('Rejected:', len(dF)-(len(AD)+len(TU)+len(NAd)))
print('===========')

# -------------------------------------------------------
# prepare countings for the stacked bar plots
if SPECIF:
    dFSpecif=dF[dF['Intern']==WHO]
    TUSpecif=TU[TU['Intern']==WHO]
    TISpecif=TI[TI['Intern']==WHO]
    ADSpecif=AD[AD['Intern']==WHO]
    NAdSpecif=NAd[NAd['Intern']==WHO]
    
    counts=dFSpecif.groupby(['ConditionFull', 'Type.of.event..corrected.','Time.break..s.']).agg(['count'])
    tubecounts=TUSpecif.groupby(['ConditionFull',"Tube.cassé..",'Time.break..s.']).agg(['count'])
    
    uniquecondition=dFSpecif['ConditionFull'].unique()
    timescondition=dFSpecif['Time.break..s.'].unique()
    dayscondition=dFSpecif['Date'].unique()

else:
    counts=dF.groupby(['ConditionFull', 'Type.of.event..corrected.','Time.break..s.']).agg(['count'])
    tubecounts=TU.groupby(['ConditionFull',"Tube.cassé..",'Time.break..s.']).agg(['count'])
    
    uniquecondition=dF['ConditionFull'].unique()
    timescondition=dF['Time.break..s.'].unique()
    dayscondition=dF['Date'].unique()

# experimental conditions
unique=sorted(uniquecondition, reverse=True)
categ=len(unique) 
experiment=Convert(unique) # convert to list

# experimental contact times
times=sorted(timescondition, reverse=False)
durees=len(times)
temps=Convert(times)

# recover the days (repeats)
daysunique=sorted(dayscondition, reverse=False)
ndays=len(daysunique)
daily=Convert(daysunique)

# prepare for barplots as category plots for counts and for names fbottom=[i+j for i,j in zip(tub, adh)]or other plotsprint(no)
r=Convert([i for i in range(categ)])
names=experiment

fig3a, ax = plt.subplots(1,durees, figsize=(3*durees,5))
k=0
for t in times:
    FractionPlot(r, counts, 'NAd', 'AD', 'TU', durees,k, unique, names, t, angle)
    k=k+1
plt.tight_layout()
plt.savefig(outputpath+'Fractions.png')
fig3b, ax = plt.subplots(1,durees, figsize=(3*durees,5))
l=0
for t in times:
    FractionPlotTubes(r, tubecounts, 'yes', 'no', durees,  l, unique, names, t, angle)
    l=l+1
plt.tight_layout()
plt.savefig(outputpath+'Fractions-tubes.png')

# -------------------------------------------------------
# print to screen basic descriptives
print('Summary :')
if SPECIF: print(WHO)
print(categ, " conditions, ", durees, " contact times, ", ndays, "days of exp.")
print('===========')
print('Antibodies -> ', unique)
print('Contact times -> ',times)
print("Days --> ", daysunique)
print('===========')

# ----------------------------------------------
# classical boxplots
fig4a = plt.figure("Tubes, Ffin - all pooled", dpi=100)
BoxDotsPlotTwoColors(TU, TI, "Jump.force..pN.", conditionsref)
plt.ylabel('Tube force (pN)')
plt.ylim(0,70)
plt.savefig(outputpath+'Tubes_Ffin.png')

fig5a = plt.figure("Tubes, L - all pooled", dpi=100)
BoxDotsPlotTwoColors(TU, TI, "Jump.end..nm.", conditionsref)
plt.ylabel('Tube length (nm)')
plt.ylim(0,20000)
plt.savefig(outputpath+'Tubes_L.png')

fig4b = plt.figure("Adh, F - all pooled", dpi=100)
BoxDotsPlot(AD, "Jump.force..pN.", conditionsref)
plt.ylim(0,70)
plt.ylabel('Adhesion force (pN)')
plt.savefig(outputpath+'Adh_F.png')

fig5b = plt.figure("Adh, L - all pooled", dpi=100)
BoxDotsPlot(AD, "Jump.end..nm.", conditionsref)
plt.ylabel('Adhesion length (nm)')
plt.ylim(0,2500)
plt.savefig(outputpath+'Adh_L.png')

# fig4c = plt.figure("Tubes, Fmin - all pooled", dpi=100)
# BoxDotsPlotTwoColors(TU, TI, "Min.Force..pN.", conditionsref)
# plt.ylabel('Min force (pN)')
# plt.ylim(0,70)
# plt.savefig(outputpath+'Tubes_Fmin.png')

chart1=ChartBoxPlot(Events, "Jump.force..pN.",  "Type.of.event..corrected.", conditionsref, angle)
chart1.set_xticklabels(rotation=angle, horizontalalignment='right')
plt.tight_layout
plt.savefig(outputpath+'AdhvsTubes_Ffin.png')

chart1b=ChartBoxPlot(Events, "Min.Force..pN.", "Type.of.event..corrected.", conditionsref, angle)
chart1b.set_xticklabels(rotation=angle, horizontalalignment='right')
plt.tight_layout
plt.savefig(outputpath+'AdhvsTubes_Fmin.png')

# Events['ForceDifference']=Events['Min.Force..pN.']-Events['Jump.force..pN.']
# chart1c=ChartBoxPlot(Events, "ForceDifference", "Type.of.event..corrected.", conditionsref, angle)
# chart1c.set_xticklabels(rotation=angle, horizontalalignment='right')
# plt.tight_layout
# plt.savefig(outputpath+'AdhvsTubes_Fmin-Ffin.png')

chart2=ChartBoxPlot(Events, "Jump.end..nm.",  "Type.of.event..corrected.", conditionsref, angle)
chart2.set_xticklabels(rotation=angle, horizontalalignment='right')
plt.tight_layout
plt.savefig(outputpath+'AdhvsTubes_L.png')    

# TI vs TU
chart3=ChartBoxPlot(Events, "Jump.force..pN.",  "Tube.cassé..", conditionsref, angle)
chart3.set_xticklabels(rotation=angle, horizontalalignment='right')
plt.savefig(outputpath+'TubesFfinisvsnon_Ffin.png') 

chart3=ChartBoxPlot(Events, "Min.Force..pN.",  "Tube.cassé..", conditionsref, angle)
chart3.set_xticklabels(rotation=angle, horizontalalignment='right')
plt.savefig(outputpath+'TubesFfinisvsnon_Fmin.png') 

# ------------------------------------------------------------------------------------
# Dabest estimate plots
CoupledDabestPlot("F tube (pN) - all contact times", TU, "Jump.force..pN.", conditionsnorefLat, orderLAT,pooled)
plt.savefig(outputpath+'FtubesdabestLAT-alltimes.png')
  
CoupledDabestPlot("F adh (pN) - all contact times", AD, "Jump.force..pN.", conditionsnorefLat, orderLAT, pooled)
plt.savefig(outputpath+'FadhdabestLAT-alltimes.png')

DabestPlot("F tube (pN) - all contact times", TU, "Jump.force..pN.", conditionsref, control, pooled)
plt.savefig(outputpath+'Ftubesdabest-alltimes.png')
  
DabestPlot("F adh (pN) - all contact times", AD, "Jump.force..pN.", conditionsref, control, pooled)
plt.savefig(outputpath+'Fadhdabest-alltimes.png')

#--------------------------------------------------------------------------------
# fig4d = plt.figure("Tubes, Delta=Min-Fin - all pooled", dpi=100)
# BoxDotsPlot(Events, 'ForceDifference', conditions)
# plt.ylabel('Force difference (cut below 20pN) (pN)')
# plt.ylim(0,20)
# plt.savefig(outputpath+'Tubes_Fmin-Ffin.png')

# DabestPlot("Fmin-Ffin, dabest", Events, 'ForceDifference', [conditions[0], conditions[2]], conditions[0], pooled, time=0)
# plt.savefig(outputpath+'Fmin-Ffin-dabest.png')CoupledDabestPlot("F tube (pN) - all contact times", TU, "Jump.force..pN.", conditionsnorefLat, orderLAT,pooled, 0)
#plt.savefig(outputpath+'FtubesdabestLAT-alltimes.png')

# CoupledDabestPlot("Fmin-Ffin, dabest", Events, 'ForceDifference', conditionsnorefLat, orderLAT, pooled, time=0)
# plt.savefig(outputpath+'Fmin-Ffin-dabest-coupled.png')
#--------------------------------------------------------------------------------
# tension calcuLation, only used if , time=0Latru[nculin is present
TU['Tension']=TU["Jump.force..pN."]**2/(8*np.pi**2*2e-19)*1e-21
INFTU=TU[TU["Tube.cassé.."]=='no']

fig6 = plt.figure("Tension", dpi=100)
BoxDotsPlotTwoColors(TU, INFTU, 'Tension', conditionsnorefLat)
plt.savefig(outputpath+'Tension.png')

# taking the medians of the cases w/o (barplot) and w/ lat (boxplot)
medTUbasaltension_aCD11aopen = TU[TU['ConditionFull']==conditionsnorefLat[1]].median()['Tension']
medTUbasaltension_aCD11aclosed = TU[TU['ConditionFull']==conditionsnorefLat[3]].median()['Tension']
medTUbasaltension_aCD3 = TU[TU['ConditionFull']==conditionsnorefLat[5]].median()['Tension']
medTUbasaltension_aCD45 = TU[TU['ConditionFull']==conditionsnorefLat[7]].median()['Tension']

# medTUtension_aCD11aopen = TU[TU['ConditionFull']==conditionsnorefLat[0]].median()['Tension']
# medTUtension_aCD11aclosed = TU[TU['ConditionFull']==conditionsnorefLat[2]].median()['Tension']
# medTUtension_aCD3 = TU[TU['ConditionFull']==conditionsnorefLat[4]].median()['Tension']
# medTUtension_aCD45 = TU[TU['ConditionFull']==conditionsnorefLat[6]].median()['Tension']

# print('Median tension, no Lat/Lat/diff')
# print('aCD3:', medTUtension_aCD3, medTUbasaltension_aCD3, medTUtension_aCD3- medTUbasaltension_aCD3)
# #print('aCD45:', medTUtension_aCD45, medTUbasaltension_aCD45,medTUtension_aCD45- medTUbasaltension_aCD45 )
# print('aCD11a', medTUtension_aCD11a, medTUbasaltension_aCD11a, medTUtension_aCD11a- medTUbasaltension_aCD11a)
# ("-------------------")

#yvalues=[medTUtension_aCD11aopen- medTUbasaltension_aCD11aopen, medTUtension_aCD11aclosed- medTUbasaltension_aCD11aclosed, medTUtension_aCD3- medTUbasaltension_aCD3, medTUtension_aCD45- medTUbasaltension_aCD45]
# fig7 = plt.figure("W0, bar", dpi=100)
# sns.barplot(x=[conditionsnorefLat[0], conditionsnorefLat[2], conditionsnorefLat[4], conditionsnorefLat[6]], y=yvalues)
# plt.ylabel('W0')
# plt.savefig(outputpath+'W0bar.png')
#----
#energy calcuLation
TU['W0']=TU['Tension']
TU['W0']=np.where(TU['ConditionFull']=='aCD3-Lat', TU['Tension']-medTUbasaltension_aCD3, TU['W0'])
TU['W0']=np.where(TU['ConditionFull']=='aCD45-Lat', TU['Tension']-medTUbasaltension_aCD45, TU['W0'])
TU['W0']=np.where(TU['ConditionFull']=='aCD11a open-Lat', TU['Tension']-medTUbasaltension_aCD11aopen, TU['W0'])
TU['W0']=np.where(TU['ConditionFull']=='aCD11a closed-Lat', TU['Tension']-medTUbasaltension_aCD11aclosed, TU['W0'])
#TU['W0']=np.where(TU['Condition']=='IgG2a', np.nan, TU['W0'] )
W0=TU[TU['ConditionFull'].isin(Convert(conditionsnoref))]#&(TU['Condition']!='IgG2a')]

# bloxplot
fig8 = plt.figure("W0, box", dpi=100)
BoxDotsPlot(W0, 'W0', conditionsnoref)
plt.ylabel('W0')
plt.tight_layout()
plt.savefig(outputpath+'W0box.png')

# dabest estimation
DabestPlot("W0, dabest", W0, 'W0',conditionsnoref, reference, pooled)
plt.savefig(outputpath+'W0-dabest.png')

#-------
# stat plot of p values
# print("Conover test + Bonferonni")
# pc=sp.posthoc_conover(W0, val_col='W0', group_col='Condition', p_adjust = 'bonferroni')
# print(pc)
# fig9 = plt.figure("Significance", dpi=100)
# heatmap_args = {'linewidths': 0.25, 'linecolor': '0.5', 'clip_on': False, 'square': True, 'cbar_ax_bbox': [0.80, 0.35, 0.04, 0.3]}
# sp.sign_plot(pc, **heatmap_args)
# plt.savefig(outputpath+'Tests.png')
#-----------
# ADW0dabest=W0.pivot(columns='Condition', values='W0') 
# #print(W0dabest)
# shared_control=dabest.load(W0dabest, idx=convert(A))
# shared_control.median_diff.plot(swarm_label="W0",
#                                      contrast_label="Median difference")
# plt.tight_layout()
# plt.savefig(outputpath+'W0dabest.png')
     
if showplots:
    plt.show()
else:
    plt.close(fig='all')