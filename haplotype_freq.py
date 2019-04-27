#!/usr/bin/python
from __future__ import division
import pandas as pd
import sys,os
import numpy as np
import scipy.stats as stats
from scipy.stats import chi2_contingency
import argparse
def OR_CI_abcd(a,b,c,d):
   print a,b,c,d
   chi2, pvalue_chi2, dof, ex= stats.chi2_contingency([[a,c],[b,d]],correction=False)
   chi2_corrected, pvalue_corrected, dof_corrected, ex_corrected= stats.chi2_contingency([[a,c],[b,d]])
   odd_ratio_fisher,pvalue_fisher=stats.fisher_exact([[a,c],[b,d]])
   if a==0 or b==0 or c==0 or d==0: 
     a+=0.5
     b+=0.5
     c+=0.5
     d+=0.5
   odd_ratio=(a*d)/(b*c)
   se=np.sqrt(1/a+1/b+1/c+1/d)
   ci_l=np.exp(np.log(odd_ratio)-1.96*se)
   ci_h=np.exp(np.log(odd_ratio)+1.96*se)
   return odd_ratio,se,ci_l,ci_h,pvalue_fisher,chi2,pvalue_chi2,chi2_corrected,pvalue_corrected
#print OR_CI_abcd(28,178,442,2192)
#sys.exit()

def OR_CI(row):
   a=row['Patient Count']
   b=row['Healthy Count']
   c=row['Remaining Patients']
   d=row['Remaining Healthy']
   ###print a,b,c,d
   odd_ratio_fisher,pvalue_fisher=stats.fisher_exact([[a,c],[b,d]])
   if a==0 or b==0 or c==0 or d==0: 
     a+=0.5
     b+=0.5
     c+=0.5
     d+=0.5
   odd_ratio=(a*d)/(b*c)
   chi2, pvalue_chi2, dof, ex= stats.chi2_contingency([[a,c],[b,d]],correction=False)
   chi2_corrected, pvalue_corrected, dof_corrected, ex_corrected= stats.chi2_contingency([[a,c],[b,d]])
   se=np.sqrt(1/a+1/b+1/c+1/d)
   ci_l=np.exp(np.log(odd_ratio)-1.96*se)
   ci_h=np.exp(np.log(odd_ratio)+1.96*se)
   return odd_ratio,se,ci_l,ci_h,pvalue_fisher,chi2,pvalue_chi2,chi2_corrected,pvalue_corrected

parser = argparse.ArgumentParser(description='Uses disease and control data in tsv files and the output file to print the resutl')
parser.add_argument('-p','--patient_data', dest='pfile',type=argparse.FileType('r'),help='The file containing the patient data along with their alleles')
parser.add_argument('-c', '--control_data', dest='cfile',type=argparse.FileType('r'),help='The file containing the data from control/healthy individuals along with their alleles')
parser.add_argument('-o','--output', dest='ofile',type=argparse.FileType('w'),default='output.tsv',help='The output file to save the results in a file, default is output.tsv')
args = parser.parse_args()

dis=pd.read_csv(args.pfile,sep="\t")
con=pd.read_csv(args.cfile,sep="\t")

print con[con['Age']== ""]
print con['Age'].unique().tolist()


def run_param(dis,con):
  #parameters=['A_ag','B_ag','C_ag','DR_ag','DQ_ag']
  parameters=['A_ag','B_ag','DR_ag']
  final_pd=pd.DataFrame()
  for param in parameters:
#  param='DQ_ag'
    dis_list=[]
    dis[param+"_1"].replace(" ","",inplace=True,regex=True)
    con[param+"_1"].replace(" ","",inplace=True,regex=True)
    dis[param+"_2"].replace(" ","",inplace=True,regex=True)
    con[param+"_2"].replace(" ","",inplace=True,regex=True)
    
    dis[param+"_2"].replace("-",dis[param+"_1"],inplace=True)
    con[param+"_2"].replace("-",con[param+"_1"],inplace=True)
    dis[param+"_1"].replace("-",dis[param+"_2"],inplace=True)
    con[param+"_1"].replace("-",con[param+"_2"],inplace=True)
    
    #for i in ["_1","_2"]:
    #  dis[param+i].dropna('Blank',inplace=True)
    #  con[param+i].fillna('Blank',inplace=True)
    col1=list(dis)
    col1.remove(param+"_2")
    dis1=dis[col1].copy()
    dis1.rename(columns={param+"_1":param},inplace=True)
    col2=list(dis)
    col2.remove(param+"_1")
    dis2=dis[col2].copy()
    dis2.rename(columns={param+"_2":param},inplace=True)
    dis=dis1.append(dis2,sort=True)
    col1=list(con)
    col1.remove(param+"_2")
    con1=con[col1].copy()
    con1.rename(columns={param+"_1":param},inplace=True)
    col2=list(con)
    col2.remove(param+"_1")
    con2=con[col2].copy()
    con2.rename(columns={param+"_2":param},inplace=True)
    con=con1.append(con2,sort=True)
  con=con.drop_duplicates()
  con.dropna(subset=parameters,inplace=True)
  dis=dis.drop_duplicates()
  dis.dropna(subset=parameters,inplace=True)
  con_gp=con.groupby(parameters).size().to_frame('Healthy Count').reset_index()
  dis_gp=dis.groupby(parameters).size().to_frame('Patient Count').reset_index()
  data_pd=pd.merge(con_gp,dis_gp,on=parameters,how='outer').fillna(0)
  dsum=data_pd['Patient Count'].sum()
  csum=data_pd['Healthy Count'].sum()
  data_pd['Remaining Patients']=dsum-data_pd['Patient Count']
  data_pd['Remaining Healthy']=csum-data_pd['Healthy Count']
  data_pd['Patient frequency (%)']=100*data_pd['Patient Count']/dsum
  data_pd['Healthy frequency (%)']=100*data_pd['Healthy Count']/csum
  data_pd=data_pd.reset_index()
  data_pd[['odds ratio','se','ci lower','ci higher','pvalue_fisher','chi2','pvalue_chi2','chi2_corrected','pvalue_corrected']]=data_pd.apply(OR_CI,axis=1,result_type="expand")
  final_pd=final_pd.append(data_pd,sort=True)
  #dis.to_csv('cml_mod.tsv',sep="\t")
  #con.to_csv('control_mod.tsv',sep="\t")
  return final_pd

df_output=run_param(dis,con)
df_output['Category']='Overall'
#print df_out

#### Gender wise analysis
dis['Sex'].replace(" ","",inplace=True,regex=True)
con['Sex'].replace(" ","",inplace=True,regex=True)
##print dis['Sex'].unique()
##print con['Sex'].unique()
#df_output=pd.DataFrame()
print con['Age'].unique()
print con['Age'].unique()
for gender in ['M','F']:
  dis_gender=dis.loc[dis['Sex']==gender].copy()
  con_gender=con.loc[con['Sex']==gender].copy()
  df_gen=run_param(dis_gender,con_gender)
  df_gen['Category']=gender
  df_output=df_output.append(df_gen,ignore_index=True,sort=True)
dis_younger=dis.loc[dis['Age']<=18].copy()
con_younger=con.loc[con['Age']<=18].copy()
df_gen=run_param(dis_younger,con_younger)
df_gen['Category']='Younger <=18 years'
df_output=df_output.append(df_gen,ignore_index=True,sort=True)
dis_adults=dis.loc[dis['Age']>18].copy()
con_adults=con.loc[con['Age']>18].copy()
df_gen=run_param(dis_adults,con_adults)
df_gen['Category']='Adults >18 years'
df_output=df_output.append(df_gen,ignore_index=True,sort=True)
df_output.to_csv(args.ofile,index=False)
