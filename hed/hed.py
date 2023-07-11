__Author__ = "Min Li"
# !usr/bin/env python3.7
# coding:utf8

import argparse
import pandas as pd
import pysam as ps
import os
import sys
import numpy as np

__description__="""Package zhuxl's HED. 
                   Hed refers to the Grantham distance between two alleles 
                   located at the same locus on a pair of homologous chromosomes."""


# parameter
def opt():
    args = argparse.ArgumentParser(description=__description__)
    args.add_argument('-t', '--hlatype', dest='hlatype', required=True, 
                     help="hla typing results,[.tsv].")
    args.add_argument("-o", "--output", dest="output", required=True,
                     help="output dir.[string]")
    args.add_argument("-s", "--sample", dest="sample", required=True,
                     help="sample id[string]")
    args.add_argument("-f", "--hlaFasta", dest="hlaFasta", required=True,
                     help="db.fa,[db]")
    args.add_argument("-d", "--matrix", dest="matrix", required=True,
                     help="matrix file,to CalculatePairwiseDistances.pl's -d parameter[.cnv]")
    return args.parse_args()
    
# read hlatype
def read_hlatype(hlatype):
    hlatype_df = pd.read_csv(hlatype, sep="\t")
    columns = hlatype_df.columns
    if "Unnamed: 0" in columns:
        hlatype_df.rename(columns={"Unnamed: 0": "HLA"},inplace=True)
        if "QC" not in columns:
            hlatype_df["QC"]="Failed"
    else:
        # header ABC is lose
        if hlatype_df["HLA"][0]!=0:
            hlatype_df["A1"]=hlatype_df["HLA"]
            hlatype_df["HLA"]=0
        if "QC" not in columns:
            hlatype_df["QC"]="Failed"
        # header exist,value not exist, QC not exit !!!
        hlatype_df = hlatype_df[["HLA", "A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective", "QC"]]
        # nas = np.where((hlatype_df=="Qualification")|(hlatype_df=="Failed"))[1] # python2.7.13ok, python3 no
        nas = np.where((hlatype_df.loc[0]=="Qualification")|(hlatype_df.loc[0]=="Failed"))[0] # python2.7.13ok, python3 ok
        nas0 = nas[0]
        print(nas0) 
        if nas0<9:
            # Location Reads,Objective and QC
            for i in range(nas0 , nas0-3, -1):
                hlatype_df[columns[i+9-nas0]] = hlatype_df[columns[i]]
                if i<7:
                    hlatype_df[columns[i]] = np.nan
            # hla list
            hlalist = [hlatype_df[c][0] for c in columns[1:7] if str(hlatype_df[c][0])!="nan"]
            hlatype_df[["A1", "B1", "C1", "A2", "B2", "C2"]] = np.nan
            hladic = {"A":"A1", "B":"B1", "C":"C1"}
            for hla in hlalist:
                for i in ["A", "B", "C"]:
                    if i in hla:
                        hlatype_df[hladic[i]]=hla
                        hladic[i] = i+"2"
                        break
    return hlatype_df
    
# change_hlatype
def change_hlatype(hlatype_df, wb_output):
    if hlatype_df["QC"][0] == "Failed":
        print("hla typing result's QC is Failed")
    df = hlatype_df.copy(deep=True)
    for i in ["A", "B", "C"]:
        hla1 = df[i+"1"][0]
        hla2 = df[i+"2"][0]
        if str(hla1)=="nan" or str(hla2)=="nan":
            df[i+"1"] = np.nan
            df[i+"2"] = np.nan
    NAs = df.columns[np.where(df.isnull())[1]]
    NA_length = len(NAs)
    # NA_length {6:no hed, all hom:no hed, no all hom:hed}
    if NA_length>=6:
        print("all hla typing result is NA")
        return (NAs, None) # no hed
    heds = []
    for i in ["A", "B", "C"]:
        hla1 = df[i+"1"][0]
        hla2 = df[i+"2"][0]
        if i+"1" not in NAs and hla1!=hla2:
            heds.append("hla_"+hla1.replace("*","_").replace(":","_").lower())
            heds.append("hla_"+hla2.replace("*","_").replace(":","_").lower())
    if len(heds)>0:
        heds = pd.DataFrame(heds)
        wb_hlatype = os.path.join(wb_output, "hla_type_change.txt")
        heds.to_csv(wb_hlatype, sep="\t", header=False, index=False)
        return NAs, wb_hlatype
    else:
        print("hla typing result is all hom")
        return NAs, None


# merge db.tsv
def mergedb(sample, wb_output, wb_hlatype, hlaFasta):
    df = pd.read_csv(wb_hlatype, sep="hla_", header=None, engine='python')
    df[0] = df[1].apply(lambda x: x.replace("_","").upper())
    df = df[[0]]
    db = pd.read_csv(hlaFasta, sep="\t", header=None)
    db[0] = db[0].apply(lambda x: x.replace(">",""))
    db=db[[0,1]]
    df = pd.merge(df, db, on=0, how="left")
    df[0] = df[0].apply(lambda x: ">{}_{}".format(sample,x))
    hlafa = os.path.join(wb_output, "{sampleid}.HLA.fa".format(sampleid=sample))
    df.to_csv(hlafa, sep="\n", index=False, header=False)
    return hlafa


def CalculatePairwiseDistances(sample, wb_output, matrix, hlafa):
    wb_hlahed = os.path.dirname(os.path.abspath(__file__))
    pdl = os.path.join(wb_output, "{sampleid}.HLA.fa_PairwiseDistanceList.txt".format(sampleid=sample))
    mod = os.path.join(wb_output, "{sampleid}.HLA.fa_PairwiseDistanceList.txt.mod".format(sampleid=sample))
    t_cmd = "/usr/local/bin/perl {wb_hlahed}/bin/CalculatePairwiseDistances.pl "\
            "-d {aamatrix} "\
            "-f {wb_hlafa} && "\
            "sed '1d;s/ _ /\t/' {pdl} "\
            "> {s_mod}".format(wb_hlahed = wb_hlahed,
                       aamatrix=matrix,
                      wb_hlafa=hlafa,
                     pdl=pdl,
                    s_mod=mod)
    print(t_cmd)
    os.system(t_cmd)
    return mod
    
    
def generate_result(sample, hlatype_df, wb_output, mod, NAs):
    dfleft = hlatype_df
    if mod is not None:
        df = pd.read_csv(mod, names=['Allele1', 'Allele2', 'Distance'], sep='\t')
        df['Sample1'] = sample
        df['Sample2'] = sample
        df[['Allele1']] = df.Allele1.str.replace(sample+"_", "")
        df[['Allele2']] = df.Allele2.str.replace(sample+"_", "")    

        df = df[(df.Sample1 == df.Sample2) & (df.Allele1.str[:1] == df.Allele2.str[:1])]
        df['Allele'] = df["Allele1"].apply(lambda x:"HED-"+x[:1])
        df = df [['Allele', 'Distance']]
        df = df.set_index('Allele').T
        df.reset_index(drop=True, inplace=True)
        df = pd.concat([dfleft, df], axis=1)
        total = 0
        count = 0
        
        for i in ["A", "B", "C"]:
            if i+"1" in NAs:
                df["HED-"+i] = np.nan
            elif "HED-"+i in list(df.columns):
                df["HED-"+i] = round(df["HED-"+i][0], 2)
                total += df["HED-"+i][0]
                count += 1
            else:
                df["HED-"+i] = 0
                count += 1
        df["HED-TOTAL"] = round(total/count, 2)
    else:
        df = dfleft
        for i in ["A", "B", "C"]:
            if i+"1" in NAs:
                df["HED-"+i] = np.nan
            else:
                df["HED-"+i] = 0
        if len(NAs)>=6:
            df["HED-TOTAL"] = np.nan
        else:
            df["HED-TOTAL"] = 0
    df.fillna("NA", inplace=True)
    return df


def main():
    # get parameter
    args = opt()
    sample = args.sample
    output = args.output
    hlatype = args.hlatype
    hlaFasta = args.hlaFasta
    matrix = args.matrix
    columns = ["HLA", "A1", "A2", "B1", "B2", "C1", "C2", "Reads", "Objective", "QC", "HED-A", "HED-B", "HED-C", "HED-TOTAL"]
    fail = os.path.join(output, "{sampleid}.hed.FAIL".format(sampleid=sample))
    succ = os.path.join(output, "{sampleid}.hed.SUCCESS".format(sampleid=sample))
    current_work_dir = os.path.split(os.path.realpath(__file__))[0]
    t_cmd = "python {work_dir}/hed.py "\
            "-t {hla_typing} "\
            "-o {hed_dir} "\
            "-s {sampleid} "\
            "-f {work_dir}/db/db.fa "\
            "-d {work_dir}/bin/HED/AAdistMatrix_Grantham.cnv".format(work_dir=current_work_dir,
                       hla_typing=hlatype,
                      hed_dir=output,
                     sampleid=sample,)
    print(t_cmd)
    if not os.path.exists(output):
        os.mkdir(output)	
    try:
        hlatype_df = read_hlatype(hlatype)
    except:
        df = pd.DataFrame([{i:"Failed" for i in columns}])
        os.system("touch {0}".format(fail))
        os.system("rm -rf {0}".format(succ))
    else:
        try:
            # change hla format
            NAs, wb_hlatype = change_hlatype(hlatype_df, output)
            # merge hla db sequence
            hlafa = None
            if wb_hlatype is not None:
                hlafa = mergedb(sample, output, wb_hlatype, hlaFasta)
            mod = None    
            if hlafa is not None:
                mod = CalculatePairwiseDistances(sample, output, matrix, hlafa)
            df = generate_result(sample, hlatype_df, output, mod, NAs)
            os.system("touch {0}".format(succ))
            os.system("rm -rf {0}".format(fail))
        except:
            df =  hlatype_df
            df.fillna("NA", inplace=True)
            for i in ["A", "B", "C", "TOTAL"]:
                df["HED-"+i]="Failed"
            os.system("touch {0}".format(fail))
            os.system("rm -rf {0}".format(succ))
    finally:
        hedtsv = os.path.join(output, "{sampleid}.hed.tsv".format(sampleid=sample))
        df.to_csv(hedtsv, sep='\t', index=False, columns=columns)

if __name__ == "__main__":
    main()
