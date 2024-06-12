from pysam import VariantFile
import numpy as np
import pandas as pd
import os

vcf_path = 'BayesInteGration-main/Data/vcfFile/'
save_path = 'BayesInteGration-main/Data/csvFile/'
chro = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']

vcf_list = os.listdir(vcf_path)
for vcf_name in vcf_list:
    key = 'pbsv.vcf'
    if key in vcf_name:
        print(vcf_name)
        
        vcf = VariantFile(vcf_path+vcf_name,"r")
        arr = []
        
        for rec in vcf:
            if rec.chrom in chro:
                method = 'pbsv'
                
                SVtype = rec.info['SVTYPE']
                if SVtype =='INS' or SVtype == 'DEL' or SVtype == 'DUP' or SVtype == 'INV':
                    SVlen = abs(rec.info['SVLEN'][0])
                else:
                    continue
                
                print(method, rec.chrom, rec.start+1, rec.stop, SVlen, SVtype, rec.qual)
                arr.append([method, rec.chrom, rec.start+1, rec.stop, SVlen, SVtype, rec.qual])
        
        if not os.path.exists(save_path): os.mkdir(save_path)

        df = pd.DataFrame(arr, columns=['Method','chro','start','end','length','type','qual'])
        df.to_csv(save_path + key + ".csv", index = False)