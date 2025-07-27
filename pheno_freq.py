import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf")
parser.add_argument("-p", "--population")
args = parser.parse_args()

file_name = args.vcf
pop_name = args.population
#读取VCF和人群分组文件
pol_list = pop_name.split("/")[-1] #分组信息文件名
df = pd.read_table(file_name)
pop = pd.read_table(pop_name)

#将VCF文件表头信息添加到人群分组中
header = ["FORMAT", "INFO", "FILTER", "QUAL", "ALT", "REF", "ID", "POS", "#CHROM"]
pop_sample = pop['sample'].values.tolist()
for i in header:
    pop_sample.insert(0,i)

#从总VCF文件中摘取目标群体及根据VCF文件ID列进行分组
pop_genotype = df[pop_sample]
grouped_df = pop_genotype.groupby("ID")

#定义表型频率计算函数
def freq(x):
    df = x.transpose()
    df = df.reset_index()
    array = np.array(df)
    list = array.tolist()
    phenotypelist = []
    for i in list[9:len(list)]:
        phenotype = ""
        for genotype in i[1:len(i)]:
            genotype = genotype.split("|")
            genotype.sort()
            phenotype = phenotype + "_" + genotype[0] + "|" + genotype[1]
        phenotypelist.append(phenotype)
    out = set(phenotypelist)
    for allele in out:
        name = allele
        count = phenotypelist.count(allele)
        freq = count/len(phenotypelist)
        print(f"{name}\t{freq}\t{count}")

for name, group in grouped_df:
    print(">" + name + "&" + pol_list)
    freq(group) 
