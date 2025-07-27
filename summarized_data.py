import pandas as pd
import argparse

def parse_eas_file(filename):
    """ 解析EAS文件，提取基因名、表型、频率和个数 """
    with open(filename, "r", encoding="utf-8") as file:
        data = {}
        gene = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                gene = line[1:].split("&")[0]  # 提取基因名
                data[gene] = {}
            elif gene and line:
                parts = line.split("\t")
                if len(parts) == 3:
                    phenotype, frequency, count = parts
                    data[gene][phenotype] = (float(frequency), int(count))
    return data

def parse_mixture_file(filename):
    """ 解析mixture文件，提取基因名和表型 """
    with open(filename, "r", encoding="utf-8") as file:
        data = []
        gene = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                gene = line[1:].split("&")[0]  # 提取基因名
            elif gene and line:
                parts = line.split("\t")
                if len(parts) == 3:
                    phenotype = parts[0]
                    data.append((gene, phenotype))
    return data

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(description="匹配基因表型频率和个数")
    parser.add_argument("-e", "--eas", required=True, help="EAS文件路径")
    parser.add_argument("-m", "--mixture", required=True, help="mixture文件路径")
    parser.add_argument("-o", "--output", required=True, help="输出Excel文件路径")

    args = parser.parse_args()

    # 读取数据
    eas_data = parse_eas_file(args.eas)
    mixture_data = parse_mixture_file(args.mixture)

    # 匹配数据
    result = []
    for gene, phenotype in mixture_data:
        if gene in eas_data and phenotype in eas_data[gene]:
            freq, count = eas_data[gene][phenotype]
            result.append([gene, phenotype, freq, count])
        else:
            result.append([gene, phenotype, "Not Found", "Not Found"])

    # 转换为DataFrame
    df = pd.DataFrame(result, columns=["Gene", "Phenotype", "Frequency", "Count"])

    # 保存为Excel
    df.to_excel(args.output, index=False)

    print(f"匹配完成，结果已保存到 {args.output}")

if __name__ == "__main__":
    main()

