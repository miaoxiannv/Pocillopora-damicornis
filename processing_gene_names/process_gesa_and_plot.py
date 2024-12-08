"""
读取go-bp,go-mf,go-cc
格式如下
Term	Term_name	Adjusted P-value
GO:0006412	translation	1.6417118478617708e-11
GO:0007049	cell cycle	0.0246210461042345
将Term和Term_name合并，并且存储为字典，加上一列Type，表示bp,mf,cc
Term和Term_name合并使用"~"连接
加上一行Count，值是每个Term的P-value的-log10
读取kegg的tsv文件，格式如下：
Term	Term_name	Adjusted P-value
ko03010	 Genetic Information Processing; Translation; Ribosome	9.011042377397073e-14
ko04915	 Organismal Systems; Endocrine system; Estrogen signaling pathway	0.00043279677427967525
ko04612	 Organismal Systems; Immune system; Antigen processing and presentation	0.0023276366422452476
ko04659	 Organismal Systems; Immune system; Th17 cell differentiation	0.0023276366422452476
ko04933	 Human Diseases; Endocrine and metabolic diseases; AGE-RAGE signaling pathway in diabetic complications	0.007612498330968531
ko04657	 Organismal Systems; Immune system; IL-17 signaling pathway	0.01596281406922652
ko05164	 Human Diseases; Infectious diseases: Viral; Influenza A	0.021337110819139034
ko05418	 Human Diseases; Cardiovascular diseases; Fluid shear stress and atherosclerosis	0.02450115935747163
ko04621	 Organismal Systems; Immune system; NOD-like receptor signaling pathway	0.03239905400617683
ko05166	 Human Diseases; Infectious diseases: Viral; HTLV-I infection	0.033448502996175
ko05169	 Human Diseases; Infectious diseases: Viral; Epstein-Barr virus infection	0.04173370814529228



加上一列Type，表示kegg
将Term分割为ko和和数字，然后用：合并，合并后的数据和Term_name合并，合并规则是~
加上一列Count，值是每个Term的P-value的-log10

"""
from pathlib import Path

import numpy as np
import pandas as pd

"""
处理GO注释文件并合并结果

Parameters:
bp_file, mf_file, cc_file: GO注释文件的路径

Returns:
DataFrame: 合并后的数据框
"""


def process_single_file(file_path, go_type):
    try:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"找不到文件: {file_path}")

        # 读取文件
        df = pd.read_csv(file_path, sep='\t')

        # 合并Term和Term_name
        df['Term'] = df['Term'] + '~' + df['Term_name']

        # 删除原Term_name列
        df = df.drop('Term_name', axis=1)

        # 添加Type列
        df['Type'] = go_type

        # 计算Count (-log10(p-value))
        df['Count'] = -np.log10(df['Adjusted P-value'])

        return df

    except Exception as e:
        raise Exception(f"处理{go_type}文件时出错: {str(e)}")


def process_go_files(bp_file, mf_file, cc_file):
    # 处理每个文件
    bp_df = process_single_file(bp_file, 'BP')
    mf_df = process_single_file(mf_file, 'MF')
    cc_df = process_single_file(cc_file, 'CC')

    # 合并所有数据框
    result_df = pd.concat([bp_df, mf_df, cc_df], axis=0)

    return result_df


def process_kegg_file(kegg_file):
    """
    处理KEGG注释文件
    
    Parameters:
    kegg_file: KEGG注释文件的路径
    
    Returns:
    DataFrame: 处理后的数据框
    """
    try:
        # 检查文件是否存在
        if not Path(kegg_file).exists():
            raise FileNotFoundError(f"找不到文件: {kegg_file}")

        # 读取文件
        df = pd.read_csv(kegg_file, sep='\t')

        # 验证必要的列是否存在
        required_columns = ['Term', 'Term_name', 'Adjusted P-value']
        if not all(col in df.columns for col in required_columns):
            raise ValueError("文件缺少必要的列：Term, Term_name, Adjusted P-value")
            
        # 处理Term列：检查格式并转换
        def process_term(term):
            if not term.startswith('ko'):
                raise ValueError(f"Term格式错误: {term}，应以'ko'开头")
            return 'ko:' + term[2:]
            
        df['Term'] = df['Term'].apply(process_term)
        
        # 去除Term_name开头的空格
        df['Term_name'] = df['Term_name'].str.strip()
        
        # 合并Term和Term_name
        df['Term'] = df['Term'] + '~' + df['Term_name']
        
        # 删除原Term_name列
        df = df.drop('Term_name', axis=1)
        
        # 添加Type列
        df['Type'] = 'KEGG'
        
        # 计算Count (-log10(p-value))，处理异常值
        df['Count'] = df['Adjusted P-value'].apply(
            lambda x: -np.log10(x) if x > 0 else 0
        )
        
        return df
        
    except Exception as e:
        raise Exception(f"处理KEGG文件时出错: {str(e)}")


def get_enrichment_dataframe(bp_file, mf_file, cc_file, kegg_file):
    """
    直接获取合并后的富集分析数据框
    
    Parameters:
    -----------
    bp_file : str
        GO生物过程文件路径
    mf_file : str
        GO分子功能文件路径
    cc_file : str
        GO细胞组分文件路径
    kegg_file : str
        KEGG通路文件路径
        
    Returns:
    --------
    pd.DataFrame
        合并后的富集分析结果数据框
    """
    try:
        # 处理GO文件
        go_result = process_go_files(bp_file, mf_file, cc_file)

        # 处理KEGG文件
        kegg_result = process_kegg_file(kegg_file)

        # 合并结果
        final_result = pd.concat([go_result, kegg_result], axis=0)

        return final_result

    except Exception as e:
        raise Exception(f"获取富集分析数据框时出错: {str(e)}")



def save_enrichment_results(df, output_path):
    """
    保存富集分析结果到txt文件，处理重复的Term
    
    Parameters:
    -----------
    df : pd.DataFrame
        富集分析结果数据框
    output_path : str
        输出文件路径
    """
    try:
        # 检查是否有重复的Term
        duplicated_terms = df[df['Term'].duplicated()]['Term'].tolist()
        if duplicated_terms:
            print(f"发现重复的Term: {duplicated_terms}")
            
            # 为重复的Term添加编号
            df['Term'] = df.groupby('Term').cumcount().astype(str) + '_' + df['Term']
            # 只为重复项添加编号（第一个保持原样）
            df['Term'] = df.apply(lambda x: x['Term'][2:] if x['Term'].startswith('0_') else x['Term'], axis=1)
            
        # 确保输出目录存在
        output_dir = Path(output_path).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存为tab分隔的txt文件
        df.to_csv(output_path, sep='\t', index=False)
        print(f"结果已保存到: {output_path}")
        
    except Exception as e:
        raise Exception(f"保存结果时出错: {str(e)}")

def process_kog_pfam_file(file_path, file_type):
    """
    处理KOG或PFAM注释文件
    
    Parameters:
    -----------
    file_path : str
        文件路径
    file_type : str
        文件类型，'KOG'或'PFAM'
        
    Returns:
    --------
    DataFrame: 处理后的数据框
    """
    try:
        # 检查文件是否存在
        if not Path(file_path).exists():
            raise FileNotFoundError(f"找不到文件: {file_path}")
            
        # 读取文件
        df = pd.read_csv(file_path, sep='\t')
        
        # 验证必要的列是否存在
        required_columns = ['Term', 'Term_name', 'Adjusted P-value']
        if not all(col in df.columns for col in required_columns):
            raise ValueError("文件缺少必要的列：Term, Term_name, Adjusted P-value")
        
        # 处理Term列：分离字母和数字，用:连接
        def process_term(term):
            # 对于KOG，将字母和数字分开并用:连接
            if file_type == 'KOG':
                letters = ''.join(filter(str.isalpha, term))
                numbers = ''.join(filter(str.isdigit, term))
                return f"{letters}:{numbers}"
            # 对于PFAM，保持原样并添加:
            else:
                return f"PFAM:{term}"
                
        df['Term'] = df['Term'].apply(process_term)
        
        # 去除Term_name的首尾空格
        df['Term_name'] = df['Term_name'].str.strip()
        
        # 合并Term和Term_name
        df['Term'] = df['Term'] + '~' + df['Term_name']
        
        # 删除原Term_name列
        df = df.drop('Term_name', axis=1)
        
        # 添加Type列
        df['Type'] = file_type
        
        # 计算Count (-log10(p-value))
        df['Count'] = df['Adjusted P-value'].apply(
            lambda x: -np.log10(x) if x > 0 else 0
        )
        
        return df
        
    except Exception as e:
        raise Exception(f"处理{file_type}文件时出错: {str(e)}")

def get_kog_pfam_dataframe(kog_file, pfam_file):
    """
    获取合并后的KOG和PFAM富集分析数据框
    
    Parameters:
    -----------
    kog_file : str
        KOG注释文件路径
    pfam_file : str
        PFAM注释文件路径
        
    Returns:
    --------
    pd.DataFrame
        合并后的富集分析结果数据框
    """
    try:
        # 处理KOG文件
        kog_result = process_kog_pfam_file(kog_file, 'KOG')
        
        # 处理PFAM文件
        pfam_result = process_kog_pfam_file(pfam_file, 'PFAM')
        
        # 合并结果
        final_result = pd.concat([kog_result, pfam_result], axis=0)
        
        return final_result
        
    except Exception as e:
        raise Exception(f"获取KOG和PFAM数据框时出错: {str(e)}")

# 使用示例：
def save_kog_pfam_results(kog_file, pfam_file, output_path):
    """
    处理并保存KOG和PFAM的富集分析结果
    
    Parameters:
    -----------
    kog_file : str
        KOG注释文件路径
    pfam_file : str
        PFAM注释���件路径
    output_path : str
        输出文件路径
    """
    try:
        # 获取合并的数据框
        df = get_kog_pfam_dataframe(kog_file, pfam_file)
        
        # 保存结果
        save_enrichment_results(df, output_path)
        
    except Exception as e:
        print(f"保存KOG和PFAM结果时出错: {str(e)}")
# 使用示例：

"""df = get_enrichment_dataframe(
    bp_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\go_bp_res.tsv",
    mf_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\go_mf_res.tsv",
    cc_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\go_cc_res.tsv",
    kegg_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\kegg_res.tsv"
)

output_path = r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\go-kegg_results.txt"
save_enrichment_results(df, output_path)



# 示例使用
save_kog_pfam_results(
    kog_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\kog_res.tsv",
    pfam_file=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\pfam_res.tsv",
    output_path=r"D:\nextcloud\pd论文\result\GSEA\Cluster0\positive\kog-pfam_results.txt"
)"""

def process_cluster_enrichment(base_dir, cluster_list, direction_list=None):
    """
    批量处理多个聚类的富集分析结果
    
    Parameters:
    -----------
    base_dir : str
        基础目录路径
    cluster_list : list
        聚类编号列表，如 ['0', '1', '2']
    direction_list : list, optional
        方向列表，如 ['positive', 'negative']，默认为 None
        如果为 None，则只处理一个方向
    
    Returns:
    --------
    dict: 处理结果的字典，包含每个聚类的处理状态
    """
    results = {}
    
    try:
        for cluster in cluster_list:
            cluster_dir = f"Cluster{cluster}"
            
            # 如果没有指定方向，则创建空列表
            if direction_list is None:
                direction_list = ['']
            
            for direction in direction_list:
                try:
                    # 构建完整路径
                    current_dir = Path(base_dir) / cluster_dir
                    if direction:
                        current_dir = current_dir / direction
                    
                    # 构建文件路径
                    file_paths = {
                        'go': {
                            'bp': current_dir / "go_bp_res.tsv",
                            'mf': current_dir / "go_mf_res.tsv",
                            'cc': current_dir / "go_cc_res.tsv"
                        },
                        'kegg': current_dir / "kegg_res.tsv",
                        'kog': current_dir / "kog_res.tsv",
                        'pfam': current_dir / "pfam_res.tsv"
                    }
                    
                    # 处理GO和KEGG结果
                    go_kegg_df = get_enrichment_dataframe(
                        bp_file=str(file_paths['go']['bp']),
                        mf_file=str(file_paths['go']['mf']),
                        cc_file=str(file_paths['go']['cc']),
                        kegg_file=str(file_paths['kegg'])
                    )
                    
                    # 处理KOG和PFAM结果
                    kog_pfam_df = get_kog_pfam_dataframe(
                        kog_file=str(file_paths['kog']),
                        pfam_file=str(file_paths['pfam'])
                    )
                    
                    # 构建输出路径
                    output_dir = current_dir
                    output_dir.mkdir(parents=True, exist_ok=True)
                    
                    # 保存结果
                    go_kegg_output = output_dir / "go-kegg_results.txt"
                    kog_pfam_output = output_dir / "kog-pfam_results.txt"
                    
                    save_enrichment_results(go_kegg_df, str(go_kegg_output))
                    save_enrichment_results(kog_pfam_df, str(kog_pfam_output))
                    
                    # 记录处理状态
                    result_key = f"Cluster{cluster}"
                    if direction:
                        result_key += f"_{direction}"
                    results[result_key] = "成功"
                    
                except Exception as e:
                    result_key = f"Cluster{cluster}"
                    if direction:
                        result_key += f"_{direction}"
                    results[result_key] = f"失败: {str(e)}"
                    
    except Exception as e:
        print(f"处理过程中出错: {str(e)}")
        
    return results

# 使用示例：
if __name__ == "__main__":
    # 基础目录
    base_dir = r"D:\nextcloud\pd论文\result\GSEA"
    
    # 要处理的聚类列表
    clusters = ['0', '1', '2', '3', '6', '7']
    
    # 方向列表
    directions = ['positive', 'negative']
    
    # 执行批量处理
    results = process_cluster_enrichment(base_dir, clusters, directions)
    
    # 打印处理结果
    print("\n处理结果汇总:")
    for cluster, status in results.items():
        print(f"{cluster}: {status}")