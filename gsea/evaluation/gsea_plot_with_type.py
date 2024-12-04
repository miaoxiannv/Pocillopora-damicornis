import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from typing import Dict, Optional
import seaborn as sns
import matplotlib.colors as mcolors
import numpy as np

'''
sample tsv file:

Term    -log-pv    Type
GO:0000001   1.0    Biological Process
GO:0000002   2.0    Biological Process
GO:0000003   3.0    Biological Process
GO:0000004   4.0    Biological Process
GO:0000005   5.0    Biological Process
GO:0000006   6.0    Biological Process
GO:0000007   7.0    Cellular Component
GO:0000008   8.0    Cellular Component
GO:0000009   9.0    Cellular Component
GO:0000010   10.0   Cellular Component
GO:0000011   11.0   Cellular Component
GO:0000012   12.0   Cellular Component
'''

class EnrichmentPlot:
    def __init__(self, 
                 color_scheme: Optional[Dict[str, str]] = None,
                 figsize: tuple = (12, 10),
                 font_size: dict = None,
                 font_family: str = 'Times New Roman'
                 ):
        """
        Initialize the EnrichmentPlot class.
        
        Args:
            color_scheme (dict): Custom color mapping for different types
            figsize (tuple): Figure size for the plot
            font_size (dict): Dictionary containing font sizes for different elements
                             Keys: 'title', 'label', 'tick', 'legend'
        """
        #set default front family
        if font_family:
            plt.rcParams['font.family'] = font_family
        else:
            plt.rcParams['font.family'] = 'Times New Roman'
        self.figsize = figsize
        self.color_scheme = None
        self.df = None
        # Set default font sizes
        self.font_size = {
            'title': 14,
            'label': 12,
            'tick': 10,
            'legend': 10
        }
        # Update font sizes based on user input
        if font_size:
            self.font_size.update(font_size)

    def _generate_color_scheme(self, types: list) -> None:
        """
        Generate a color scheme for the unique types in the data.
        Uses seaborn's color palette for distinct colors.
        
        Args:
            types (list): List of unique types in the data
        """
        if len(types) > 6:
            raise ValueError("Maximum 6 different types are supported")
        
        # Generate distinct colors using seaborn's color palette
        colors = sns.color_palette("deep", n_colors=len(types))
        # Convert RGB to hex colors
        hex_colors = [mcolors.rgb2hex(color) for color in colors]
        # Create color mapping
        self.color_scheme = dict(zip(types, hex_colors))
        
    #load data from tsv file or pandas dataframe
    def load_data(self, file_path, 
                  transformed: bool = False,
                  trans_para_dic: dict = None,
                  top_n: Optional[int] = None
                  ) -> None:
        """
        Load data from TSV file.
        
        Args:
            file_path (str): Path to the TSV file
            transformed (bool): Whether the data is already transformed
            trans_para_dic (dict): Parameters for transform_pvalue
            top_n (int, optional): Number of top terms to keep for each type
        """
        if isinstance(file_path, pd.DataFrame):
            self.df = file_path
        else:
            try:
                self.df = pd.read_csv(file_path, sep='\t', header=0)
                if not transformed:
                    if trans_para_dic is not None:
                        self.transform_pvalue(**trans_para_dic)
                else:
                    required_columns = {'Term', '-log-pv', 'Type'}
                    if not all(col in self.df.columns for col in required_columns):
                        raise ValueError(f"TSV file must contain columns: {required_columns}")
                    
                    # change column "-log-pv" to float
                    self.df['-log-pv'] = self.df['-log-pv'].astype(float)
                    
                    # Generate color scheme based on unique types
                    unique_types = sorted(self.df['Type'].unique())
                    self._generate_color_scheme(unique_types)
                
                if top_n is not None:
                    # 对每个类型选择前N个值
                    self.df = (self.df.groupby('Type')
                              .apply(lambda x: x.nlargest(top_n, '-log-pv'))
                              .reset_index(drop=True))
                
            except Exception as e:
                print(f"Error loading data: {str(e)}")
                raise

    def sort_data(self, ascending: bool = True) -> None:
        """
        Sort the data by -log-pv.
        
        Args:
            ascending (bool): Sort order
        """
        # first sort by "Type" and then by "-log-pv"
        if self.df is not None:
            self.df = self.df.sort_values(by=['Type', '-log-pv'], ascending=[True, ascending])

    def create_plot(self, 
                   title: str = 'Enrichment Analysis',
                   xlabel: str = '-log(p-value)',
                   ylabel: str = 'Term') -> None:
        """
        Create and display the horizontal bar plot.
        
        Args:
            title (str): Plot title
            xlabel (str): X-axis label
            ylabel (str): Y-axis label
        """
        if self.df is None:
            raise ValueError("No data loaded. Please load data first using load_data()")

        # Create figure
        plt.figure(figsize=self.figsize)

        # Create horizontal bar plot
        bars = plt.barh(y=self.df['Term'], 
                       width=self.df['-log-pv'],
                       color=[self.color_scheme.get(t, '#000000') for t in self.df['Type']])

        # Apply font size settings
        plt.title(title, fontsize=self.font_size['title'])
        plt.xlabel(xlabel, fontsize=self.font_size['label'])
        plt.ylabel(ylabel, fontsize=self.font_size['label'])
        plt.xticks(fontsize=self.font_size['tick'])
        plt.yticks(fontsize=self.font_size['tick'])

        # Add legend with custom font size
        legend_elements = [
            Patch(facecolor=color, label=type_)
            for type_, color in self.color_scheme.items()
            if type_ in self.df['Type'].unique()
        ]
        plt.legend(handles=legend_elements, 
                  title='Type', 
                  bbox_to_anchor=(1.05, 1), 
                  loc='upper left',
                  fontsize=self.font_size['legend'],
                  title_fontsize=self.font_size['legend'])

        # Adjust layout
        plt.tight_layout()

        # Show plot
        plt.show()

    def save_plot(self, output_path: str, dpi: int = 300) -> None:
        """
        Save the plot to a file.
        
        Args:
            output_path (str): Path where to save the plot
            dpi (int): DPI for the output image
        """
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')

    def transform_pvalue(self, pvalue_column: str = 'Adjusted P-value', 
                        category_column: str = 'category') -> None:
        """
        转换 P-value 为 -log(P-value)，并调整数据格式以适配绘图需求
        
        Args:
            pvalue_column (str): P-value 列的名称
            category_column (str): 分类列的名称
        """
        if self.df is None:
            raise ValueError("请先加载数据")
            
        # 计算 -log10(p-value)
        self.df['-log-pv'] = -np.log10(self.df[pvalue_column])
        
        # 重命名 category 列为 Type
        self.df = self.df.rename(columns={category_column: 'Type'})
        
        # 确保有必要的列
        if 'Term_name' in self.df.columns:
            #if Term in self.df.columns, rename it to Term_id
            if 'Term' in self.df.columns:
                self.df = self.df.rename(columns={'Term': 'Term_id'})
            # reanme Term_name to Term
            self.df = self.df.rename(columns={'Term_name': 'Term'})
        elif 'Term' not in self.df.columns:
            self.df['Term'] = self.df.index
            
        # 生成颜色方案
        unique_types = sorted(self.df['Type'].unique())
        self._generate_color_scheme(unique_types)

# Example usage:
if __name__ == "__main__":
    # 自定义字体大小
    font_sizes = {
        'title': 20,
        'label': 16,
        'tick': 16,
        'legend': 12
    }

    # 初始化绘图器时设置字体大小
    plotter = EnrichmentPlot(figsize=(15, 12), 
                             font_size=font_sizes, 
                             font_family='Arial')
    
    # Load data
    plotter.load_data(r'V:\DATA\nextcloud\pd论文\data\test\gsea\plot\gsea.tsv', 
                      transformed=False,
                      trans_para_dic={
                           'pvalue_column': 'Adjusted P-value',
                           'category_column': 'category'
                      },
                      top_n=15)

    
    # Optionally sort data
    plotter.sort_data(ascending=True)
    
    # Create and show plot
    plotter.create_plot()
    
    # Optionally save plot
    # plotter.save_plot('output.png')