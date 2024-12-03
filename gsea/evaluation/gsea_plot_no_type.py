import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional
import numpy as np

class SimpleEnrichmentPlot:
    def __init__(self, figsize: tuple = (12, 10),
                 font_family: str = 'Times New Roman'
                 ):
        """
        Initialize the SimpleEnrichmentPlot class.
        
        Args:
            figsize (tuple): Figure size for the plot
        """
        self.figsize = figsize
        self.df = None
        
        # 设置全局字体
        if font_family:
            plt.rcParams['font.family'] = font_family
        else:
            plt.rcParams['font.family'] = ['Times New Roman', 'DejaVu Serif']
        
        # 添加预定义的颜色方案
        self.color_schemes = {
            'blue': '#1f77b4',      # 经典蓝色
            'red': '#d62728',       # 醒目红色
            'green': '#2ca02c',     # 森林绿
            'purple': '#9467bd',    # 优雅紫色
            'orange': '#ff7f0e',    # 活力橙色
            'brown': '#8c564b',     # 深棕色
            'pink': '#e377c2',      # 粉红色
            'gray': '#7f7f7f'       # 中性灰色
        }

    def load_data(self, file_path, 
                  transformed: bool = False,
                  pvalue_column: str = 'Adjusted P-value',
                  top_n: Optional[int] = None
                  ) -> None:
        """
        Load data from TSV file.
        
        Args:
            file_path (str): Path to the TSV file
            transformed (bool): Whether the data is already transformed
            pvalue_column (str): Column name for p-values
            top_n (int, optional): Number of top terms to keep
        """
        if isinstance(file_path, pd.DataFrame):
            self.df = file_path
        else:
            try:
                self.df = pd.read_csv(file_path, sep='\t', header=0)
                
                if not transformed:
                    # 计算 -log10(p-value)
                    self.df['-log-pv'] = -np.log10(self.df[pvalue_column])
                    
                    # 确保有Term列
                    if 'Term_name' in self.df.columns:
                        if 'Term' in self.df.columns:
                            self.df = self.df.rename(columns={'Term': 'Term_id'})
                        self.df = self.df.rename(columns={'Term_name': 'Term'})
                    elif 'Term' not in self.df.columns:
                        self.df['Term'] = self.df.index
                
                if top_n is not None:
                    # 选择前N个值
                    self.df = self.df.nlargest(top_n, '-log-pv')
                
            except Exception as e:
                print(f"Error loading data: {str(e)}")
                raise

    def sort_data(self, ascending: bool = True) -> None:
        """
        Sort the data by -log-pv.
        
        Args:
            ascending (bool): Sort order
        """
        if self.df is not None:
            self.df = self.df.sort_values(by='-log-pv', ascending=ascending)

    def create_plot(self, 
                   title: str = 'Enrichment Analysis',
                   xlabel: str = '-log(p-value)',
                   ylabel: str = 'Term',
                   color: str = 'blue',
                   title_fontsize: int = 12,
                   label_fontsize: int = 10,
                   tick_fontsize: int = 10,
                   ) -> None:
        """
        Create and display the horizontal bar plot.
        
        Args:
            title (str): Plot title
            xlabel (str): X-axis label
            ylabel (str): Y-axis label
            color (str): Bar color, can be one of: 'blue', 'red', 'green', 'purple',
                        'orange', 'brown', 'pink', 'gray' or a hex color code
            title_fontsize (int): Font size for the title
            label_fontsize (int): Font size for axis labels
            tick_fontsize (int): Font size for tick labels
        """
        if self.df is None:
            raise ValueError("No data loaded. Please load data first using load_data()")

        # 获取颜色值
        plot_color = self.color_schemes.get(color, color)

        # Create figure
        plt.figure(figsize=self.figsize)

        # Create horizontal bar plot
        plt.barh(y=self.df['Term'], width=self.df['-log-pv'], color=plot_color)

        # Customize the plot with font sizes
        plt.xlabel(xlabel, fontsize=label_fontsize)
        plt.ylabel(ylabel, fontsize=label_fontsize)
        plt.title(title, fontsize=title_fontsize)
        
        # 设置刻度标签的字体大小
        plt.xticks(fontsize=tick_fontsize)
        plt.yticks(fontsize=tick_fontsize)

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

# Example usage:
if __name__ == "__main__":
    # Initialize the plotter
    plotter = SimpleEnrichmentPlot(figsize=(15, 12),
                                   font_family='Arial')
    
    # Load data
    plotter.load_data(r'V:\DATA\nextcloud\pd论文\data\test\gsea\plot\gsea_no_type.tsv', 
                      transformed=False,
                      pvalue_column='Adjusted P-value',
                      top_n=10)
    
    # Optionally sort data
    plotter.sort_data(ascending=True)

    # 自定义字体大小
    font_sizes = {
        'title': 20,
        'label': 16,
        'tick': 16
    }
    
    # Create and show plot
    plotter.create_plot(
        title_fontsize=font_sizes['title'],
        label_fontsize=font_sizes['label'],
        tick_fontsize=font_sizes['tick']
    )
    
    # Optionally save plot
    # plotter.save_plot('output.png')