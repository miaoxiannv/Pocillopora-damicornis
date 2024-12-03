
import pandas as pd
import matplotlib.pyplot as plt

# Create sample data
data = {
    'Term': [
        'Immune response',
        'Defense response to bacterium',
        'Cell chemotaxis',
        'Cell adhesion',
        'Complement activation',
        'Cytokine-mediated signaling pathway',
        'Phagocytosis, engulfment',
        'Negative regulation of JAK-STAT cascade',
        'Epoxygenase P450 pathway',
        'Chemokine-mediated signaling pathway',
        'Negative regulation of leukocyte apoptotic process',
        'B cell receptor signaling pathway',
        'Cellular response to tumor necrosis factor',
        'Positive regulation of chemotaxis',
        'Positive regulation of angiogenesis',
        'Collagen fibril organization',
        'Positive regulation of homotypic cell-cell adhesion',
        'Regulation of cell projection assembly',
        'Prostate epithelial cord elongation',
        'Bile acid catabolic process',
        'Cellular response to drug',
        'Glycosaminoglycan metabolic process'
    ],
    'Count': [
        20, 12, 10, 18, 8, 8, 5, 5, 4, 6, 4, 5, 6, 3, 6, 4, 20, 15, 10, 9, 8, 4
    ],
    'Type': [
        'BP', 'BP', 'BP', 'BP', 'BP',
        'MF', 'MF', 'MF', 'MF', 'MF', 'MF', 'MF',
        'CC', 'CC', 'CC', 'CC',
        'KEGG', 'KEGG', 'KEGG', 'KEGG', 'KEGG', 'KEGG'
    ]
}

# Create DataFrame
df = pd.DataFrame(data)

# Set figure size
plt.figure(figsize=(12, 10))

# Create color mapping
colors = {'BP': '#FF6B6B', 'MF': '#4ECDC4', 'CC': '#45B7D1', 'KEGG': '#2C3E50'}

# Create horizontal bar plot
bars = plt.barh(y=df['Term'], width=df['Count'], 
                color=[colors[t] for t in df['Type']])

# Customize the plot
plt.xlabel('Count')
plt.ylabel('Term')
plt.title('Biological Process Analysis')

# Add legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=colors[key], label=key) for key in colors]
plt.legend(handles=legend_elements, title='Type', bbox_to_anchor=(1.05, 1), loc='upper left')

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Show the plot
plt.show()






