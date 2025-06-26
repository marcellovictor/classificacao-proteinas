# Análise Exploratória de Proteínas usando BioPython e dados do UniProt

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import ExPASy
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import warnings
import time

warnings.filterwarnings('ignore')

print("=== ANÁLISE EXPLORATÓRIA DE PROTEÍNAS - UNIPROT ===")
print("Coletando dados do UniProt usando BioPython...")

def extrair_caracteristicas_proteina(record):
    """
    Extrai características bioquímicas de um record do SwissProt
    """
    try:
        sequencia = str(record.sequence)
        
        if len(sequencia) == 0:
            return None
            
        protein_analysis = ProteinAnalysis(sequencia)
        
        # Características básicas
        comprimento = len(sequencia)
        massa_molecular = protein_analysis.molecular_weight()
        
        # Ponto isoelétrico
        try:
            ponto_isoeletrico = protein_analysis.isoelectric_point()
        except:
            ponto_isoeletrico = np.nan
            
        # Composição de aminoácidos
        composicao = protein_analysis.get_amino_acids_percent()
        
        # Hidrofobicidade
        try:
            hidrofobicidade = protein_analysis.gravy()
        except:
            hidrofobicidade = np.nan
            
        # Instabilidade
        try:
            instabilidade = protein_analysis.instability_index()
        except:
            instabilidade = np.nan
            
        # Estrutura secundária
        try:
            estrutura_sec = protein_analysis.secondary_structure_fraction()
            helix_fraction = estrutura_sec[0]
            turn_fraction = estrutura_sec[1] 
            sheet_fraction = estrutura_sec[2]
        except:
            helix_fraction = turn_fraction = sheet_fraction = np.nan
            
        # Proporção de aminoácidos hidrofóbicos
        aminoacidos_hidrofobicos = ['A', 'I', 'L', 'M', 'F', 'W', 'Y', 'V']
        prop_hidrofobicos = sum([composicao.get(aa, 0) for aa in aminoacidos_hidrofobicos])
        
        # Proporção de aminoácidos carregados
        aminoacidos_positivos = ['K', 'R', 'H']
        aminoacidos_negativos = ['D', 'E']
        prop_positivos = sum([composicao.get(aa, 0) for aa in aminoacidos_positivos])
        prop_negativos = sum([composicao.get(aa, 0) for aa in aminoacidos_negativos])
        
        # Classificação por organismo
        organismo = str(record.organism)
        if 'Homo sapiens' in organismo:
            classe = 'Humana'
        elif 'Escherichia coli' in organismo:
            classe = 'Bacteria'
        elif 'Saccharomyces cerevisiae' in organismo:
            classe = 'Levedura'
        elif 'Mus musculus' in organismo:
            classe = 'Camundongo'
        else:
            classe = 'Outros'
            
        return {
            'ID': record.entry_name,
            'Comprimento_Sequencia': comprimento,
            'Massa_Molecular': massa_molecular,
            'Ponto_Isoeletrico': ponto_isoeletrico,
            'Hidrofobicidade': hidrofobicidade,
            'Instabilidade': instabilidade,
            'Helix_Fraction': helix_fraction,
            'Turn_Fraction': turn_fraction,
            'Sheet_Fraction': sheet_fraction,
            'Prop_Hidrofobicos': prop_hidrofobicos,
            'Prop_Positivos': prop_positivos,
            'Prop_Negativos': prop_negativos,
            'Carga_Liquida': prop_positivos - prop_negativos,
            'Classe': classe,
            'Organismo': organismo
        }
        
    except Exception as e:
        print(f"Erro ao processar proteína: {e}")
        return None

# Lista de IDs de proteínas do UniProt
protein_ids = [
    # Proteínas humanas
    'P04637',  # p53
    'P01308',  # Insulina
    'P02768',  # Albumina
    'P01009',  # Alpha-1-antitrypsin
    'P00748',  # Fator XII
    'P68871',  # Hemoglobina
    'P02647',  # Apolipoprotein A-I
    
    # Proteínas de E. coli
    'P0A6F5',  # Carbamoyl-phosphate synthase
    'P0AG63',  # Elongation factor Tu
    'P0A796',  # ATP synthase
    'P0A6Y8',  # Chaperone protein DnaK
    
    # Proteínas de levedura
    'P00924',  # Enolase
    'P00330',  # Alcohol dehydrogenase
    'P00817',  # Pyruvate kinase
    
    # Proteínas de camundongo
    'P01942',  # Hemoglobina alpha
    'P07724',  # Serum albumin
]

print(f"Coletando dados de {len(protein_ids)} proteínas...")

# Coletar dados das proteínas
dados_proteinas = []
proteinas_processadas = 0

for i, protein_id in enumerate(protein_ids):
    try:
        print(f"Processando proteína {i+1}/{len(protein_ids)}: {protein_id}")
        
        handle = ExPASy.get_sprot_raw(protein_id)
        record = SwissProt.read(handle)
        handle.close()
        
        caracteristicas = extrair_caracteristicas_proteina(record)
        
        if caracteristicas:
            dados_proteinas.append(caracteristicas)
            proteinas_processadas += 1
            
        time.sleep(0.5)
        
    except Exception as e:
        print(f"Erro ao processar {protein_id}: {e}")
        continue

print(f"\nTotal de proteínas processadas: {proteinas_processadas}")

# Criar DataFrame
df = pd.DataFrame(dados_proteinas)

if len(df) == 0:
    print("Erro: Nenhuma proteína foi processada.")
    exit()

print(f"\nDataFrame criado com {len(df)} proteínas e {len(df.columns)} colunas")

# Configurar display
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

print("\n=== ANÁLISE EXPLORATÓRIA DOS DADOS ===")

print("\nPrimeiras linhas do dataset:")
print(df.head())

print(f"\nShape do dataset: {df.shape}")
print(f"\nColunas: {df.columns.tolist()}")
print(f"\nTipos de dados:\n{df.dtypes}")
print(f"\nEstatísticas descritivas:\n{df.describe()}")

# Classes disponíveis
lista_classes = df['Classe'].unique()
print(f"\nClasses encontradas: {lista_classes}")
print(f"\nDistribuição por classe:\n{df['Classe'].value_counts()}")

# Valores nulos
print(f"\nValores nulos por coluna:\n{df.isna().sum()}")

# Colunas numéricas
colunas_numericas = df.select_dtypes(include=[np.number]).columns.tolist()
print(f"\nColunas numéricas: {colunas_numericas}")

print("\n=== VISUALIZAÇÕES ===")

# 1. Histogramas por classe
if len(lista_classes) > 1 and len(colunas_numericas) > 0:
    print("\nCriando histogramas...")
    
    n_cols = min(3, len(colunas_numericas))
    n_rows = (len(colunas_numericas) + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    cores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, coluna in enumerate(colunas_numericas):
        row = i // n_cols
        col = i % n_cols
        
        if n_rows == 1:
            ax = axes[col] if n_cols > 1 else axes[0]
        else:
            ax = axes[row, col] if n_cols > 1 else axes[row]
        
        for j, classe in enumerate(lista_classes):
            dados_classe = df[df['Classe'] == classe][coluna]
            dados_limpos = dados_classe.dropna()
            if len(dados_limpos) > 0:
                ax.hist(dados_limpos, bins=10, alpha=0.6, 
                       label=f'{classe} (n={len(dados_limpos)})', 
                       color=cores[j % len(cores)], density=True)
        
        ax.set_xlabel(coluna)
        ax.set_ylabel('Densidade')
        ax.set_title(f'Distribuição de {coluna}')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Remover subplots extras
    for i in range(len(colunas_numericas), n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        if n_rows == 1:
            fig.delaxes(axes[col] if n_cols > 1 else axes[0])
        else:
            fig.delaxes(axes[row, col] if n_cols > 1 else axes[row])
    
    plt.tight_layout()
    plt.show()

# 2. Boxplots
print("\nCriando boxplots...")
if len(colunas_numericas) > 0:
    fig, axes = plt.subplots(1, len(colunas_numericas), 
                            figsize=(4*len(colunas_numericas), 6))
    
    if len(colunas_numericas) == 1:
        axes = [axes]
    
    for i, coluna in enumerate(colunas_numericas):
        dados_coluna = df[coluna].dropna()
        if len(dados_coluna) > 0:
            sns.boxplot(data=df, y=coluna, x='Classe', ax=axes[i])
            axes[i].set_title(f'Boxplot - {coluna}')
            axes[i].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()

# 3. Matriz de correlação
print("\nCriando matriz de correlação...")
if len(colunas_numericas) > 1:
    dados_num = df[colunas_numericas]
    correlacao = dados_num.corr()
    
    plt.figure(figsize=(10, 8))
    mask = np.triu(np.ones_like(correlacao, dtype=bool))
    sns.heatmap(correlacao, annot=True, cmap='coolwarm', center=0,
                square=True, fmt='.3f', mask=mask)
    plt.title('Matriz de Correlação')
    plt.tight_layout()
    plt.show()
    
    print("\nCorrelações mais fortes (>0.5):")
    for i in range(len(correlacao.columns)):
        for j in range(i+1, len(correlacao.columns)):
            corr_val = correlacao.iloc[i, j]
            if abs(corr_val) > 0.5:
                print(f"{correlacao.columns[i]} vs {correlacao.columns[j]}: {corr_val:.3f}")

# 4. Scatter plot
if len(colunas_numericas) >= 2:
    print("\nCriando scatter plot...")
    var_x = colunas_numericas[0]
    var_y = colunas_numericas[1]
    
    plt.figure(figsize=(10, 6))
    cores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for j, classe in enumerate(lista_classes):
        dados_classe = df[df['Classe'] == classe]
        dados_x = dados_classe[var_x].dropna()
        dados_y = dados_classe[var_y].dropna()
        
        # Pegar apenas dados onde ambas as variáveis não são NaN
        indices_validos = dados_classe[[var_x, var_y]].dropna().index
        x_vals = dados_classe.loc[indices_validos, var_x]
        y_vals = dados_classe.loc[indices_validos, var_y]
        
        if len(x_vals) > 0:
            plt.scatter(x_vals, y_vals, 
                       label=f'{classe} (n={len(x_vals)})', 
                       alpha=0.7, color=cores[j % len(cores)], s=50)
    
    plt.xlabel(var_x)
    plt.ylabel(var_y)
    plt.title(f'Scatter Plot: {var_x} vs {var_y}')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

print("\n=== ANÁLISE CONCLUÍDA ===")
print(f"Dataset final: {len(df)} proteínas")
print(f"Classes: {list(lista_classes)}")
print(f"Variáveis numéricas: {len(colunas_numericas)}")

# Salvar dados
df.to_csv('proteinas_uniprot_processadas.csv', index=False)
print("\nDados salvos em 'proteinas_uniprot_processadas.csv'")
print("Análise concluída!")
