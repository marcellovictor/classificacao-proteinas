import pandas as pd
from Bio.SeqUtils import ProtParam
from Bio.Seq import Seq
import numpy as np

# Carregar arquivo .tsv exportado do UniProt
caminho_arquivo = "uniprotkb.tsv"  # Altere para o nome do seu arquivo
df = pd.read_csv(caminho_arquivo, sep="\t")


# Função para calcular propriedades da sequência
def calcular_propriedades(seq_str):
    try:
        seq = Seq(seq_str)
        analyser = ProtParam.ProteinAnalysis(str(seq))

        pi = analyser.isoelectric_point()
        gravy = analyser.gravy()
        charge = analyser.charge_at_pH(7.0)
        aa_count = analyser.count_amino_acids()

        total = sum(aa_count.values())
        polar = sum([aa_count.get(aa, 0) for aa in ['Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'W']])
        apolar = sum([aa_count.get(aa, 0) for aa in ['A', 'V', 'L', 'I', 'P', 'F', 'M', 'G']])

        proporcao_polar = polar / total if total > 0 else 0
        proporcao_apolar = apolar / total if total > 0 else 0

        return pd.Series([pi, gravy, charge, proporcao_polar, proporcao_apolar])

    except Exception as e:
        print(f"Erro na sequência: {e}")
        return pd.Series([np.nan]*5)

# Aplicar função
df[['Ponto_Isoeletrico', 'Hidrofobicidade', 'Carga_Total', 
    'Proporcao_Polar', 'Proporcao_Apolar']] = df['Sequence'].apply(calcular_propriedades)

print(df.head())
print("\n\n\n")
print(df.describe)
print("\n\n\n")
print(df.info())
print("\n\n\n")
print(df.columns)
print("\n\n\n")

# Exibir colunas selecionadas

df_selecionado = df[['Sequence', 'Mass', 'Ponto_Isoeletrico', 'Hidrofobicidade', 'Carga_Total', 'Proporcao_Polar', 'Proporcao_Apolar', 'Length', 'Gene Ontology IDs']]

print(df_selecionado.head())
print("\n\n\n")

# Salvar como novo CSV
df.to_csv("uniprot_com_features.csv", index=False)
print("Arquivo salvo como uniprot_com_features.csv")
