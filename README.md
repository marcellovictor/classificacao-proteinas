# Base de Dados UniProt

## Obtenção dos Dados

A base de dados do UniProt utilizada neste projeto não está incluída no repositório devido ao seu tamanho, que ultrapassa os limites permitidos pelo sistema de controle de versão.

### Como obter a base de dados

Para obter a base de dados necessária, siga os passos abaixo:

1. Acesse o link: https://www.uniprot.org/uniprotkb?query=*&facets=reviewed%3Atrue

2. Faça o download do dataset **Reviewed (Swiss-Prot)** em formato `.tsv`

3. Durante o processo de download, certifique-se de selecionar as seguintes colunas:
   - **Sequence**
   - **Gene Ontology IDs**
   - **Length**
   - **Mass**

### Configuração

Após o download, coloque o arquivo `.tsv` no diretório do projeto.
