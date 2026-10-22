# Flamingo: metamaterial acústico com ressoadores de Helmholtz

Análise das medições de um metamaterial acústico: um tubo de PVC com cinco seringas, que funcionam como ressoadores de Helmholtz, e uma terminação de lã ou de LEGO. Estudamos a ressonância local de cada seringa, o efeito de Bragg da rede periódica e a fusão das duas bandas.

**Autoras:** Ana Sousa (ist1118764), Matilde Salgado (ist1118517), Matilde Falcão (ist1118714), Yasmin Altoe (ist1117577)

**Disciplina:** Laboratório de Física Experimental (LFE), projeto de Acústica, IST, 2026

## Pastas

| Pasta | Experiência |
|---|---|
| `teto_sistema/` | Teto do sistema: maior atenuação mensurável, limitada pelo ruído e pela fuga de som por fora do tubo |
| `correcao_gargalo/` | Calibração do gargalo: f0 em função do volume da seringa, comprimento efetivo L_ef e correção ΔL |
| `numero_elementos/` | Atenuação (TL) e largura da banda em função do número de seringas ativas (N = 1 a 5) |
| `distancia_altifalante/` | Ressonância de uma seringa em função da distância ao altifalante |
| `bragg_vs_local/` | Ressonância local vs. Bragg: fusão a 4,1 mL e separação a 5 mL |
| `reflexao_2_microfones/` | Reflexão de uma seringa e ondas estacionárias com a terminação de lã e com a de LEGO |
| `decaimento/` | Decaimento da ressonância depois de cortar o tom (Q pelo tempo) |
| `incertezas/` | Balanço de incertezas (Apêndice C do relatório) |

`estilo_tl.py` tem o estilo comum dos gráficos de TL, o modelo de matrizes de transferência e o cálculo do teto, que vários scripts usam.

Cada pasta tem um `README.md` com o resumo da experiência, as gravações e o passo a passo do código.

## Dependências

Python 3.10 ou mais recente.

```
pip install -r requirements.txt
```

## Como correr

Dentro da pasta da experiência:

```
cd numero_elementos
python tl_vs_N.py
```

Cada script escreve os resultados no terminal e guarda as figuras `.png` na mesma pasta. As figuras não estão no repositório: estão no relatório e são geradas pelos scripts.
