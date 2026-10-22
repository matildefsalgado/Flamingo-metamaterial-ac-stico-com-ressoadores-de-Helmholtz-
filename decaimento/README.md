# Decaimento da ressonância

**O que mede:** o tempo de decaimento τ da seringa 1 (14,6 mL) depois de se cortar um tom a 430 Hz, e daí Q = π·f0·τ. Compara com um controlo, com todos os êmbolos no fundo.

**Script:** `decaimento.py`. Lê as gravações e gera `decaimento.png`, e escreve τ da seringa e do controlo, a razão entre os dois e Q.
- `--teste` corre o script com sinais sintéticos de τ conhecido.
- `--f0`, `--banda` e `--pasta` mudam a frequência, o filtro e a pasta das gravações.

**Gravações (22/10):** `decaimento 1.wav` (seringa) e `controlo (1).wav` (controlo).

**Parâmetros:** filtro de f0 ± 150 Hz, envolvente de Hilbert suavizada a 1 ms e ajuste de −3 a −25 dB. Resultado: τ da seringa ≈ τ do controlo ≈ 19 ms, por isso o decaimento é do tubo e só se pode dizer que Q < 26.

## Passo a passo do código

```text
DECAIMENTO DA RESSONÂNCIA (Q pelo tempo, seringa 1 a 14,6 mL)
==============================================================

Uso:  python decaimento.py
      python decaimento.py --f0 435              (se o tom foi posto noutra frequência)
      python decaimento.py --banda 50            (filtro mais estreito, se houver outras frequências perto)
      python decaimento.py --pasta <outra pasta> (se as gravações estiverem noutro sítio)
      python decaimento.py --teste               (sinais sintéticos com τ conhecido, para ver se o
                                                  código funciona; figura decaimento_teste.png)

MONTAGEM
  A de sempre: lã na ponta, microfone na posição 2 (jusante), Audacity a 48 kHz, mono.
  Gerador em SENO contínuo a f0 = 430 Hz (medido: 430 Hz na aula 3, 435 Hz a 14/10).
  Não precisa de ser exato: depois do corte a seringa vibra à frequência dela.
  NÃO mexer no gerador, no volume nem no ganho entre as duas gravações.

AS GRAVAÇÕES (nesta pasta, mesma sessão)
  dec_seringa_1.wav    seringa 1 a 14,6 mL, as outras com o êmbolo no fundo
  dec_controlo_1.wav   todos os êmbolos no fundo (mesmo tom, mesma frequência)
  Em cada gravação: tom ~2 s ligado, ~1 s desligado, 5 vezes seguidas (~15 s).
  O corte tem de ser SECO (botão de saída do gerador, ou modo burst se houver). Rodar o volume
  para baixo não serve: o decaimento dura só ~10 ms.
  Mais repetições: dec_seringa_2.wav, ... O código junta todos os cortes de todas as gravações.
  Também aceita "-" em vez de "_" e o " (1)" que o Windows acrescenta.

PORQUÊ O CONTROLO
  O sistema (altifalante, tubo, lã) também demora a calar-se: uma ida e volta em 1,75 m são ~10 ms.
  O decaimento medido nunca é mais rápido do que o do sistema. Se τ_seringa não for claramente
  maior que τ_controlo (pelo menos o dobro), o Q que sai é só um limite superior.

O QUE O CÓDIGO FAZ
  1. Filtra cada gravação à volta de f0 (f0 ± 150 Hz, fase zero), tira a envolvente (Hilbert) e
     suaviza-a 1 ms. Em dB. O filtro tira o ruído e o estalo do altifalante no corte.
  2. Encontra os cortes: envolvente suavizada a 50 ms, limiar 20 dB abaixo do tom ligado
     (os silêncios do Audacity podem ser zeros exatos: a envolvente tem chão a −120 dB); só conta se o tom esteve ligado ≥ 0,5 s antes e desligado ≥ 0,3 s depois.
  3. Em cada corte: nível ligado (mediana de −250 a −50 ms), ruído (mediana de +150 a +300 ms).
     Ajusta uma reta à queda em dB, de 3 dB abaixo do nível ligado até 25 dB abaixo (ou até 10 dB
     acima do ruído, o que vier primeiro).
       τ = −20 / (declive · ln 10)      Q = π · f0 · τ      T60 = 6,9 · τ
  4. Média dos τ de todos os cortes, incerteza = desvio padrão da média. O mesmo no controlo.
     Razão τ_seringa/τ_controlo e z da diferença; aviso se a razão for < 2.
  5. Incerteza de Q: τ e f0 (± 5 Hz) em quadratura.
  6. Compara com o Q pela largura do pico: 13,2 (N1, aula 3) e 15,1 ± 1,1 (seringa 1, 14/10).
  7. Limite do método: o filtro e a suavização também "tocam". O código mede o τ que vê num tom
     cortado instantaneamente (~1,1 ms) e avisa se algum τ medido estiver a menos de 3× disso.
  8. Figura decaimento.png: à esquerda as quedas de todos os cortes (seringa e controlo) com os
     ajustes; à direita a primeira gravação inteira com os cortes encontrados.

TESTE (--teste): seringa τ = 10 ms, controlo τ = 3 ms, ruído a −50 dB.
  Recupera 10,00 ms e 2,89 ms, Q = 13,5 (verdadeiro π·430·0,010 = 13,5).

O QUE ESPERAR
  Com Q = 13,2 e f0 = 430 Hz: τ = Q/(π f0) ≈ 10 ms, T60 ≈ 70 ms.
  O decaimento é uma medição independente do Q: não muda as conclusões já tiradas, só confirma
  (ou não) o Q da largura do pico.

GRAVAÇÕES DE 22/10 (nomes usados na aula; o código aceita-os)
  decaimento 1.wav    seringa 1 a 14,6 mL, tom a 430 Hz cortado 5 vezes   -> entra como seringa
  controlo (1).wav    êmbolos todos no fundo, mesmo tom, 5 cortes           -> entra como controlo
  decaimento.wav      primeira tentativa, tom a 600 Hz (frequência errada)  -> não entra
  A_referencia.wav    varrimento normal com os êmbolos no fundo             -> não entra (não há
                      varrimento com a seringa nesta sessão para fazer um TL)

RESULTADO (22/10)
  τ seringa = 19,2 ± 1,7 ms, τ controlo = 19,1 ± 0,5 ms (razão 1,00, z = 0,0).
  O decaimento é o do próprio tubo, não o da seringa: Q = π f0 τ = 25,9 ± 2,3 é só um limite superior
  (Q < 26), compatível com o Q = 13–15 da largura do pico, mas não é uma segunda medida independente.
```
