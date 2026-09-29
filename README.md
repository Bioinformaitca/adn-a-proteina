# Del ADN a la proteína: replicación, transcripción y traducción con Biopython

**Asignatura:** Bioinformática

**Autores:** Aimar Alejandro Santana y Daniel Perdomo Medina

**Repositorio:** <https://github.com/Bioinformaitca/adn-a-proteina>

## Resumen

Resolvemos seis ejercicios sobre el dogma central de la biología molecular. Cada uno lo realizamos primero a mano y posteriormente lo comprobamos mediante un script en Python con Biopython (`src/`). Nuestro trabajo concluye con un pipeline (`ej6_pipeline.py`) que replica, transcribe y traduce la CDS del gen humano HBB (β-globina, NCBI RefSeq NM_000518.5), informando progresivamente de cada paso.

**Ejecución:**

```bash
pip install -r requirements.txt
python src/ej6_pipeline.py --demo-mutacion
```

*Salida completa:* `resultados/salida_consola.txt`

## Ejercicio 1. Replicación del ADN

**Manual:** La helicasa separa las hebras en una replicación semiconservativa. La nueva hebra se sintetiza siempre en sentido 5’→3’ y resulta antiparalela a su hebra molde:

| Molécula hija | Hebra parental | Hebra nueva | 
| :--- | :--- | :--- | 
| 1 | `5’-ATG CCG TTA GCT-3’` | `3’-TAC GGC AAT CGA-5’` | 
| 2 | `3’-TAC GGC AAT CGA-5’` | `5’-ATG CCG TTA GCT-3’` | 

**Enzimas involucradas:**

* **Helicasa:** Rompe los puentes de hidrógeno y abre la horquilla de replicación.
* **Primasa:** Sintetiza un cebador corto de ARN con un grupo 3’-OH libre necesario para que la polimerasa pueda iniciar la síntesis.
* **ADN polimerasa:** Añade nucleótidos complementarios en sentido 5’→3’ (de forma continua en la hebra adelantada y mediante fragmentos de Okazaki en la retrasada). Su actividad exonucleasa 3’→5’ permite corregir errores y reemplazar los cebadores por ADN.
* **Ligasa:** Sella los cortes entre fragmentos formando enlaces fosfodiéster.

**Reflexión:** Si la polimerasa comete un error no corregido por la lectura de prueba o los sistemas de reparación, en la siguiente ronda la base incorrecta servirá de molde, fijando la mutación en una de las células hijas. Sus efectos pueden ser nulos (mutación silenciosa), cambiar un aminoácido (sentido erróneo), generar un codón de paro prematuro (sin sentido) o alterar el marco de lectura en caso de inserciones o deleciones.

**Biopython:** Nuestro script `ej1_replicacion.py` emplea `Seq.complement()` y genera `TACGGCAATCGA`, coincidiendo con nuestro resultado manual (True).

## Ejercicio 2. Transcripción

**Cadena molde:** La secuencia `5’-ATG CCT GAA TGC-3’` inicia con el codón ATG, por lo que es la hebra codificante (sentido). La hebra molde corresponde a `3’-TAC GGA CTT ACG-5’`. La ARN polimerasa lee la hebra molde en sentido 3’→5’ y sintetiza el ARNm en sentido 5’→3’, reemplazando T por U:

* **ARNm:** `5’-AUG CCU GAA UGC-3’` (idéntica a la hebra codificante pero con U).

**Promotor y región codificante:** El promotor se ubica aguas arriba (hacia 5’) del inicio de transcripción (en eucariotas incluye la caja TATA a ≈ -25 o -30, donde se unen los factores de transcripción y la ARN polimerasa II). La región codificante inicia en el ATG (+1 del marco de lectura) y abarca los 12 nucleótidos provistos.

**Biopython:** `ej2_transcripcion.py` lee un archivo FASTA y obtiene el mismo ARNm mediante el complemento del molde y mediante `Seq.transcribe()`. Al tomar la hebra opuesta como codificante (complemento inverso), obtenemos `5’-GCA UUC AGG CAU-3’`, que al no iniciar en AUG produciría el péptido AFRH en lugar de MPEC, evidenciando la importancia de la orientación.

## Ejercicio 3. Traducción

**Secuencia:** `5’-AUG UAU GCU UAA-3’`. AUG actúa como codón de inicio (Met) y UAA como codón de paro (ocre).

* **Resultado:** Met-Tyr-Ala (MYA), confirmado por nosotros con `Bio.Seq.translate(to_stop=True)`.

**Mutación AUG→GUG:** En eucariotas, el ribosoma omite este inicio y escanea el ARNm hasta el siguiente AUG en otro marco (UAUGCU…), obteniendo Met-Leu sin codón de paro (proteína distinta o no producida). En bacterias, GUG puede operar como inicio alternativo e incorporar fMet con menor eficiencia.

**Pérdida del codón de paro (UAA→CAA):** Causa read-through, traduciendo la región 3’ no traducida para formar una proteína alargada e inestable (Met-Tyr-Ala-Gln-Gly-Phe…), activando la vía de degradación non-stop decay.

## Ejercicio 4. Splicing alternativo

A partir de un gen modelo con 5 exones (E1 con AUG, E5 con codón de paro y E3 de 40 nt, no múltiplo de 3), en nuestro script `ej4_splicing.py` analizamos tres isoformas:

| Isoforma | Longitud (nt) | Marco de lectura | Proteína | 
| :--- | :--- | :--- | :--- | 
| 1-2-3-4-5 (canónica) | 147 nt | Conservado | 48 aa | 
| 1-2-4-5 (salto de E3) | 107 nt | Desplazado | 28 aa (extremo C distinto) | 
| 1-3-5 (salto de E2 y E4) | 88 nt | Desplazado | 29 aa (sin dominio A) | 

**Diferencias esperadas:** Si el exón omitido es múltiplo de 3, la proteína únicamente pierde el dominio codificado. Si no lo es, se altera el marco de lectura, generando una secuencia distinta a partir de la unión y provocando frecuentemente codones de paro prematuros que activan la degradación por NMD.

**Diversidad proteica:** $n$ exones opcionales pueden generar hasta $2^n$ combinaciones desde un solo locus, permitiendo que ≈ 20.000 genes humanos produzcan más de 100.000 proteínas distintas.

**FGFR2 en Ensembl (ENSG00000066468, cromosoma 10, hebra −):** Ensembl registra 25 transcritos (`data/FGFR2_transcritos_ensembl.tsv`).

* **Isoforma IIIb (FGFR2-215, 822 aa):** Incluye el exón IIIb, se expresa en epitelios y une FGF7/FGF10.
* **Isoforma IIIc (FGFR2-206, canónica, 821 aa):** Incluye el exón IIIc, se expresa en mesénquima y une FGF2.

El cambio entre estas isoformas se relaciona con la transición epitelio-mesénquima en cáncer. Isoformas cortas (371-593 aa) sin dominios extracelulares o quinasa presentan señalización reducida o nula.

## Ejercicio 5. Proteínas

**Extremos:** H₂N-Met-Ile-Ser-Gly-Val-Lys-His-COOH

* **Extremo N-terminal:** Met (grupo amino libre).
* **Extremo C-terminal:** His (grupo carboxilo libre).

**Valores ProtParam:** GRAVY +0,33, pI ≈ 8,5, masa 770,9 Da.

**El orden de los residuos:** La estructura primaria determina los puentes de hidrógeno, interacciones hidrofóbicas, puentes salinos y ángulos φ/ψ, guiando la estructura secundaria y el plegamiento terciario (hipótesis de Anfinsen).

**Núcleo hidrofóbico:** Introducir un residuo polar o cargado en el interior hidrofóbico supone un elevado coste energético al no poder formar puentes de hidrógeno con el agua, provocando desestabilización, agregación o pérdida de función.

**PDB Ubiquitina (1UBQ, 1.8 Å):** Posee una hélice α principal (residuos 23-34), una hélice 3₁₀ corta (56-59) y una lámina β mixta de 5 hebras (1-7, 10-17, 40-45, 48-50, 64-72), totalizando 21 % de hélice y 43 % de lámina.

**Mutación V26K:** Cambia una valina del núcleo por lisina. La hidropatía del residuo pasa de +4,2 a -3,9. Aunque el GRAVY global apenas cambia (de -0,49 a -0,60), la introducción de una carga en el núcleo desestabiliza el plegamiento.

## Ejercicio 6. Pipeline integrador (ADN → proteína)

**Secuencia objetivo:** CDS de HBB (NM_000518.5, 444 nt, 56,1 % GC), almacenada en `data/HBB_cds.fasta`.

**Pasos de nuestro `ej6_pipeline.py`:**

1. **Validación:** Comprueba nucleótidos A, C, G, T, longitud y %GC.
2. **Replicación:** Genera las dos hebras hijas y valida la exactitud con la original.
3. **Transcripción:** Utiliza la hebra antisentido como molde para generar el ARNm.
4. **Traducción:** Detecta el primer AUG, traduce la secuencia hasta el codón de paro e informa de los parámetros.
5. **Salida:** Exporta los archivos FASTA resultantes en `resultados/`.

**Resultados:** Inicio en el nt 1, 148 codones y codón de paro UAA. Proteína obtenida de 147 aa (MVHLTPEEKSAVTALWGKVNVDEVGG…KYH), coincidente con NP_000509.1 (146 aa tras la eliminación de la Met inicial).

**Vulnerabilidad de procesos:** La replicación es el paso más vulnerable. Errores en transcripción o traducción son transitorios y afectan únicamente a moléculas individuales de ARN o proteína que se degradan. Un error de replicación no corregido se transmite a todas las células descendientes.

**Anemia falciforme (`--demo-mutacion`):** La mutación A20T convierte el codón GAG en GTG, produciendo el cambio E7V (Glu6Val en la proteína madura), alterando la carga superficial de la hemoglobina y provocando su polimerización.

## Conclusiones

Nuestros análisis manuales y computacionales coinciden plenamente. Para nuestras automatizaciones bioinformáticas precisas es indispensable especificar la orientación 5’→3’, la cadena molde utilizada y el marco de lectura.

## Estructura del repositorio

* `src/`: Scripts de los ejercicios 1 a 6.
* `data/`: Archivos FASTA de HBB, tablas de Ensembl para FGFR2 y estructura secundaria de 1UBQ.
* `resultados/`: Salidas FASTA y registros de consola.

## Fuentes bibliográficas

* NCBI RefSeq NM_000518.5
* Ensembl ENSG00000066468 (FGFR2)
* RCSB PDB 1UBQ (Vijay-Kumar et al., 1987)
* Cock et al. (2009) Biopython
* Kyte & Doolittle (1982)