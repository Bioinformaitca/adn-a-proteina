# Del ADN a la proteína: replicación, transcripción y traducción con Biopython

**Asignatura:** Bioinformática · **Autor:** Aimar 

**Repositorio:** <https://github.com/Bioinformaitca/adn-a-proteina>

## Resumen

Se resuelven seis ejercicios sobre el dogma central de la biología molecular. Cada uno se hace primero a mano y después se comprueba con un script en Python con Biopython (`src/`). El trabajo termina con un *pipeline* (`ej6_pipeline.py`) que replica, transcribe y traduce una secuencia real: la CDS del gen humano *HBB* (β-globina, NCBI RefSeq NM_000518.5). El programa va informando de cada paso.

**Ejecución:** `pip install -r requirements.txt` y, por ejemplo, `python src/ej6_pipeline.py --demo-mutacion`. La salida completa de todos los scripts está en `resultados/salida_consola.txt`.

## Ejercicio 1. Replicación del ADN

**Manual.** La helicasa separa las hebras y cada una sirve de molde (replicación semiconservativa). La nueva hebra se sintetiza siempre 5'→3' y es antiparalela a su molde:

| Molécula hija | Hebra parental | Hebra nueva |
|---|---|---|
| 1 | 5'-ATG CCG TTA GCT-3' | 3'-TAC GGC AAT CGA-5' |
| 2 | 3'-TAC GGC AAT CGA-5' | 5'-ATG CCG TTA GCT-3' |

**Enzimas.** La *helicasa* rompe los puentes de hidrógeno y abre la horquilla. La *primasa* sintetiza un cebador corto de ARN con un 3'-OH libre, porque la polimerasa no puede empezar de cero. La *ADN polimerasa* añade nucleótidos complementarios en sentido 5'→3': de forma continua en la hebra adelantada y en fragmentos de Okazaki en la retrasada. Además, con su actividad exonucleasa 3'→5' corrige errores y sustituye los cebadores por ADN. La *ligasa* sella los cortes entre fragmentos mediante enlaces fosfodiéster.

**Reflexión.** Si la polimerasa se equivoca y el error no se corrige (ni por la lectura de prueba ni por la reparación de desapareamientos), en la siguiente ronda la base incorrecta sirve de molde y la mutación queda fijada en una de las dos células hijas. Su efecto puede ser nulo (mutación silenciosa), cambiar un aminoácido (de sentido erróneo), crear un codón de paro prematuro (sin sentido) o, si hay inserciones o deleciones, desplazar el marco de lectura.

**Biopython.** `ej1_replicacion.py` usa `Seq.complement()` y obtiene `TACGGCAATCGA`, igual que el resultado manual (`True`).

## Ejercicio 2. Transcripción

**Cadena molde.** La hebra 5'-ATG CCT GAA TGC-3' empieza por el codón de inicio ATG, así que es la **codificante** (sentido). La **molde** es la 3'-TAC GGA CTT ACG-5'. La ARN polimerasa lee el molde 3'→5' y sintetiza el ARN 5'→3', con U en lugar de T:

**ARNm: 5'-AUG CCU GAA UGC-3'**, que tiene la misma secuencia que la hebra codificante pero con U.

**Promotor y región codificante.** El promotor no aparece en el fragmento. Estaría *aguas arriba* (hacia 5') del inicio de la transcripción; en eucariotas incluye la caja TATA (≈ −25 a −30), a la que se unen los factores de transcripción y la ARN polimerasa II. La región codificante empieza en el ATG (+1 del marco de lectura) y abarca los 12 nt dados. Como no hay codón de paro, el fragmento sigue más allá.

**Biopython.** `ej2_transcripcion.py` lee un FASTA y obtiene el mismo ARNm por dos vías: el complemento del molde y `Seq.transcribe()`. En el experimento de orientación se toma la otra hebra como codificante (el complemento inverso): sale 5'-GCA UUC AGG CAU-3', que no empieza por AUG y daría un péptido distinto (AFRH en vez de MPEC). Por eso la orientación de la hebra es crítica.

## Ejercicio 3. Traducción

5'-**AUG** UAU GCU **UAA**-3': AUG es el codón de inicio (Met) y UAA el de paro (ocre). Resultado: **Met-Tyr-Ala** (MYA). `Bio.Seq.translate(to_stop=True)` da lo mismo.

**Mutación AUG→GUG.** En eucariotas el ribosoma no reconoce el inicio y sigue explorando el ARNm hasta el siguiente AUG. Aquí aparece uno en otro marco (U**AUG**CU…) y el script obtiene *Met-Leu* sin codón de paro, es decir, una proteína distinta o que no se produce. En bacterias GUG puede actuar como inicio alternativo y se incorporaría igualmente fMet, aunque con menos eficiencia.

**Pérdida del codón de paro (UAA→CAA).** El ribosoma no se detiene y traduce la región 3' no traducida (*read-through*). Se forma una proteína alargada (en la simulación, Met-Tyr-Ala-Gln-Gly-Phe…), que suele ser inestable. Además, un ARNm sin codón de paro activa la vía de degradación *non-stop decay*.

## Ejercicio 4. Splicing alternativo

Se diseñó un gen modelo de cinco exones: E1 contiene el AUG, E5 el codón de paro y E3 mide 40 nt (no es múltiplo de 3). Se compararon tres isoformas (`ej4_splicing.py`):

| Isoforma | nt | Marco | Proteína |
|---|---|---|---|
| 1-2-3-4-5 (canónica) | 147 | conservado | 48 aa |
| 1-2-4-5 (salto de E3) | 107 | desplazado | 28 aa, extremo C distinto |
| 1-3-5 (salto de E2 y E4) | 88 | desplazado | 29 aa, sin dominio A |

**Diferencias esperadas.** Si el exón omitido es múltiplo de 3, la proteína solo pierde el dominio que codificaba (por ejemplo, un sitio de unión o una región transmembrana). Si no lo es, cambia el marco de lectura: la secuencia es distinta desde el punto de unión y suele aparecer un codón de paro prematuro, que puede provocar la degradación del ARNm por NMD. Las isoformas pueden variar en su localización, su afinidad por ligandos o su regulación.

**Diversidad sin más genes.** Con *n* exones opcionales se pueden obtener hasta 2ⁿ combinaciones a partir de un único locus. Así, unos 20.000 genes humanos dan lugar a más de 100.000 proteínas distintas, y cada tejido puede expresar las isoformas que necesita.

**FGFR2 en Ensembl** (ENSG00000066468, cromosoma 10, hebra −; consultado en septiembre de 2026). Ensembl lista 25 transcritos: isoformas codificantes de 371 a 822 aa, transcritos con intrón retenido y transcritos degradados por NMD (`data/FGFR2_transcritos_ensembl.tsv`). El caso clásico son los exones **IIIb y IIIc**, mutuamente excluyentes, que forman la mitad del dominio Ig-III de unión al ligando. La isoforma IIIb (FGFR2-215, 822 aa) se expresa en epitelios y une FGF7/FGF10. La IIIc (FGFR2-206, canónica, 821 aa) se expresa en el mesénquima y une FGF2. Un solo exón cambia qué factores de crecimiento reconoce el receptor, y un cambio de isoforma se ha relacionado con la transición epitelio-mesénquima en cáncer. Las isoformas cortas (371-593 aa), que no tienen parte del dominio extracelular o del dominio quinasa, tendrían una señalización reducida o nula.

## Ejercicio 5. Proteínas

**Extremos.** H₂N-**Met**-Ile-Ser-Gly-Val-Lys-**His**-COOH. El extremo N-terminal es Met (grupo amino libre, primer residuo que se traduce) y el C-terminal es His (grupo carboxilo libre). El péptido mezcla residuos apolares (M, I, V, G) y polares o básicos (S, K, H). Valores de `ProtParam`: GRAVY +0,33, pI ≈ 8,5 y masa 770,9 Da.

**El orden importa.** La secuencia (estructura primaria) determina qué puentes de hidrógeno, interacciones hidrofóbicas, puentes salinos y ángulos φ/ψ son posibles. Por eso fija la estructura secundaria y el plegamiento terciario (hipótesis de Anfinsen). Con los mismos aminoácidos en otro orden se obtiene otra proteína.

**Cambio hidrofóbico→hidrofílico en el interior.** El núcleo de una proteína globular es hidrofóbico; el efecto hidrofóbico es la principal fuerza que impulsa el plegamiento. Colocar ahí un residuo polar o cargado supone un gran coste energético: no puede formar puentes de hidrógeno con el agua y rompe el empaquetamiento. Las consecuencias son pérdida de estabilidad, plegamiento incorrecto, agregación o degradación, y por tanto pérdida de función.

**PDB: ubiquitina (1UBQ, rayos X, 1,8 Å).** Según los registros HELIX/SHEET, tiene una hélice α principal (residuos 23-34), una hélice 3₁₀ corta (56-59) y una lámina β mixta de cinco hebras (1-7, 10-17, 40-45, 48-50, 64-72): 21 % hélice y 43 % lámina. La hélice se apoya sobre la cara cóncava de la lámina, y entre ambas queda un núcleo hidrofóbico con residuos como Ile3, Val5, Ile23, Val26 y Leu67. En la simulación de la mutación **V26K** (valina de la hélice → lisina), la hidropatía del residuo pasa de +4,2 a −3,9. El GRAVY global apenas cambia (de −0,49 a −0,60), pero introducir una carga en el núcleo desestabilizaría el plegamiento. Esto muestra que un solo cambio puntual puede ser decisivo aunque la proteína «promedio» casi no varíe.

## Ejercicio 6. Pipeline integrador (ADN → proteína)

**Secuencia.** CDS de *HBB* (NM_000518.5, 444 nt, 56,1 % GC), descargada de NCBI en formato FASTA (`data/HBB_cds.fasta`).

**Pipeline (`ej6_pipeline.py`).** Registra cada paso con `logging`:

1. **Validación:** solo admite A, C, G y T, y calcula la longitud y el %GC.
2. **Replicación:** genera las dos hebras nuevas (complemento inverso de cada parental) y comprueba que las dos moléculas hijas son idénticas a la original.
3. **Transcripción:** toma como molde la hebra antisentido y obtiene el ARNm (`reverse_complement_rna`). Comprueba que coincide con la codificante cambiando T por U.
4. **Traducción:** busca el primer AUG, traduce hasta el codón de paro e informa de su posición, del número de codones y del codón de paro.
5. **Salida:** guarda un FASTA con las hebras nuevas, el ARNm y la proteína en `resultados/`.

**Resultado.** Inicio en el nucleótido 1, 148 codones y parada en UAA. Se obtiene una proteína de **147 aa** (MVHLTPEEKSAVTALWGKVNVDEVGG…KYH), que coincide con NP_000509.1. En la proteína madura se elimina la Met inicial y quedan 146 aa.

**¿Qué paso es más vulnerable?** La **replicación**. Los errores de transcripción y traducción afectan a una sola molécula de ARN o de proteína, que se degrada y se sustituye; son transitorios. Un error de replicación no corregido, en cambio, pasa a todas las células descendientes y a todos los ARNm y proteínas que se produzcan a partir de él. La opción `--demo-mutacion` lo ilustra: el cambio A20T en el ADN convierte el codón GAG en GTG y produce **E7V** en la cadena naciente (Glu6Val en la proteína madura). Es la mutación de la **anemia falciforme**: un solo nucleótido cambia la carga de la superficie de la hemoglobina y provoca que polimerice. Dentro de la traducción, lo más delicado es la selección del marco: un error en el inicio o una inserción o deleción cambian toda la proteína que viene después.

## Conclusiones

El resultado manual y el computacional coinciden en todos los ejercicios. Para un análisis automatizado correcto hay que dejar explícitas la orientación de las hebras (5'→3'), la elección de la cadena molde y el marco de lectura. Biopython (`Seq`, `SeqIO`, `CodonTable`, `ProtParam`) permite hacer estas comprobaciones de forma reproducible.

**Estructura del repositorio:** `src/` (scripts ej1-ej6) · `data/` (FASTA de HBB, tabla de FGFR2 de Ensembl, estructura secundaria de 1UBQ) · `resultados/` (FASTA generados y salida de consola).

**Fuentes:** NCBI RefSeq NM_000518.5 · Ensembl ENSG00000066468 (FGFR2) · RCSB PDB 1UBQ (Vijay-Kumar *et al.*, 1987, *J. Mol. Biol.* 194:531) · Cock *et al.* (2009) Biopython, *Bioinformatics* 25:1422 · Kyte & Doolittle (1982), *J. Mol. Biol.* 157:105.
