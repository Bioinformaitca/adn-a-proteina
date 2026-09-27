"""Ejercicio 4 - Splicing alternativo con un gen modelo de 5 exones.

Se define un gen de juguete con exones de longitudes distintas (en nt) y se
comparan varias isoformas: longitud, conservación del marco de lectura y
dominios presentes. Después se resumen los transcritos reales de FGFR2
descargados de Ensembl (data/FGFR2_transcritos_ensembl.tsv).
"""
import csv
from pathlib import Path
from Bio.Seq import Seq

RAIZ = Path(__file__).resolve().parent.parent

# Gen modelo: E1 lleva el AUG, E5 el codón de paro. Longitudes elegidas para
# que omitir E3 (40 nt, no múltiplo de 3) desplace el marco de lectura.
EXONES = {
    1: "ATGGCTAGCAAAGGA",                                   # 15 nt  péptido señal
    2: "GAACTGTTCACCGGGGTGGTGCCCATCCTG",                    # 30 nt  dominio A
    3: "GTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCG",          # 40 nt  dominio B (no múltiplo de 3)
    4: "TGTCCGGCGAGGGCGAGGGCGATGCCACC",                     # 29 nt
    5: "TACGGCAAGCTGACCCTGAAGTTCATCTGCTAA",                 # 33 nt  cola + STOP
}
ISOFORMAS = {"canónica 1-2-3-4-5": [1, 2, 3, 4, 5],
             "salto E3: 1-2-4-5": [1, 2, 4, 5],
             "salto E2+E4: 1-3-5": [1, 3, 5]}


def analizar(nombre, orden):
    cds = "".join(EXONES[e] for e in orden)
    prot = Seq(cds[:len(cds) - len(cds) % 3]).translate(to_stop=True)
    marco = "conservado" if len(cds) % 3 == 0 else "DESPLAZADO"
    print(f"{nombre:<22} {len(cds):>4} nt  marco {marco:<10} {len(prot):>3} aa  {prot}")
    return prot


if __name__ == "__main__":
    print("== Ejercicio 4: splicing alternativo (gen modelo) ==")
    for n, o in ISOFORMAS.items():
        analizar(n, o)

    tsv = RAIZ / "data" / "FGFR2_transcritos_ensembl.tsv"
    print("\n== FGFR2 (Ensembl ENSG00000066468, cromosoma 10, hebra -) ==")
    with open(tsv, encoding="utf-8") as f:
        filas = list(csv.DictReader(f, delimiter="\t"))
    for r in filas:
        print(f"{r['transcrito']}  {r['nombre']:<10} {r['biotipo']:<24} exones={r['n_exones']:>2}  prot={r['longitud_proteina_aa']:>4} aa  canónico={r['canonico']}")
    cod = [int(r["longitud_proteina_aa"]) for r in filas if r["biotipo"] == "protein_coding"]
    print(f"\nIsoformas codificantes mostradas: {len(cod)} | proteína de {min(cod)} a {max(cod)} aa")
