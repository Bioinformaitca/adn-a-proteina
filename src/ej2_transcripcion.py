"""Ejercicio 2 - Transcripción ADN -> ARNm a partir de un FASTA.

Uso: python ej2_transcripcion.py [fichero.fasta]
La secuencia del FASTA se interpreta como hebra CODIFICANTE 5'->3'.
El script muestra la cadena molde, el ARNm y qué pasa si se usa la
orientación equivocada de la hebra.
"""
import sys
from pathlib import Path
from Bio import SeqIO

RAIZ = Path(__file__).resolve().parent.parent


def transcribir_desde_molde(molde_3_5: str) -> str:
    """La ARN polimerasa lee el molde 3'->5' y sintetiza ARN 5'->3' (U en vez de T)."""
    par = {"A": "U", "T": "A", "G": "C", "C": "G"}
    return "".join(par[b] for b in molde_3_5)


def fmt(s):
    return " ".join(s[i:i + 3] for i in range(0, len(s), 3))


if __name__ == "__main__":
    fasta = Path(sys.argv[1]) if len(sys.argv) > 1 else RAIZ / "data" / "ej2_adn.fasta"
    for rec in SeqIO.parse(fasta, "fasta"):
        cod = rec.seq.upper()
        molde = cod.complement()                     # 3'->5', alineada
        print(f"== {rec.id} ==")
        print(f"Hebra codificante: 5'-{fmt(str(cod))}-3'")
        print(f"Hebra molde      : 3'-{fmt(str(molde))}-5'")

        arn_manual = transcribir_desde_molde(str(molde))
        arn_bio = str(cod.transcribe())              # Biopython: codificante con T->U
        print(f"ARNm (desde molde): 5'-{fmt(arn_manual)}-3'")
        print(f"ARNm (Bio.Seq)    : 5'-{fmt(arn_bio)}-3'  coincide={arn_manual == arn_bio}")

        # Experimento: orientación cambiada -> usar la otra hebra como codificante
        otra = cod.reverse_complement()
        arn_otra = str(otra.transcribe())
        print("\nExperimento: si tomamos la otra hebra como codificante (reverse complement):")
        print(f"  ARNm alternativo: 5'-{fmt(arn_otra)}-3'")
        print(f"  Traducción      : {otra.transcribe().translate()}  (vs correcta: {cod.transcribe().translate()})")
        print(f"  ¿Empieza por AUG? {arn_otra.startswith('AUG')} -> producto distinto / sin sentido\n")
