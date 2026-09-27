"""Ejercicio 6 - Pipeline del dogma central: replicación -> transcripción -> traducción.

Uso:
    python ej6_pipeline.py [fichero.fasta] [--salida DIR]

Por defecto usa data/HBB_cds.fasta (CDS del gen HBB humano, beta-globina,
NCBI RefSeq NM_000518.5). Informa por consola (logging) de cada paso y
guarda los resultados en FASTA en la carpeta de salida.
"""
import argparse
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

RAIZ = Path(__file__).resolve().parent.parent
log = logging.getLogger("dogma")
VALIDAS = set("ACGT")


def validar(seq: Seq) -> Seq:
    seq = seq.upper()
    raras = set(str(seq)) - VALIDAS
    if raras:
        raise ValueError(f"Bases no válidas en la secuencia: {sorted(raras)}")
    log.info("Validación OK: %d nt, GC=%.1f%%", len(seq),
             100 * (seq.count("G") + seq.count("C")) / len(seq))
    return seq


def replicacion(cod: Seq):
    """Paso 1: cada hebra parental es molde de una hebra hija nueva."""
    log.info("[1/3 REPLICACIÓN] Helicasa abre la doble hélice (%d pb)", len(cod))
    molde = cod.reverse_complement()                 # hebra antisentido, leída 5'->3'
    log.info("  Hebra parental 2 (antisentido, 5'->3'): %s...", molde[:30])
    nueva_1 = cod.reverse_complement()               # sintetizada sobre la codificante
    nueva_2 = molde.reverse_complement()             # sintetizada sobre la antisentido
    log.info("  ADN polimerasa sintetiza 5'->3' la hebra nueva sobre la parental 1: %s...", nueva_1[:30])
    log.info("  ADN polimerasa sintetiza 5'->3' la hebra nueva sobre la parental 2: %s...", nueva_2[:30])
    assert nueva_1 == molde and nueva_2 == cod, "Las hijas no son idénticas a la parental"
    log.info("  Ligasa une fragmentos de Okazaki. 2 moléculas hijas idénticas (semiconservativa) ✔")
    return nueva_1, nueva_2, molde


def transcripcion(cod: Seq, molde: Seq) -> Seq:
    """Paso 2: la ARN polimerasa lee el molde 3'->5' y produce ARNm 5'->3'."""
    log.info("[2/3 TRANSCRIPCIÓN] Cadena molde = hebra antisentido")
    arnm = molde.reverse_complement_rna()            # complemento del molde con U
    assert arnm == cod.transcribe()
    log.info("  ARNm (%d nt): %s...", len(arnm), arnm[:30])
    log.info("  Comprobación: ARNm == codificante con T->U ✔")
    return arnm


def traduccion(arnm: Seq) -> Seq:
    """Paso 3: el ribosoma lee codones desde el primer AUG hasta un stop."""
    log.info("[3/3 TRADUCCIÓN] Buscando codón de inicio AUG")
    ini = arnm.find("AUG")
    if ini == -1:
        log.warning("  No hay AUG: no se traduce")
        return Seq("")
    orf = arnm[ini:]
    prot = orf.translate(to_stop=True)
    fin = ini + 3 * len(prot)
    stop = arnm[fin:fin + 3]
    log.info("  Inicio en posición %d; %d codones leídos; codón de paro: %s",
             ini + 1, len(prot) + 1, stop if len(stop) == 3 else "NINGUNO")
    log.info("  Proteína (%d aa): %s", len(prot), prot)
    return prot


def demo_mutacion(cod: Seq, pos: int = 20, nueva: str = "T"):
    """Introduce un error puntual en el ADN (p. ej. fallo no corregido de la
    ADN polimerasa) y lo propaga por transcripción y traducción."""
    mut = cod[:pos - 1] + nueva + cod[pos:]
    p0, p1 = cod.transcribe().translate(to_stop=True), mut.transcribe().translate(to_stop=True)
    diffs = [f"{a}{i + 1}{b}" for i, (a, b) in enumerate(zip(p0, p1)) if a != b]
    log.info("[DEMO MUTACIÓN] ADN pos %d: %s->%s | codón %s->%s | cambio proteico: %s",
             pos, cod[pos - 1], nueva, cod[(pos - 1) // 3 * 3:(pos - 1) // 3 * 3 + 3],
             mut[(pos - 1) // 3 * 3:(pos - 1) // 3 * 3 + 3], diffs or "ninguno (silenciosa)")
    return diffs


def pipeline(fasta: Path, salida: Path, demo: bool = False):
    salida.mkdir(parents=True, exist_ok=True)
    for rec in SeqIO.parse(fasta, "fasta"):
        log.info("=== Registro %s ===", rec.id)
        cod = validar(rec.seq)
        h1, h2, molde = replicacion(cod)
        arnm = transcripcion(cod, molde)
        prot = traduccion(arnm)
        regs = [SeqRecord(h1, id=f"{rec.id}_nueva_hebra1", description="replicación"),
                SeqRecord(h2, id=f"{rec.id}_nueva_hebra2", description="replicación"),
                SeqRecord(arnm, id=f"{rec.id}_ARNm", description="transcripción"),
                SeqRecord(prot, id=f"{rec.id}_proteina", description="traducción")]
        out = salida / f"{rec.id.split('|')[-1]}_resultados.fasta"
        SeqIO.write(regs, out, "fasta")
        if demo:
            demo_mutacion(cod)
        log.info("Resultados guardados en %s", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", nargs="?", default=RAIZ / "data" / "HBB_cds.fasta", type=Path)
    ap.add_argument("--salida", default=RAIZ / "resultados", type=Path)
    ap.add_argument("--demo-mutacion", action="store_true",
                    help="propaga una mutación puntual A20T (en HBB = mutación de la anemia falciforme)")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s | %(message)s")
    pipeline(a.fasta, a.salida, a.demo_mutacion)
