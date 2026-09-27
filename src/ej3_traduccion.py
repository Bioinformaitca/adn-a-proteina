"""Ejercicio 3 - Traducción ARNm -> proteína y efecto de mutaciones."""
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable

TABLA = CodonTable.unambiguous_rna_by_id[1]      # código genético estándar


def codones(arn):
    return [arn[i:i + 3] for i in range(0, len(arn) - len(arn) % 3, 3)]


def traducir_manual(arn):
    """Busca el primer AUG y traduce codón a codón hasta un codón de paro."""
    ini = arn.find("AUG")
    if ini == -1:
        return None, "sin codón de inicio AUG: no hay traducción (eucariota)"
    prot = []
    for c in codones(arn[ini:]):
        if c in TABLA.stop_codons:
            return "-".join(prot), f"paro en {c}"
        prot.append(seq3(TABLA.forward_table[c]))
    return "-".join(prot), "sin codón de paro: lectura continuaría (read-through)"


if __name__ == "__main__":
    arn = "AUGUAUGCUUAA"
    print("== Ejercicio 3: traducción ==")
    print("Codones:", " ".join(codones(arn)))
    print("Inicio :", codones(arn)[0], "| Paro:", codones(arn)[-1])
    manual, nota = traducir_manual(arn)
    bio = Seq(arn).translate(to_stop=True)
    print(f"Manual : {manual}  ({nota})")
    print(f"Bio.Seq: {seq3(str(bio))}  -> 1 letra: {bio}  coincide={seq3(str(bio)) == manual.replace('-', '')}")

    print("\nMutación AUG -> GUG en el inicio:")
    mut = "GUG" + arn[3:]
    print("  Eucariota (busca AUG):", traducir_manual(mut))
    print("  Si GUG se usa como inicio (bacterias), el iniciador pone fMet igualmente;"
          " traducción 'cruda' por Bio.Seq:", Seq(mut).translate(to_stop=True))

    print("\nMutación del codón de paro UAA -> CAA (pérdida de stop):")
    mut2 = arn[:-3] + "CAA" + "GGCUUUAG"   # continúa en la región 3' no traducida
    print("  ", traducir_manual(mut2))
