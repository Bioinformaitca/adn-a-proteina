"""Ejercicio 5 - Secuencia, estructura y función de proteínas.

1) Péptido del enunciado: extremos N/C e hidropatía (escala Kyte-Doolittle).
2) Ubiquitina (PDB 1UBQ): mapa de hélices alfa y láminas beta y simulación
   de una mutación hidrofóbico -> hidrofílico en el núcleo (Val26, residuo enterrado de la hélice).
"""
from pathlib import Path
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis

RAIZ = Path(__file__).resolve().parent.parent
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


def leer_ss(ruta):
    seq, ss = "", []
    for linea in open(ruta, encoding="utf-8"):
        p = linea.split()
        if not p or p[0].startswith("#"):
            continue
        if p[0] == "SEQ":
            seq = p[1]
        else:
            ss.append((p[0], int(p[1]), int(p[2])))
    mapa = ["-"] * len(seq)
    for tipo, a, b in ss:
        for i in range(a - 1, b):
            mapa[i] = "H" if tipo == "HELIX" else "E"
    return seq, "".join(mapa)


if __name__ == "__main__":
    pep3 = ["Met", "Ile", "Ser", "Gly", "Val", "Lys", "His"]
    pep = seq1("".join(pep3))
    print("== Ejercicio 5: péptido ==")
    print("H2N-(N-terminal) " + "-".join(pep3) + " (C-terminal)-COOH")
    for a3, a1 in zip(pep3, pep):
        tipo = {"M": "apolar/hidrofóbico", "I": "apolar/hidrofóbico", "V": "apolar/hidrofóbico",
                "G": "apolar (cadena lateral = H)", "S": "polar sin carga",
                "K": "básico, carga +", "H": "básico (carga + parcial)"}[a1]
        print(f"  {a3}: KD={KD[a1]:+.1f}  {tipo}")
    pa = ProteinAnalysis(pep)
    print(f"GRAVY={pa.gravy():.2f}  pI={pa.isoelectric_point():.2f}  MW={pa.molecular_weight():.1f} Da")

    seq, ss = leer_ss(RAIZ / "data" / "1UBQ_estructura_secundaria.txt")
    print("\n== Ubiquitina 1UBQ (H=hélice alfa, E=lámina beta) ==")
    for i in range(0, len(seq), 40):
        print(f"{i+1:>3} {seq[i:i+40]}\n    {ss[i:i+40]}")
    print(f"Hélice: {ss.count('H')} res ({ss.count('H')/len(seq):.0%}) | Lámina: {ss.count('E')} res ({ss.count('E')/len(seq):.0%})")

    # Mutación en el núcleo hidrofóbico: Val26 (hélice) -> Lys
    pos = 26
    mut = seq[:pos - 1] + "K" + seq[pos:]
    print(f"\nMutación puntual V{pos}K (residuo en '{ss[pos-1]}'): hidropatía {KD['V']:+.1f} -> {KD['K']:+.1f}")
    print(f"GRAVY proteína: {ProteinAnalysis(seq).gravy():.3f} -> {ProteinAnalysis(mut).gravy():.3f}")
    print("Una carga (+) enterrada en el núcleo desestabiliza el plegamiento aunque el cambio global sea pequeño.")
