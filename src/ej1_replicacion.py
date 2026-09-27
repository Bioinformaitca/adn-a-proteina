"""Ejercicio 1 - Replicación del ADN (modelo semiconservativo).

Dada una hebra 5'->3', genera su complementaria y simula una ronda de
replicación: cada hebra parental sirve de molde para una hebra nueva.
Compara el resultado automático con el obtenido a mano.
"""
from Bio.Seq import Seq


def complementaria(hebra_5_3: str) -> str:
    """Devuelve la hebra complementaria escrita 3'->5' (alineada base a base)."""
    return str(Seq(hebra_5_3).complement())


def replicar(hebra_superior: str):
    """Simula una ronda de replicación y devuelve las dos moléculas hijas."""
    hebra_inferior = complementaria(hebra_superior)            # 3'->5'
    nueva_sobre_superior = complementaria(hebra_superior)      # 3'->5'
    nueva_sobre_inferior = complementaria(hebra_inferior)      # 5'->3'
    hija1 = (hebra_superior, nueva_sobre_superior)             # parental + nueva
    hija2 = (nueva_sobre_inferior, hebra_inferior)             # nueva + parental
    return hija1, hija2


def fmt(s: str) -> str:
    return " ".join(s[i:i + 3] for i in range(0, len(s), 3))


if __name__ == "__main__":
    superior = "ATGCCGTTAGCT"
    manual_inferior = "TACGGCAATCGA"   # la dada en el enunciado

    print("== Ejercicio 1: replicación ==")
    print(f"Molde parental 1: 5'-{fmt(superior)}-3'")
    auto = complementaria(superior)
    print(f"Complementaria  : 3'-{fmt(auto)}-5'")
    print(f"¿Coincide con la del enunciado/manual? {auto == manual_inferior}")
    print("Reverse complement (leída 5'->3'):", fmt(str(Seq(superior).reverse_complement())))

    h1, h2 = replicar(superior)
    print("\nMolécula hija 1 (parental arriba, NUEVA abajo):")
    print(f"  5'-{fmt(h1[0])}-3'\n  3'-{fmt(h1[1])}-5'  <- nueva")
    print("Molécula hija 2 (NUEVA arriba, parental abajo):")
    print(f"  5'-{fmt(h2[0])}-3'  <- nueva\n  3'-{fmt(h2[1])}-5'")
    print("\nCada hija conserva una hebra parental -> replicación semiconservativa.")
