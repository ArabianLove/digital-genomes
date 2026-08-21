#!/usr/bin/env python3
"""
Karyon — Core genomico eseguibile
Traslitterazione funzionale (non mera simulazione) del genoma umano
in unità relazionali digitali.

Alfabeto DNA, codice degenerato, sintenia landmark, motivi cis reali,
GRN con cinetica di Hill, epigenoma, organismo con metabolismo,
patologie, senescenza e collasso autopoietico.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import random
import math

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

MOTIFS = {
    "HOX": "TAAT",
    "MYC": "CACGTG",
    "TP53": "AGGCAT",
    "SRY": "AACAAT",
}

def random_dna(length: int, gc: float = 0.42) -> str:
    p_gc = gc / 2
    p_at = (1 - gc) / 2
    weights = [p_at, p_gc, p_gc, p_at]
    return "".join(random.choices(BASES, weights=weights, k=length))

def motif_count(seq: str, motif: str) -> int:
    count = 0
    i = 0
    m = motif.upper()
    s = seq.upper()
    while True:
        i = s.find(m, i)
        if i == -1:
            break
        count += 1
        i += len(m)
    return count

def hill(x: float, k: float = 0.5, n: float = 2.0) -> float:
    x = max(0.0, x)
    return (x ** n) / (k ** n + x ** n + 1e-12)

@dataclass
class Gene:
    name: str
    promoter: str
    cds: str
    motif: str = "TAAT"

@dataclass
class Organism:
    id: str
    genes: Dict[str, Gene]
    protein: Dict[str, float] = field(default_factory=dict)
    energy: float = 12.0
    energy_cap: float = 22.0
    age: float = 0.0
    telomere: float = 1.0
    pathology: str = ""

    def ensure_proteins(self):
        for name in self.genes:
            self.protein.setdefault(name, 0.3)

    def step(self, dt: float = 1.0, resource: float = 1.0):
        self.ensure_proteins()
        new_p = {}
        for name, gene in self.genes.items():
            tf = self.protein.get(name, 0.3)
            rate = 0.05 + hill(tf) * motif_count(gene.promoter, gene.motif) * 0.12
            dp = rate - 0.08 * self.protein[name]
            new_p[name] = max(0.0, self.protein[name] + dt * 0.2 * dp)
        self.protein.update(new_p)

        mt = self.protein.get("MT-ATP6", 0.3)
        self.energy_cap = 8.0 + 14.0 * min(1.2, mt)
        self.energy += 0.4 * resource * mt
        self.energy -= 0.15 * dt
        self.energy = max(0.0, min(self.energy, self.energy_cap))

        self.age += dt
        self.telomere = max(0.0, self.telomere - 0.0015 * dt)

        # pathologies
        if self.protein.get("TP53", 1.0) < 0.15 and self.age > 20:
            self.pathology = "guardian_loss"
        htt = self.genes.get("HTT")
        if htt and motif_count(htt.cds, "CAG") >= 6:
            self.pathology = "polyQ"
        if self.energy < 0.1 or self.telomere < 0.05:
            self.pathology = "autopoietic_collapse"

    def summary(self):
        return {
            "id": self.id,
            "age": round(self.age, 1),
            "energy": round(self.energy, 2),
            "telomere": round(self.telomere, 3),
            "pathology": self.pathology,
            "TP53": round(self.protein.get("TP53", 0), 2),
            "MYC": round(self.protein.get("MYC", 0), 2),
        }

def make_gene(name: str, motif: str, gc: float = 0.45) -> Gene:
    promoter = random_dna(40, gc) + motif + random_dna(30, gc)
    cds = random_dna(90, gc)
    if name == "HTT":
        cds = "CAG" * 4 + cds  # base polyQ
    return Gene(name=name, promoter=promoter, cds=cds, motif=motif)

def genesis(n: int = 6, seed: int = 42) -> List[Organism]:
    rng = random.Random(seed)
    random.seed(seed)
    landmarks = [
        ("HOXA1", "TAAT"),
        ("TP53", "AGGCAT"),
        ("MYC", "CACGTG"),
        ("HTT", "CAG"),
        ("MT-ATP6", "ATGC"),
        ("HLA-A", "TGGAAA"),
    ]
    orgs = []
    for i in range(n):
        genes = {name: make_gene(name, motif) for name, motif in landmarks}
        org = Organism(id=f"K-{i:03d}", genes=genes)
        # warmup
        for _ in range(8):
            org.step(dt=0.3, resource=1.0)
        orgs.append(org)
    return orgs

def run_demo(ticks: int = 15):
    print("=" * 60)
    print("KARYON — Core genomico eseguibile")
    print("Traslitterazione (non simulazione) — Demo")
    print("=" * 60)
    population = genesis(6)
    print(f"Genesi: {len(population)} organismi fondatori")
    for t in range(ticks):
        for org in population:
            org.step(dt=1.0, resource=0.9 + 0.2 * random.random())
        alive = [o for o in population if o.energy > 0.05]
        mean_e = sum(o.energy for o in alive) / max(1, len(alive))
        print(f"t={t:02d}  vivi={len(alive)}  meanE={mean_e:.2f}")
    print("\nFenotipi finali:")
    for org in population:
        s = org.summary()
        print(f"  {s['id']}  age={s['age']:5.1f}  E={s['energy']:5.2f}  "
              f"tel={s['telomere']:.3f}  path={s['pathology'] or '-'}  "
              f"TP53={s['TP53']:.2f} MYC={s['MYC']:.2f}")

if __name__ == "__main__":
    run_demo()
