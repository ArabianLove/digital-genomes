#!/usr/bin/env python3
"""
Karyon Core – Digital Genomes & Communities
Executable genomic engine (baseline v0.9)
Transliteration of human genetic architecture into autonomous digital units.
Apache-2.0
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List
import random
import math

# ── Alphabet & genetic code ──
BASES = "ACGT"
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate(cds: str) -> str:
    protein = []
    for i in range(0, len(cds) - 2, 3):
        aa = CODON_TABLE.get(cds[i:i+3].upper(), 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)

def hill(x: float, k: float = 0.5, n: float = 2.0) -> float:
    return (x ** n) / (k ** n + x ** n) if (k ** n + x ** n) > 0 else 0.0

def motif_count(seq: str, motif: str) -> int:
    return seq.upper().count(motif.upper()) if motif else 0

# ── Real cis-motifs & landmark synteny ──
MOTIFS = {
    "HOX": "TAAT",
    "MYC": "CACGTG",
    "TP53": "AGGCAT",
    "SRY": "AACAAT",
}

@dataclass
class Gene:
    name: str
    chrom: str
    promoter: str
    cds: str
    motif: str = ""
    role: str = "effector"

SINTENIA: List[Gene] = [
    Gene("HOXA1",  "7",  "TAATTAATGC" * 3, "ATGGCCGAG" + "GCG" * 20 + "TAA", "TAAT", "TF"),
    Gene("HOXA13", "7",  "TAATGC" * 4,     "ATGAGCGAG" + "CCG" * 15 + "TGA", "TAAT", "TF"),
    Gene("TP53",   "17", "AGGCATAGGCAT",   "ATG" + "GAG" * 30 + "TAA",       "AGGCAT", "checkpoint"),
    Gene("MYC",    "8",  "CACGTGCACGTG",   "ATG" + "CCC" * 25 + "TAG",       "CACGTG", "TF"),
    Gene("HLA-A",  "6",  "GGGGCGGGG",      "ATG" + "GCG" * 40 + "TAA",       "", "immune"),
    Gene("HTT",    "4",  "CGCGCG",         "ATG" + ("CAG" * 28) + "CCG" * 10 + "TAA", "", "polyQ"),
    Gene("IGF2",   "11", "GGGG",           "ATG" + "GAA" * 20 + "TGA",       "", "imprinted"),
    Gene("SNRPN",  "15", "CGCG",           "ATG" + "AGC" * 15 + "TAA",       "", "imprinted"),
    Gene("UBE3A",  "15", "TAAT",           "ATG" + "GAC" * 18 + "TAG",       "", "imprinted"),
    Gene("MT-ATP6","MT", "",               "ATG" + "TTA" * 20 + "TAA",       "", "metabolic"),
]

@dataclass
class Allele:
    gene: Gene
    promoter: str
    cds: str
    methylation: float = 0.15
    accessibility: float = 0.85

@dataclass
class Organism:
    id: str
    alleles: List[Allele]
    age: float = 0.0
    energy: float = 12.0
    energy_cap: float = 15.0
    proteome: Dict[str, float] = field(default_factory=dict)
    alive: bool = True
    pathology: List[str] = field(default_factory=list)

    def step(self, dt: float = 0.2):
        if not self.alive:
            return
        tfs = {k: v for k, v in self.proteome.items() if k in ("HOXA1", "MYC", "TP53")}
        for al in self.alleles:
            rate = 0.05
            for tf_name, level in tfs.items():
                motif = MOTIFS.get("HOX" if "HOX" in tf_name else tf_name, "")
                if motif and motif in al.promoter.upper():
                    rate += hill(level) * motif_count(al.promoter, motif) * al.accessibility * (1.0 - al.methylation)
            current = self.proteome.get(al.gene.name, 0.1)
            delta = rate * 0.3 - current * 0.1
            self.proteome[al.gene.name] = max(0.0, current + delta * dt)

        mt = self.proteome.get("MT-ATP6", 0.5)
        self.energy_cap = 8.0 + 10.0 * min(mt, 1.5)
        self.energy = min(self.energy + mt * 0.4 * dt - 0.15 * dt, self.energy_cap)
        self.age += dt

        self.pathology = []
        htt = next((a.cds for a in self.alleles if a.gene.name == "HTT"), "")
        cag = htt.upper().count("CAG")
        if cag > 35:
            self.pathology.append(f"polyQ_HTT(CAG={cag})")
        if self.proteome.get("TP53", 0) < 0.15 and self.age > 5:
            self.pathology.append("guardian_loss")
        if self.energy < 1.0 or self.age > 80:
            self.alive = False

def make_founder(oid: str) -> Organism:
    alleles = [
        Allele(g, g.promoter, g.cds,
               methylation=random.uniform(0.05, 0.25),
               accessibility=random.uniform(0.7, 0.95))
        for g in SINTENIA
    ]
    org = Organism(id=oid, alleles=alleles)
    for _ in range(24):
        org.step(0.2)
    return org

def main():
    random.seed(42)
    print("=" * 60)
    print("KARYON v0.9 – Digital Genomes & Communities")
    print("Executable genomic engine (transliteration baseline)")
    print("=" * 60)

    pop = [make_founder(f"K{i:03d}") for i in range(6)]
    print(f"\nGenesis: {len(pop)} founders created")
    print(f"Founder0 energy={pop[0].energy:.2f}  age={pop[0].age:.1f}")
    print(f"Proteome sample: { {k: round(v, 2) for k, v in list(pop[0].proteome.items())[:5]} }")

    for t in range(40):
        for o in pop:
            o.step(0.25)

    alive = sum(1 for o in pop if o.alive)
    print(f"\nAfter 40 ticks: {alive}/{len(pop)} alive")
    print(f"Sample pathologies: {pop[0].pathology}")
    print(f"HTT CAG count (founder0): {pop[0].alleles[5].cds.count('CAG')}")
    print("\nCore OK – genome, GRN (Hill), metabolism, pathology running.")
    print("=" * 60)

if __name__ == "__main__":
    main()
