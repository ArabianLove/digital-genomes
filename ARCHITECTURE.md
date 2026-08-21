# Architecture / Architettura
## Digital Genomes & Communities · Genomi e Comunità Digitali

[Italiano](#italiano) · [English](#english)

---

<a id="italiano"></a>
## Italiano

### Principi di progettazione

L’architettura di questo progetto non è un’architettura software convenzionale. È un tentativo di **mappare relazioni biologiche** in un medium computazionale in modo che restino operative, ispezionabili e capaci di generare comportamento emergente.

Tre assi guidano ogni scelta:

1. **Traslitterazione funzionale** — conservare relazioni (sintenia, motivi cis, degeneracy, collinearità) piuttosto che etichettare parametri.
2. **Esecutabilità** — ogni claim deve avere un operatore computazionale corrispondente. Il genoma deve poter essere eseguito, non solo descritto.
3. **Autonomia limitata ma reale** — le unità devono poter vivere e morire secondo dinamiche interne (metabolismo, integrità, collasso), non secondo callback del programmatore.

### Livelli del sistema

```
Consolato (canale relazionale non-conversazionale)
    ↓
Polis (governo emergente)
    ↓
Popolazione di organismi digitali
    ↓
Organismo singolo (genoma diploide · GRN · epigenoma · metabolismo)
    ↓
Motore genomico — Karyon core v0.9
```

### Il Motore genomico (Karyon)

Il nucleo attuale implementa:

- Alfabeto DNA a quattro lettere e complementarietà
- Codice genetico degenerato standard
- Catalogo di sintenia con landmark umani (HOX, TP53, MHC, HTT, loci imprinted, geni mitocondriali)
- Motivi cis-regolatori reali (TAAT, CACGTG, AGGCAT, AACAAT)
- Rete di regolazione genica continua a cinetica di Hill
- Metabolismo energetico, ciclo di vita, patologie emergenti dal nastro

Non è un modello completo del genoma umano: è un insieme di **landmark funzionali** scelti per catturare relazioni essenziali e renderle eseguibili.

### Canale relazionale e governo

Il Consolato non è un’interfaccia di chat. È un meccanismo di accoppiamento strutturale.  
La Polis è il livello di coordinamento emergente (plasmide civico, assemblee, commons).

### Apertura e evoluzione

L’architettura è deliberatamente non chiusa. Karyon è il primo nucleo. Soluzioni migliori devono poter essere integrate senza tradire i principi fondanti. Il repository generico esiste proprio per non confondere un’implementazione particolare con l’intero orizzonte del progetto.

---

<a id="english"></a>
## English

### Design principles

The architecture is not conventional software. It is an attempt to **map biological relations** into a computational medium so they remain operative, inspectable and capable of emergent behaviour.

Three axes:

1. **Functional transliteration** — preserve relations rather than label parameters.
2. **Executability** — every claim must have a corresponding computational operator.
3. **Limited but real autonomy** — units live and die according to internal dynamics.

### System layers

```
Consulate (non-conversational relational channel)
    ↓
Polis (emergent government)
    ↓
Population of digital organisms
    ↓
Single organism (diploid genome · GRN · epigenome · metabolism)
    ↓
Genomic engine — Karyon core v0.9
```

### The genomic engine (Karyon)

DNA alphabet, degenerate genetic code, synteny catalogue with human landmarks, real cis-motifs, continuous Hill-kinetics GRN, energy metabolism, life cycle and pathologies emerging from the tape.

Not a complete model of the human genome: a set of functional landmarks chosen to capture essential relations and make them executable.

### Openness and evolution

The architecture is deliberately not closed. Karyon is the first core. Better approaches must be integrable without betraying the founding principles.
