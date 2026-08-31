# Generate a single RaCInG graph realization

Generate a single RaCInG graph realization

## Usage

``` r
model1(
  N,
  avdeg,
  cellLigList,
  cellRecList,
  Dcelltype,
  Dligrec,
  Signmatrix = NULL
)
```

## Arguments

- N:

  Number of vertices (cells) in the graph.

- avdeg:

  Target average degree.

- cellLigList:

  Cell-by-ligand compatibility matrix.

- cellRecList:

  Cell-by-receptor compatibility matrix.

- Dcelltype:

  Cell-type abundance probabilities.

- Dligrec:

  Ligand-by-receptor probability matrix.

- Signmatrix:

  Optional ligand-receptor sign matrix.

## Value

A list with vertex labels, an edge list, and ligand-receptor types.
