# Extract multiple communication feature types from the same simulated graphs

Generates the requested number of Monte Carlo graph realizations *once*
per patient and extracts every requested feature family from each
realization, instead of the alternative of calling
[`countDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countDirect.md),
[`countWedges()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countWedges.md),
etc. separately – which would independently re-simulate a brand new set
of random graphs per feature type. Graph generation
([`model1()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/model1.md)'s
vertex/edge sampling), not feature extraction, is the expensive part of
a Monte Carlo run, and it does not depend on which feature you
eventually want – so sharing one set of realizations across every
requested type avoids paying that cost once per type.

## Usage

``` r
countAllTypes(
  Dcell,
  Dconn,
  lig,
  rec,
  cellnames,
  N,
  av,
  itNo,
  types = c("D", "W", "TT", "CT", "GSCC")
)
```

## Arguments

- Dcell:

  Cell-type abundance vector for one patient.

- Dconn:

  Ligand-receptor probability matrix for one patient.

- lig:

  Cell-by-ligand compatibility matrix.

- rec:

  Cell-by-receptor compatibility matrix.

- cellnames:

  Character vector of cell-type names.

- N:

  Number of cells per simulated graph.

- av:

  Target average degree.

- itNo:

  Number of Monte Carlo iterations.

- types:

  Character vector of feature families to extract from each simulated
  graph: any combination of `"D"`, `"W"`, `"TT"`, `"CT"`, `"GSCC"`.

## Value

A named list, one element per requested entry in `types`, each with the
same structure the corresponding single-type `count*()` function returns
(e.g. `out$D` matches
[`countDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/countDirect.md)'s
return value).
