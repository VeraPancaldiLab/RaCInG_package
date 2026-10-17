# Convert raw simulation outputs into feature matrices

Convert raw simulation outputs into feature matrices

## Usage

``` r
compute_results_processing(
  celltypes,
  patient_names,
  communication_type,
  remove_direction = TRUE,
  normalized = !is.null(sim_norm_file),
  sim_raw_file = NULL,
  sim_norm_file = NULL,
  output_folder = NULL,
  file.name = NULL
)
```

## Arguments

- celltypes:

  Character vector of cell-type names.

- patient_names:

  Character vector of patient identifiers: either one per row of the
  simulated `Cmatrix` (matched through the patient numbers stored in the
  `.out` file, so subset runs get the right names), or exactly one per
  patient in the file, in file order.

- communication_type:

  Communication type of the simulation files (`"D"` or `"GSCC"`).

- remove_direction:

  Logical; if `TRUE`, merge the two directions of each cell-type pair
  (direct communication only).

- normalized:

  Logical; if `TRUE`, divide raw values by the uniform baseline in
  `sim_norm_file`. Defaults to `TRUE` exactly when `sim_norm_file` is
  supplied. A feature that is zero in both simulations (0/0) is reported
  as `0`, as in the kernel method's
  [`calculateDirect()`](https://VeraPancaldiLab.github.io/RaCInG_package/reference/calculateDirect.md).

- sim_raw_file:

  Path to the raw simulation output file.

- sim_norm_file:

  Optional path to the normalized simulation output file.

- output_folder:

  Directory where the processed CSV files are written.

- file.name:

  Output CSV file name for the means (e.g. `"RaCInG_input.csv"`); the
  standard deviations go to the same name with an `_sd` suffix
  (`"RaCInG_input_sd.csv"`).

## Value

A list with two patient-by-feature data frames, `mean` and `sd`. Direct
features are named `"Dir_A_B"`, GSCC features `"GSCC_A"`; rows are named
by patient. Patients whose simulation produced no edges are dropped.
