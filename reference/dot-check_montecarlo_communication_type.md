# Check that only Monte Carlo feature families still in RaCInG are requested

Check that only Monte Carlo feature families still in RaCInG are
requested

## Usage

``` r
.check_montecarlo_communication_type(communication_type)
```

## Arguments

- communication_type:

  Requested feature family/families.

## Value

Invisibly `TRUE`; errors for anything other than `"D"`/`"GSCC"`.
