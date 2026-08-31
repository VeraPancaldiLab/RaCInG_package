# Reset foreach backend to sequential

Internal helper that unregisters any active parallel backend used by
foreach and restores the sequential backend via
[`foreach::registerDoSEQ()`](https://rdrr.io/pkg/foreach/man/registerDoSEQ.html),
mirroring `pipeML`'s helper of the same purpose. Called after
[`parallel::stopCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
to make sure the cluster is fully released.

## Usage

``` r
.unregister_dopar()
```
