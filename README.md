This package contains R functions to perform linear assignment problems using Jonker-Vogenant algorithms.

The package uses the C++ code from the [Python LAPJV package](https://github.com/Bram94/lapjv).

## Installation

```
devtools::install_github("dpmcsuss/lapjv")
```

## Use

Currently the `lapjv` function has been tested and compared to the output of the CLUE package and yields identical results (hopefully :-) ).

```
set.seed(12345)
library(lapjv)
x <- matrix(runif(100), 10)
lapjv(x)
#  [1]  8  3  9  1  2  4  6  7 10  5
clue::solve_LSAP(x)
# Optimal assignment:
# 1 => 8, 2 => 3, 3 => 9, 4 => 1, 5 => 2, 6 => 4, 7 => 6, 8 => 7,
# 9 => 10, 10 => 5

lapjv(x, maximize = TRUE)
# [1]  6  9  4  3  7  5  2  1  8 10
clue::solve_LSAP(x, maximum = TRUE)
# Optimal assignment:
# 1 => 6, 2 => 9, 3 => 4, 4 => 3, 5 => 7, 6 => 5, 7 => 2, 8 => 1,
# 9 => 8, 10 => 10
```



## License

BSD 2-clause.