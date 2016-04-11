## Introduction
I will write the introduction to this project later when I am free.

## Build
Change `compiler_set` and `addition_includes` in `SConstruct` according to your
actual compiler.

run `scons` to build.

## Run
The recommended way to run the programs is to create separate folders for each
cases. The following is a general work flow in normal cases:

```
  cp -r data/ ~/
  mkdir -p job/01
  cd job
  vim README.md   # write some comments about case 01
  cd 01
  cp ../../model/marm/* .
  vim SConstruct  # update the parameters is necessary
  vim run.sh      # choose the task you want to start
  ./run.sh
```

If you want to run `enfwi` or `enfwi-sw`, you will need perturbation files as
input. You can generate this file by compiling and runing the `

## Show result
In the same directory where you run the case, you can check the result including
objective function values, L1 L2 normalized model fit and inverted images.

- you can display the numeric values of absolute objective function values by
  `sfdisfil < absobjs.rsf`. The `sfdisfil` command also works for `norobjs.rsf`,
  `l1norm.rsf` and `l2norm.rsf`.

- you can display the plots by running `sfpen < absobjs.vpl`, `sfpen vel300.vpl`
  ...

- of course, you can also checkout the log manually.

## SWC++ compiler issues
- we can only pass the `-O0` optimized flags to the `mpicxx -ver 5.421-sw-437 -host`. Compile it with `-O2` will creash the program.

- the follow codes in `src/common/parabola-vertex.cpp`

```
  if (std::abs(k2 - k1) < 0.001 * (std::max(std::abs(k2), std::abs(k1))) ||
      (xv == -std::numeric_limits<double>::quiet_NaN())) {
    WARNING() << "THE SET OF POINTS DON'T FIT PARABOLIC WELL, SET y TO y3 ON PURPOSE
    xv = x3;
    yv = y3;
  }
```

  it cannot support `std::numeric_limits<double>::quiet_NaN())`, so I use
  `std::isnan`

- `posix_memalign` doesn't work, use `malloc` instead.
