## A program that calculates the rank of an elliptic curve
Using [Sutherland's library smalljac](https://math.mit.edu/~drew/) we implement a `C` program that calculates the rank of elliptic curves in a parallel manner.

For example we call `result = smalljac_parallel_Lpolys(curve, minp, maxp, flags, dum)` where `curve=[0,1,1,-347,48061]` and `minp=1` and `maxp=100` and (by default) `flags=0` (possible values are SMALLJAC_GOOD_ONLY, SMALLJAC_A1_ONLY, SMALLJAC_DEGREE1_ONLY).

The relevant code is:

```c
start_time = time(0);
result = smalljac_parallel_Lpolys (curve, minp, maxp, flags, dum)
end_time = time(0);
fclose (context.fp);
smalljac_curve_clear(curve);
```
