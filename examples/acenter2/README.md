
acenter_eq:
  The analytic centering with equality constraints example of section 9.1 of CVXBook.

To build and run. (Makefile assumes libARMAS is installed on $HOME/.local/lib)

```sh
make all
acenter_data -m 5 -n 6 | cvxsolver -v
```


