
acenter:
  The analytic centering with cone constraints example of section 9.1 of CVXBook.

To build and run. (Makefile assumes libARMAS is installed on $HOME/.local/lib)

```sh
make all
cvxsolver -F acenter -G G.json -h h.json --dims dims.json -v
```


