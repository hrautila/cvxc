
## Convex linear program 

The floorplan example from CVXOPT (CVXbook) on how to use shared object with cvxsolver.

### Building

After building the libcvxc run make on this directory to create the convex program objects.

```sh
$ make
$ ../runsolver.sh -P floorplan.json -v
```

Alternatively.

```sh
$ ../runsolver.sh -F floorplan -G G.json -h h.json -c c.json --dims dims.json -v
```
