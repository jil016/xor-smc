
## Entry files

### produce.txt:

Example:
```
6 3
0 0 1 1 2 2
```

represents: 6 nodes, 3 productions.
node 0~5 produces 0 0 1 1 2 2 correspondingly.

### demand.txt:

Example:
```
4
2 0 2
3 0 1
4 1 2
5 1 1
```

represents: 4 lines describing demands
`2 0 2`: node 2 needs raw material 0, 2 units.
`3 0 1`: node 3 needs raw material 0, 1 unit.
...


