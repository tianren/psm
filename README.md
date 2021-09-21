# Multi-party PSM, Revisited

## Improved Communication

Add the following code to the end of `psm.py`
to find a 9-party PSM with communication complexity O(n^4)
```python
general_solver(psm_2blocks_per_party(9), mod=19)
```

## Unbalanced Communication

Add the following code to the end of `psm.py`
to find a 2-party PSM with unbounded communication complexity O(N^{7/20}), O(N^{13/20})
```python
general_solver(psm_2party_tradeoff(20,7), mod=23)
general_solver(psm_2party_tradeoff(20,13), mod=23)
```

## Options

Pass option `verbose=True` to `general_solver` to display the intermediate steps, i.e., the decomposition of all referee-computable masked terms

Pass option `latex=True` to `general_solver` to output in LaTeX format
