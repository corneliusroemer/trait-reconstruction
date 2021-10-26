# Trait reconstruction with majority rule

## Input

- Newick tree
- Trait for each tip in TSV format

## Output

- Metadata for each internal node in TSV format

## Algorithm

Can start off with majority rule, but hierarchical pango creates problems. How to use hierarchy? Possible traits can be any of the tips, plus any hierarchically superior traits. When to choose superior trait? When there is less than x agreement?
