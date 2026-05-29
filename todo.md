### All Todos

- Arguments for input formats
- GenericKey structures
- Ciphertext Structures
- Evaluators add on output ptr -> document
- Revisit justifications:
  - lwe to rgsw uses templates
  - idea is that we want the converter to construct the blindrotator
  - then, converter doesn't need to know all possible br methods
  - and we can plug in whatever we can
  - and let others plug in whatever we want
  - alternative, reframe object+params as constructor+object
  - then pass interface for params, and construct within the clas
  - Possible advantage: compiler opt 
- Fuse blind-rotation algs