CREST input (ADIT does not run it)

  struct.xyz   starting structure (the conformer with the lowest force-field energy)
  run_crest.sh crest struct.xyz --gfn2 --chrg 0 --uhf 0 -T 1
               charge, unpaired electrons, solvent and threads are copied from the xtb settings in spec.json

Install: conda install -c conda-forge crest   (conda-forge, like xtb)
Results: crest_conformers.xyz (unique conformers), crest_best.xyz (lowest), crest.energies (relative energies). Docs: https://crest-lab.github.io/crest-docs/
