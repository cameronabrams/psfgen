# 5JYN -- Trimeric HIV-1 gp41 transmembrane domain in DMPC bilayer

## About this structure

The transmembrane domain (TM) of gp41 anchors the trimeric envelope
spike complex to the viral membrane of HIV-1.  This protein complex 
structure was solved using NMR spectroscopy on TM peptides 
reconstituted in DMPC bicelles [1].  In this PSF generation directory, 
I demonstrate how to use packmol in a staged approach to build an 
equilibrated protein-bilayer system.  In particular, I use a 
model-built template for a single DMPC molecule.  Note that this
does not involve using the CHARMM-GUI membrane builder, since
this is such a simple membrane.

## Instructions (bash)

```
$ export PSFGEN_BASEDIR=/path/to/psfgen/repository/on/your/system
$ mkdir my_5jyn
$ cd my_5jyn
$ $PSFGEN_BASEDIR/5jyn/do_test.sh [-restart] [-psfgen_args -seed <#>]
```

The `-restart` switch, if present, tells the script to skip everything up to and including packmol and just begin solvated MD.  One can also set the random number seed (used by packmol) with the `-seed` switch.

## References

1. "Structural basis for membrane anchoring of HIV-1 envelope spike." Dev, J., Park, D., Fu, Q., Chen, J., Ha, H.J., Ghantous, F., Herrmann, T., Chang, W., Liu, Z., Frey, G., Seaman, M.S., Chen, B., Chou, J.J. Science 2016;353:172-175 

2018, Cameron F Abrams
