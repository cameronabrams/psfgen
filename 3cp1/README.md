# 3CP1 -- Trimeric HIV-1 gp41 NHR/CHR six-helix bundle

## About this structure

The six-helix bundle is thought to be the end-state of the fusion cascade for HIV-1 gp41 meant to bring the cell membrane and viral membrane into close proximity to allow for fusion.  This structure was found by NMR on a 1:1 solutions of exised NHR and CHR polypeptides.  The mapping of residue number in the PDB to the reference HXB2 sequence is shown below.

```
  5         15         25         35          53         63         73         83
540        550        560        570         630        640        650        660
TVQ ARQLLSGIVQ QQNDLLRAIE AQQHLLQLTV WGIKQL   ME WDREINNYTS LIHSLIEESQ NQQEKNEQEL LE
---------------NHR-------------------------   -------------CHR----------------------
```
The gap from 577 to 627 is the location of the disulfide loop of gp41, and this is by construction not included in this structure.  For reference, its HXB2 sequence is below.

```
DISULFLIDE LOOP
 580        590        600        610        620
QARI LAVERYLKDQ QLLGIWGCSG KLICTTAVPW NASWSNKSLE QIWHTTW
```

To make a 6HB, we link the NHR and CHR with a GGGGG motif, since we don't know a good structure of the disulfide loop that can be grafted on here in a trimer complex.

Our system building script optionally adds on the membrane-proximal external region (MPER) as an alpha helix, beginning from residue 663:

```
     670        680 683
LDKWASLW NWFNITNWLW YIK
         ---D------ --- 4E10 Ab epitope (4AXW)
A---- 2F5 Ab epitope (3IDI)
```

The two epitopes shown above are for the 2F5 antibody and the 4E10 antibody, identified from theire respective FAB/peptide cocrystal structures (3IDI and 4AXW, respectively).  In a future implementation, the conformation of the MPER that is added on to this 6HB will be taken from these structures, with the aim of generating a complex sterically capable of binding these antibodies.

Our script also optionally adds on the transmembrane domain (after MPER) using the HXB2 sequence:

```
    690        700       709
LFIMIVG GLVGLRIVFA VLSIVNRVR
```
The TMD is grown in as an alpha helix continuing from the MPER.

## Instructions (bash)

```
$ export PSFGEN_BASEDIR=/path/to/psfgen/repository/on/your/system
$ mkdir my_5jyn
$ cd my_5jyn
$ $PSFGEN_BASEDIR/5jyn/do_test.sh [-restart] [-psfgen_args -seed <#> -mper-extend -tmd-extend -do-stalk]
```

Switch | Args | Notes
`-restart` | none | skip everything up to and including packmol and just begin solvated MD
`-seed` | any +int | random number generator seed
`-mper-extend` | none | grow MPER residues as an alpha-helix at C-terminus of CHR's
`-tmd-extend` | none | grow TM residues after MPER's (sets `-mper-extend`)
`-do-stalk` | none | embed in DMPC membrane in putative "stalk" conformation


## References

1. "Impact of the enfuvirtide resistance mutation N43D and the associated baseline polymorphism E137K on peptide sensitivity and six-helix bundle structure." Bai, X., Wilson, K.L., Seedorff, J.E., Ahrens, D., Green, J., Davison, D.K., Jin, L., Stanfield-Oakley, S.A., Mosier, S.M., Melby, T.E., Cammack, N., Wang, Z., Greenberg, M.L., Dwyer, J.J. (2008) Biochemistry 47: 6662-6670.

2019, Cameron F Abrams
