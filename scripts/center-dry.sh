#!/bin/bash
# 
# use VMD in text mode to center a selection and remove solvent in a DCD trajectory
#
# cameron f abrams cfa22@drexel.edu
#
# arguments [default]:
#    -vmd /path/to/vmd [/opt/vmd/1.9.a38/bin/vmd]
#    -psf myfile.psf REQUIRED
#    -dcd myfile1.dcd[,myfile2.dcd] REQUIRED; comma-separated list of dcd file names
#    -o   outfile.dcd REQUIRED
#    -sel "selection string" ["protein"]
#
if [[ -z "${PSFGEN_BASEDIR}" ]]; then
    PSFGEN_BASEDIR=${HOME}/research/psfgen
    if [[ ! -d $PSFGEN_BASEDIR ]]; then
        echo "Error: No PSFGEN_BASEDIR found."
        exit -1
    fi
fi
source $PSFGEN_BASEDIR/scripts/utils.sh

check_command vmd
check_command catdcd

LOG=center-dry.log
FULLPSF=
DRYPSF=
DCDCSL=
OUTFILE=
REFPDB=
SEL="protein or glycan"
DCDARGSTR=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -log)
    LOG="$2"
    shift # past argument
    shift # past value
    ;;
    -fullpsf)
    FULLPSF="$2"
    shift # past argument
    shift # past value
    ;;
    -drypsf)
    DRYPSF="$2"
    shift # past argument
    shift # past value
    ;;
    -o)
    OUTFILE="$2"
    shift
    shift
    ;;
    -sel)
    SEL="$2"
    shift
    shift
    ;;
    -refpdb)
    REFPDB="$2"
    shift
    shift
    ;;
    *)    # unknown option
    DCDARGSTR+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
if [ -z $OUTFILE ]; then
   echo "ERROR: Specify output DCD with -o myout.dcd"
   exit
fi
if [ -z $FULLPSF ]; then
   echo "ERROR: Specify fully-solvated psf file with -fullpsf mypsf.psf"
   exit
fi
if [ -z $DRYPSF ]; then
   echo "ERROR: Specify pre-solvated (dry) psf file with -drypsf mypsf.psf"
   exit
fi
if [ -z $REFPDB ]; then
   echo "ERROR:  You must specify a reference PDB for all alignments that must be"
   echo "        congruent with the DRYPSF"
   echo " $REFPDB: not found"
   exit
fi
set -- "${DCDARGSTR[@]}" # restore positional parameters

cat > dry.tcl << EOF
mol new $FULLPSF
animate read dcd ${DCDARGSTR[0]} beg 0 end 0
set sel [atomselect top "$SEL"]
set fp [open "index.ndx" "w"]
puts \$fp "[\$sel get index]"
close \$fp
exit
EOF

vmd -dispdev text -e dry.tcl > $LOG 2>&1
if [ $? -ne "0" ]; then
    exit 1
fi
echo "Created index.ndx"

catdcd -i index.ndx -o dry.dcd -otype dcd ${DCDARGSTR[@]}
echo "Created dry.dcd"
cat > center.tcl << EOF
mol new $DRYPSF
mol addfile dry.dcd waitfor all
set nframes [molinfo top get numframes]
set sel [atomselect top "$SEL"]
EOF
if [ ! -z $REFPDB ]; then
cat >> center.tcl << EOF
mol new $DRYPSF
mol addfile $REFPDB
mol top 0
set selref [atomselect 1 "$SEL"]
EOF
else
cat >> center.tcl << EOF
set selref [atomselect top "$SEL"]
EOF
fi
cat >> center.tcl << EOF
\$selref moveby [vecscale -1.0 [measure center \$selref]]
\$selref frame 0
for { set i 0 } { \$i < \$nframes } { incr i } {
    \$sel frame \$i
    \$sel move [measure fit \$sel \$selref]
}
animate write dcd $OUTFILE sel \$sel waitfor all 0
exit
EOF
echo "Centering against $SEL in $REFPDB..."

vmd -dispdev text -e center.tcl >> $LOG 2>&1
if [ $? -ne "0" ]; then
    exit 1
fi

rm -f dry.tcl dry.dcd center.tcl

exit


