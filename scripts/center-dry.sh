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
if [[ -z "${VMD}" ]]; then
    VMD=/opt/vmd/1.9.4a38/bin/vmd 
    if [[ ! -f $VMD ]]; then
        echo "No vmd found at $VMD"
        exit
    fi
fi
LOG=center-dry.log
PSF=
DCDCSL=
OUTFILE=
SEL="protein or glycan"
POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -vmd)
    VMD="$2"
    shift # past argument
    shift # past value
    ;;
    -log)
    LOG="$2"
    shift # past argument
    shift # past value
    ;;
    -psf)
    PSF="$2"
    shift # past argument
    shift # past value
    ;;
    -dcd)
    DCDCSL="$2"
    shift # past argument
    shift # past value
    ;;
    -o)
    OUTFILE="$2"
    shift
    shift
    ;;
    -opdb)
    OPDB="$2"
    shift
    shift
    ;;
    -sel)
    SEL="$2"
    shift
    shift
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
if [ -z $OUTFILE ]; then
   echo "ERROR: Specify output DCD with -o myout.dcd"
   exit
fi
if [ -z $PSF ]; then
   echo "ERROR: Specify psf file with -psf mypsf.psf"
   exi
fi
set -- "${POSITIONAL[@]}" # restore positional parameters

DCDARGSTR=`echo "$DCDCSL" | sed s/','/' '/g`
cat > tmp.tcl << EOF
mol new $PSF 
foreach dcd \$argv {
   mol addfile \$dcd waitfor all
}
set nf [molinfo top get numframes]
puts "\$nf frames"
set sel [atomselect top "$SEL"]
set selref [atomselect top "$SEL"]
\$selref moveby [vecscale -1.0 [measure center \$selref]]
\$selref frame 0
for { set i 0 } { \$i <= [molinfo top get numframes] } { incr i } {
    \$sel frame \$i
    \$sel move [measure fit \$sel \$selref]
}
animate write dcd $OUTFILE sel \$sel waitfor all 0
EOF
if [ -v OPDB ]; then
cat >> tmp.tcl << EOF
\$sel frame 0
\$sel writepdb $OPDB
EOF
fi
cat >> tmp.tcl << EOF
exit
EOF
$VMD -dispdev text -e tmp.tcl -args $DCDARGSTR > $LOG 2>&1
rm tmp.tcl
exit


