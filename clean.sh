

# remove from main directory any hidden files and other temporaries 
# (e.g. octave, python files)
rm -rf `find . -name "*~" -print`
#rm -rf *.m
rm -rf octave-workspace
rm -rf *.gtsp *.list *.sol *.tsp *.tour  *.par *.pi *.stats *.cand
rm -rf commands.txt *.cfg
rm -rf ./TMP/*
rm -rf *.m
#rm -rf *.scn
rm -rf cliquer*

for i in "$@"
do
case $i in
    -all)
    # remove the build and binary directories
    rm -rf build
    rm -rf bin
    rm -rf lib
    rm -rf ./TMP/*
    rm LKH
    shift # past argument=value
    ;;
    -ext)
    # remove the build and binary directories
    rm -rf build
    rm -rf bin
    rm -rf lib
    rm -rf ./TMP/*
    rm -rf ./external/armadillo-8.500.0
    rm -rf ./external/ceres-solver  
    rm -rf ./external/cliquer-1.21
    rm -rf ./external/GLKH-1.0
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done
