# Orbiting Dubins Traveling Salesman Problem
# A. Wolek, J. McMahon, Sep-2018
#
# LICENSE:
#
# The source code is in the public domain and not licensed or under
# copyright. The information and software may be used freely by the public.
# As required by 17 U.S.C. 403, third parties producing copyrighted works
# consisting predominantly of the material produced by U.S. government
# agencies must provide notice with such work(s) identifying the U.S.
# Government material incorporated and stating that such material is not
# subject to copyright protection.
#
# Derived works shall not identify themselves in a manner that implies an
# endorsement by or an affiliation with the Naval Research Laboratory.
#
# RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
# SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
# RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
# OF RECIPIENT IN THE USE OF THE SOFTWARE.
#
# (UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.


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
