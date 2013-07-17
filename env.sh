#!/bin/bash


function set_environment
{
    case $1 in
        "-r")
        gcc_version="4.7.2"
        source /afs/cern.ch/sw/lcg/external/gcc/$gcc_version/x86_64-slc5/setup.sh ""
        root_version="5.34.09"
        source /afs/cern.ch/sw/lcg/app/releases/ROOT/$root_version/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh
        #echo "GCC $gcc_version and ROOT $root_version are configured"
        ;;
        
        "-p")
        export PEC_FWK_INSTALL=`pwd`
        #echo "Framework-specific variables have been set"
        ;;
    esac
}

if [[ $# -gt 0 && $1 = "-h" ]] ; then
    echo "Usage:"
    echo " . env.sh [options]"
    echo "Options:"
    echo " -r Setup GCC and ROOT"
    echo " -p Setup framework-specific variables"
    echo " -a All the above"
    echo " -h This usage message"
else
    
    if [[ $# -eq 0 ]] || [[ $1 = "-a" ]] ; then
        set_environment -r
        set_environment -p
    else
        for arg in "$@" ; do
            set_environment $arg
        done
    fi
fi
