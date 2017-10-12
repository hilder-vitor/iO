#!/bin/bash

SAGE_PATH=../sage  # directory where binary sage is

$SAGE_PATH/sage --preparse lattice_trapdoor_SIS.sage  # create lattice_trapdoor_SIS.py
$SAGE_PATH/sage --preparse graph_based_mmaps.sage  # create graph_based_mmaps.py

ln -fs lattice_trapdoor_SIS.sage.py lattice_trapdoor_SIS.py
ln -fs graph_based_mmaps.sage.py graph_based_mmaps.py

time $SAGE_PATH/sage multipartite_key_agreement.sage  # runs the main program
