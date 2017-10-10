#!/bin/bash

SAGE_PATH=../sage  # directory where binary sage is

$SAGE_PATH/sage --preparse lattice_trapdoor_SIS.sage  # create lattice_trapdoor_SIS.py

time $SAGE_PATH/sage graph_based_mmaps.sage  # runs the main program
