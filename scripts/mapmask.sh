#!/bin/bash
mapmask mapin $1 mapout "xyz_"$1 << eof
AXIS X Y Z
MODE mapin
eof

