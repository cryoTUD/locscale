#!/bin/bash
mapmask mapin $1 mapout $2 << eof
AXIS X Y Z
MODE mapin
eof
