#!/bin/bash

echo "Setting up your LocScale paths"
CURRENT_DIR=$(pwd)
echo "Use ${CURRENT_DIR} as your default LocScale path? (y/Y)"
read USER_CHOICE
if [[ $USER_CHOICE == "y" || $USER_CHOICE == "Y" ]]; then
	echo "Thank you"
	LOCSCALE_INSTALL=${CURRENT_DIR}
else
	echo "Please enter the path to LocScale (/home/user/path/to/locscale)"
	read USER_PATH_LOCSCALE
	echo "Thank you"
	LOCSCALE_INSTALL=${USER_PATH_LOCSCALE}
fi
echo "Exporting LocScale path to your environment"
echo $LOCSCALE_INSTALL
export LOCSCALE_PATH=$LOCSCALE_INSTALL
echo "Add this path as an environment variable? (y/Y)"
read USER_CHOICE
if [[ $USER_CHOICE == "y" || $USER_CHOICE == "Y" ]]; then
        echo "Adding a new environment variable: "
	echo "LOCSCALE_PATH=${LOCSCALE_PATH}"
	echo "export LOCSCALE_PATH=${LOCSCALE_PATH}" >> ~/.bashrc
else
        echo "Remember to source this file before using LocScale"
        echo "Thank you"
fi
