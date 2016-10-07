#! /bin/bash

#Bash-Script to update pyfdap icons

icotool -c -o pyfdap_icon.ico pyfdap_icon.png

echo "Converted .ico"

png2icns pyfdap_icon.icns pyfdap_icon.png

echo "Converted .icns"




