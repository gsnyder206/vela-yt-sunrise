#/bin/bash

#rsync -av --progress --include="VELA03/" --include="*_sunrise/" --include="images/" --include="broadbandz.fits*" --include="images_*_sunrise.tar" --include="*galprops.npy" --exclude="*" pfe:Runs/VELA_v2/ .
rsync -av --progress --include="VELA31/" --include="VELA32/" --include="VELA33/" --include="VELA34/" --include="VELA35/" --include="*_sunrise/" --include="images/" --include="broadbandz.fits*" --include="images_*_sunrise.tar" --include="*galprops.npy" --exclude="*" pfe:Runs/VELA_v2/ .
