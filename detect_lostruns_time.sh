#!/bin/bash

while [ "1" = "1" ]; do
		echo $(date +"%Y.%m.%d.%S.%N") lpad detect_lostruns --refresh
        lpad detect_lostruns --refresh
        echo "Waiting 3 min"
        sleep 180
done
