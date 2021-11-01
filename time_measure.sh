#!/bin/bash
start=$(date "+%s")

./euler2D.exe

end=$(date "+%s")
time=$((end-start))
echo "time used:$time seconds"

