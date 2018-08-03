#!/bin/bash
while read line; do kill -SIGKILL $line; done < pids
sleep 1
rm -rf /tmp/ipi*
