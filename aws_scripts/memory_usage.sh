#! /bin/bash
# monitors memory usage by a process

rm -f mem_usage.log
while :
do
    ps -o rss= $1 >>$1_mem_usage.log
    sleep 0.1s
done
