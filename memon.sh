#! /bin/bash

PID=$1

top -b -p $PID -d 1 | grep $PID
