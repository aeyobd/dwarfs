#!/bin/bash

rsync -av cosma:~/paper/figures/ ./  --exclude *.sh --exclude *.jpg
