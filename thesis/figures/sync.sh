#!/bin/bash

rsync -av cosma:~/thesis/plots/figures/ ./ --delete --exclude *.sh --exclude *.jpg
