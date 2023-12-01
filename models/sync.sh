#!/bin/bash

rsync -rv . ~/dwarfs/models --exclude-from 'exclude.txt' --delete
