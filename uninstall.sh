#!/bin/bash
# This is a very simplistic uninstall scipt. Use with care!

if [ -n $1 ] && [ "$1" = "--system" ]; then
  rm -v /usr/local/bin/tr-*
  rm -vr /usr/local/lib/python*/site-packages/tracks
else
  rm -v $HOME/bin/tr-*
  rm -vr $HOME/lib/python/tracks
fi
