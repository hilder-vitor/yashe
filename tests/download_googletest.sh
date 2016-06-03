#!/bin/bash
wget https://github.com/google/googletest/archive/master.zip
unzip master.zip
mv googletest-master/googletest/ .
rm -r googletest-master/
rm master.zip
