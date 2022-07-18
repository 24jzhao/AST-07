#!/bin/env bash
mv ast07 ../
mkdir ../ast07/AST-07
mv ./* ../ast07/AST-07
echo "Setup Complete"
rm -rf $(pwd) && exit 1
