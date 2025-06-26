#!/bin/bash
echo "$(date) >>> Workflow started ... "
echo "$(date) >>> Running multi_processor.py ... "
python -u multi_processor.py >multi_processor.log 2>&1

echo "$(date) >>> Finished running multi_processor.py, start extracting and renaming xyz ... "
python -u extract-xyz.py 2>&1

echo "$(date) >>> Finished running extract_xyz.py, start filtering unique SMILES ... "
python -u unismi.py 2>&1

echo "$(date) >>> Finished running extract_xyz.py, start zipping files ... "
cd ./data || exit
zip -rq ligands.zip ligands
zip -rq ligands-xyz.zip ligands-xyz

echo "$(date) >>> Workflow finished."
