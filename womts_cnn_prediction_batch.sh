#! /bin/bash

python_script=/home/weihua/git_repo/oai/direct_super_learn_cnn_v3.py
start_id=0
end_id=10000
step=25
while [ $start_id -lt $end_id ]; do
	step_id=`expr $start_id + $step`
	echo $start_id  $step_id
	python $python_script $start_id $step_id
	start_id=`expr $start_id + $step`
done
