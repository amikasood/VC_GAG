for i in *.h *.cpp; do
	#diff ${i} /wilde_spaces/scratch1/yao/temp_vc_debug/VC_GAG/autodock_vina_1_1_2/src/lib/${i} > all_diff/diff_${i}.txt
	diff ${i} /home/yao/temp_vc_offcial/autodock_vina_1_1_2/src/lib/${i} > all_diff/diff_${i}.txt
done
