wfdbcheck -r sample_mat2wfdb_digital
wfdbcheck -r sample_mat2wfdb
wfdbcheck -r sample_wfdbpython
wfdbcheck -r sample_wrsamp
wfdbcheck -r sample_wrsampNG # -v

rdsamp -r sample_mat2wfdb_digital -c > sample_mat2wfdb_digital_rec.csv
rdsamp -r sample_mat2wfdb -c > sample_mat2wfdb_rec.csv
rdsamp -r sample_wfdbpython -c > sample_wfdbpython_rec.csv
rdsamp -r sample_wrsamp -c > sample_wrsamp_rec.csv
rdsamp -r sample_wrsampNG -c > sample_wrsampNG_rec.csv

# diff sample_mat2wfdb_rec.csv $raw_data_file > raw_mat2wfdb_diff.txt
# diff sample_wfdbpython_rec.csv $raw_data_file > raw_wfdbpython_diff.txt
# diff sample_wrsamp_rec.csv $raw_data_file > raw_wrsamp_diff.txt
# diff sample_mat2wfdb_rec.csv sample_wfdbpython_rec.csv > mat2wfdb_wfdbpython_diff.txt
# diff sample_mat2wfdb_rec.csv sample_wrsamp_rec.csv > mat2wfdb_wrsamp_diff.txt
