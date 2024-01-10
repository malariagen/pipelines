
output_dir=/lustre/scratch118/malaria/team112/personal/el10/phase3_data/windowed_accessibility
path_to_release=/lustre/scratch118/malaria/team112/personal/el10/vo_agam_release/v3
path_to_snp_sites=${path_to_release}/snp_genotypes/all/sites
path_to_gambcolu_accessibility=${path_to_release}/site_filters/dt_20200416/gamb_colu
path_to_arab_accessibility=${path_to_release}/site_filters/dt_20200416/arab

mkdir -p $output_dir

source activate cnv37 

python ../scripts/calculate_windowed_mean_accessibility.py $path_to_gambcolu_accessibility $path_to_snp_sites ${output_dir}/mean_accessibility_gambcolu.csv 300 > ${output_dir}/mean_accessibility_gambcolu.log 2>&1
python ../scripts/calculate_windowed_mean_accessibility.py $path_to_arab_accessibility $path_to_snp_sites ${output_dir}/mean_accessibility_arabiensis.csv 300 > ${output_dir}/mean_accessibility_arabensis.log 2>&1

