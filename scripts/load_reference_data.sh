#bin/bash


# data directory is always mounted in the /data


check_exists() {
    if ! [ -d $1 ] ; then
        echo "Error initializing reference data; failed on: $1"
        fail=1
    fi
}

safe_execute() {
    cmd=$1
    echo "running $cmd"
    eval $cmd
    ret_code=$?
    if [ $ret_code != 0 ]; then
        echo $2
        exit $ret_code
    fi
}

fail=0

date

# Move to /data - that's got room for the big files
cd /data

mkdir CHECKM_DATA 
cd CHECKM_DATA

echo "Downloading checkM reference data"
safe_execute "curl --location https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz > checkm_data.tar.gz" "failed to download reference data!"
safe_execute "tar -vxzf checkm_data.tar.gz -C /data/CHECKM_DATA" "failed to untar reference data!"
safe_execute "rm checkm_data.tar.gz" "failed to remove reference data!"
check_exists "/data/CHECKM_DATA/genome_tree" "failed to place reference data"


echo "%%%%%%%%%%%%%%%%%% SETTING /data/CHECKM_DATA/ AS ROOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
checkm data setRoot /data/CHECKM_DATA/
safe_execute "cat /miniconda/lib/python3.6/site-packages/checkm/DATA_CONFIG" "failed to cat DATA_CONFIG"
safe_execute "ls -a /data/CHECKM_DATA" "failed to ls -a /data/CHECKM_DATA"
echo "%%%%%%%%%%%%%%%%%% ~DONE~ SETTING /data/CHECKM_DATA/ AS ROOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%"



if [ $fail -eq 1 ] ; then
    echo "Unable to expand checkM Failing."
    exit 1
fi

echo "Done expanding checkM dataa"
date

echo "Success.  Writing __READY__ file."
touch /data/__READY__
