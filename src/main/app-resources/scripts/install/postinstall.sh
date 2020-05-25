#!/bin/bash

conda update conda -y

cd /application

env_file=$( find /application -name environment.yaml )

env_name=$( cat ${env_file} | head -n 1 | cut -d ':' -f 2 | tr -d ' ' ) 

conda env create --name ${env_name} --file ${env_file}

export PATH=/opt/anaconda/envs/${env_name}/bin:/opt/anaconda/bin:/opt/anaconda/condabin:/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:$PATH

#conda install --name ${env_name} -c terradue -c anaconda -c conda-forge -c defaults cioppy setuptools=41.6.0 pyyaml=3.11=py27_4 gitpython=2.1.3=py27_0 gitdb2=2.0.0=py27_0 nbformat nbconvert -y

#/opt/anaconda/envs/${env_name}/bin/python /application/.util/app_descriptor.py /application /application/application.xml  2>&1

/opt/anaconda/envs/${env_name}/bin/python -m ipykernel install --name ${env_name}

#rm -f /application/util/app_descriptor.py

mkdir -p /workspace/data
