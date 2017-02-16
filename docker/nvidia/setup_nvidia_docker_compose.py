#!/usr/bin/env python

import argparse
import requests
import yaml

parser = argparse.ArgumentParser(description='Enable GPU-aware dev container.')
parser.add_argument('--nvidia-docker-plugin-url', dest='nvidia_docker_plugin_url', type=str,
                    help='specify what GPU architecture (default: http://localhost:3476)')
parser.add_argument('--gpu-arch', dest='gpu_arch', type=str,
                    help='specify what GPU architecture (default: sm_30)')
parser.add_argument('--cuda-version', dest='cuda_version', type=str,
                    help='specify the CUDA version')
parser.add_argument('--cuda-pkg-version', dest='cuda_pkg_version', type=str,
                    help='specify what CUDA package version to install in the extended image')
args = parser.parse_args()

class nvidia_docker_plugin_api:
    def __init__(self, url):
        self.url = url
        # query the GPU devices status to make sure we are able to reach the API
        self.get('/gpu/status/json')
    def get(self, endpoint):
       try:
           response = requests.get(self.url + endpoint)
       except requests.exceptions.ConnectionError, e:
           print('Cannot connect to the nvidia docker plugin. Did you install it? Is the plugin daemon running on this host?')
           raise e
       return response.json()

if args.nvidia_docker_plugin_url is not None:
    client = nvidia_docker_plugin_api(args.nvidia_docker_plugin_url)
else:
    client = nvidia_docker_plugin_api('http://localhost:3476')

# query nvidia docker plugin for informations about the GPU devices
gpu_devices_info = client.get('/gpu/info/json')

if args.cuda_version is not None:
    cuda_version = args.cuda_version
else:
    cuda_version = gpu_devices_info['Version']['CUDA']

if args.cuda_pkg_version is not None:
    cuda_pkg_version = args.cuda_pkg_version
else:
    # from NVIDIA public hub repository https://hub.docker.com/r/nvidia/cuda/
    cuda_pkg_versions = {
        '8.0': '8-0=8.0.61-1',
        '7.5': '7-5=7.5-18',
        '7.0': '7-0=7.0-28',
    }
    if cuda_version not in cuda_pkg_versions:
        print('No CUDA package version registered for CUDA version \'{0}\''.format(cuda_version))
        print('Available CUDA version are:')
        for k in cuda_pkg_versions.keys(): print(k)
        raise NotImplementedError
    cuda_pkg_version = cuda_pkg_versions[cuda_version]

# I'm not quite sure how to handle the gpu architecture here
# We can extract the info for a given device dy doing
# ```
# gpu_devices_info['Devices'][#]['Arch']
# ```
# It will yield a string (e.g. '3.0') that we might be able to map to the name
# of the class of NVIDIA virtual GPU architecture for which the code must be
# compiled.
if args.gpu_arch is not None:
    gpu_arch = args.gpu_arch
else:
    gpu_arch = 'sm_30'

# query nvidia docker plugin for the command-line parameters to use with the
# `docker run` command
docker_cli_params = client.get('/docker/cli/json')
devices = docker_cli_params['Devices']
volumes = docker_cli_params['Volumes']

# load the template docker compose file to extend the configuration of our
# DTK development container and make it GPU-aware
with open('docker-compose.template.yml', 'r') as fin:
    config = yaml.load(fin)

# add devices and volumes configuration options to the template
config['services']['dtk_dev']['devices'] = devices
config['services']['dtk_dev']['volumes'] = volumes
config['services']['dtk_dev']['environment'] = ['GPU_ARCH=' + gpu_arch]
config['services']['dtk_dev']['build']['args'] = [
    'CUDA_VERSION=' + cuda_version,
    'CUDA_PKG_VERSION=' + cuda_pkg_version,
]
config['volumes'] = {}
config['volumes'][volumes[0].split(':')[0]] = {'external': True}

# write out the extension of the basic DTK docker compose file
with open('docker-compose.yml', 'w') as fout:
    fout.write(yaml.safe_dump(config, default_flow_style=False))
