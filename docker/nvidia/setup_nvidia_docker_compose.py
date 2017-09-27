#!/usr/bin/env python

import argparse
import os
import requests
import yaml

# Dockerfile template to install cuda toolchain with placeholder in the FROM instruction
dockerfile_template = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Dockerfile_nvidia')
placeholder = '@BASE_IMAGE@'
with open(dockerfile_template, 'r') as fin:
    assert placeholder in fin.readline()

parser = argparse.ArgumentParser(description='Enable GPU-aware dev container.')
parser.add_argument('--nvidia-docker-plugin-url', dest='nvidia_docker_plugin_url', type=str, default='http://localhost:3476',
                    help='nvidia-docker-plugin REST API URL (default: http://localhost:3476)')
parser.add_argument('--gpu-arch', dest='gpu_arch', type=str, default=os.getenv('GPU_ARCH'),
                    help='specify what GPU architecture (default: sm_30)')
parser.add_argument('--cuda-version', dest='cuda_version', type=str, default=os.getenv('CUDA_VERSION'),
                    help='specify the CUDA version')
parser.add_argument('--cuda-pkg-version', dest='cuda_pkg_version', type=str,
                    help='specify what CUDA package version to install in the extended image')
parser.add_argument('--extended-file', dest='compose_file', type=str, default='docker-compose.yml',
                    help='location of the Compose configuration file that defines the services to extend (default: docker-compose.yml)')
parser.add_argument('--extend-services', dest='extend_services', type=str, default='all',
                    help='comma separated list of services to extend (default: all)')
parser.add_argument('--override-file', dest='override_file', type=str, default='docker-compose.override.yml',
                    help='file written that contains configuration override for listed services (default: docker-compose.override.yml)')
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

client = nvidia_docker_plugin_api(args.nvidia_docker_plugin_url)

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
        '9.0': '9-0=9.0.176-1',
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
# Note that `gpu_devices_info['Devices']` is a list and it is unclear whether we
# sould let the user specify what device to use or if we should ensure that all
# have the same architecture.  Therefore for now, I set a default value and let
# the user pass it as an optional argument.  In any case, he `GPU_ARCH`
# environment variable is not added to the image, it is set at the moment the
# container gets launched.  Its value can easily be changed.  It is only used
# in the `docker_cuda_env.sh` script.
if args.gpu_arch is not None:
    gpu_arch = args.gpu_arch
else:
    gpu_arch = 'sm_30'

# query nvidia docker plugin for the command-line parameters to use with the
# `docker run` command
docker_cli_params = client.get('/docker/cli/json')
devices = docker_cli_params['Devices']
volumes = docker_cli_params['Volumes']
assert(len(volumes) == 1)
assert(volumes[0].split(':')[2] == 'ro')
volumes[0] += ',z'

with open(args.compose_file, 'r') as fin:
    base_config = yaml.load(fin)
    compose_file_version = base_config['version']
    if args.extend_services is 'all':
        extend_services = [service for service in base_config['services']]
    else:
        extend_services = args.extend_services.rstrip(',').split(',')
    for service in extend_services:
        assert service in base_config['services']
        # not implemented at this time
        # in that case we would have to append install of CUDA toolchain to the Dockerfile specified in the build section
        assert 'build' not in base_config['services'][service]

# configuration override for the services to make them GPU-aware
config = {}
config['version'] = compose_file_version
config['services'] = {}
for service in extend_services:
    config['services'][service] = {}
    config['services'][service]['devices'] = devices
    config['services'][service]['volumes'] = volumes
    config['services'][service]['environment'] = ['GPU_ARCH=' + gpu_arch]
    config['services'][service]['image'] = service + '_cuda_toolchain'
    dockerfile = 'Dockerfile_' + service + '_cuda_toochain_generated'
    config['services'][service]['build'] = {
        'context': '.',
        'dockerfile': dockerfile,
        'args': [
            'CUDA_VERSION=' + cuda_version,
            'CUDA_PKG_VERSION=' + cuda_pkg_version,
        ],
    }
    with open(dockerfile_template, 'r') as fin:
        with open(dockerfile, 'w') as fout:
            fout.write(fin.read().replace(placeholder, base_config['services'][service]['image']))
config['volumes'] = {}
config['volumes'][volumes[0].split(':')[0]] = {'external': True}

# write out the override of the base configuration
with open(args.override_file, 'w') as fout:
    fout.write(yaml.safe_dump(config, default_flow_style=False))
