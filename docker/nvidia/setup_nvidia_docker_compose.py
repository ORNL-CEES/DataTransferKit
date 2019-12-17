#!/usr/bin/env python

import argparse
import os
import yaml

parser = argparse.ArgumentParser(description='Enable GPU-aware dev container.')
parser.add_argument('--gpu-arch', dest='gpu_arch', type=str, default=os.getenv('GPU_ARCH', 'sm_70'),
                    help='specify what GPU architecture (default: sm_70)')
parser.add_argument('--cuda-version', dest='cuda_version', type=str, default=os.getenv('CUDA_VERSION', '9.0'),
        help='specify the CUDA version (default: 9.0)')
parser.add_argument('--extended-file', dest='compose_file', type=str, default='docker-compose.yml',
                    help='location of the Compose configuration file that defines the services to extend (default: docker-compose.yml)')
parser.add_argument('--extend-services', dest='extend_services', type=str, default='all',
                    help='comma separated list of services to extend (default: all)')
parser.add_argument('--override-file', dest='override_file', type=str, default='docker-compose.override.yml',
                    help='file written that contains configuration override for listed services (default: docker-compose.override.yml)')
args = parser.parse_args()

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

# see https://github.com/NVIDIA/nvidia-docker/wiki/Frequently-Asked-Questions#do-you-support-docker-compose
# NOTE consider setting nvidia runtime as the default in which case I don't think the version of docker-compose will matter any more
assert compose_file_version == '2.3'

# configuration override for the services to make them GPU-aware
config = {}
config['version'] = compose_file_version
config['services'] = {}
for service in extend_services:
    config['services'][service] = {}
    config['services'][service]['environment'] = ['GPU_ARCH=' + args.gpu_arch]
    config['services'][service]['image'] = base_config['services'][service]['image'] + '-cuda' + args.cuda_version.replace('.', '')
    config['services'][service]['runtime'] = 'nvidia'

# write out the override of the base configuration
with open(args.override_file, 'w') as fout:
    fout.write(yaml.safe_dump(config, default_flow_style=False))
