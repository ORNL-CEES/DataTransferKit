#!/usr/bin/env python

import requests
import yaml

# query nvidia docker plugin for the command-line parameters to use with the
# `docker run` command
response = requests.get('http://localhost:3476/docker/cli/json')
docker_cli_params = response.json()
devices = docker_cli_params['Devices']
volumes = docker_cli_params['Volumes']

# load the template docker compose file to extend the configuration of our
# DTK development container and make it GPU-aware
with open('docker-compose.template.yml', 'r') as fin:
    config = yaml.load(fin)

# add devices and volumes configuration options to the template
config['services']['dtk_dev']['devices'] = devices
config['services']['dtk_dev']['volumes'] = volumes
config['volumes'] = {}
config['volumes'][volumes[0].split(':')[0]] = {'external': True}

# write out the extension of the basic DTK docker compose file
with open('docker-compose.yml', 'w') as fout:
    fout.write(yaml.safe_dump(config, default_flow_style=False))
