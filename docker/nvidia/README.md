# Setting up your DTK development container to leverage NVIDIA GPUs on the host
```bash
[HOST]$ cd docker
[HOST]$ nvidia/setup_nvidia_docker_compose.py
[HOST]$ docker-compose build
[HOST]$ docker-compose up -d
[HOST]$ docker-compose exec dtk_dev bash
[CONTAINER]$ cd $TRILINOS_DIR/DataTransferKit
[CONTAINER]$ source scripts/set_kokkos_env.sh
[CONTAINER]$ # etc.
```

# Requirements
* NVIDIA Docker plugin. See [installation guidelines](https://github.com/NVIDIA/nvidia-docker/wiki/Installation) and consider downloading binary packages from the [release page](https://github.com/NVIDIA/nvidia-docker/releases).
* PyYAML, a YAML parser and emitter for Python.  Consider using your system package manager (e.g. `apt-get install python-yaml` or `yum install python-yaml`) or pip (`pip install pyyaml`).
