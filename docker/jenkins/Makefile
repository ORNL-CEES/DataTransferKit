BASE_FILE=docker-compose.yml
OVERRIDE_FILE=docker-compose.override.yml
SERVICE=ci
SETUP_SCRIPT=../nvidia/setup_nvidia_docker_compose.py
ENV_VARIABLES=CUDA_VERSION GPU_ARCH

all: setup_nvidia_docker

setup_nvidia_docker: ${OVERRIDE_FILE}

${OVERRIDE_FILE}: ${BASE_FILE} ${SETUP_SCRIPT} Makefile ${ENV_VARIABLES}
	${SETUP_SCRIPT} --extended-file ${BASE_FILE} --extend-services ${SERVICE} --override-file ${OVERRIDE_FILE}

clean:
	rm -vf ${OVERRIDE_FILE} ${ENV_VARIABLES}

.PHONY: always-rebuild

.SILENT: ${ENV_VARIABLES}
${ENV_VARIABLES}: always-rebuild
	echo $@=$$$@ > $@.tmp
	diff $@ $@.tmp &>/dev/null || cp $@.tmp $@
	rm -f $@.tmp
