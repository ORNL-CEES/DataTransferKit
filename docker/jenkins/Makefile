BASE_FILE=docker-compose.yml
OVERRIDE_FILE=docker-compose.override.yml
SERVICE=ci
SETUP_SCRIPT=../nvidia/setup_nvidia_docker_compose.py

all: setup_nvidia_docker

setup_nvidia_docker: ${OVERRIDE_FILE}

${OVERRIDE_FILE}: ${BASE_FILE} ${SETUP_SCRIPT} Makefile
	${SETUP_SCRIPT} --extended-file ${BASE_FILE} --extend-services ${SERVICE} --override-file ${OVERRIDE_FILE}

clean:
	rm -vf *_generated ${OVERRIDE_FILE}
