include .env
BASE_FILE=docker-compose.yml
OVERRIDE_FILE=docker-compose.override.yml
SERVICE=ci
SETUP_SCRIPT=../nvidia/setup_nvidia_docker_compose.py
export COMPOSE_FILE=${BASE_FILE}:${OVERRIDE_FILE}

all:
	make run
config:
	docker-compose config
build: ${OVERRIDE_FILE}
	docker-compose build --pull ${SERVICE}
${OVERRIDE_FILE}: ${BASE_FILE}
	${SETUP_SCRIPT} --extended-file ${BASE_FILE} --extend-services ${SERVICE} --override-file ${OVERRIDE_FILE}
run: build
	docker-compose run ${SERVICE}
down:
	docker-compose down
clean:
	rm -vf *_generated ${OVERRIDE_FILE}
